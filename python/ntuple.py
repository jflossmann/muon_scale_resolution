import ROOT
from array import array
import yaml
import time
from multiprocessing import Pool, RLock
from tqdm import tqdm
import os
import json

# function in c++ code. Finds indices of muons closest to z boson mass
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<Int_t> get_indices(
        UInt_t nMuon,
        ROOT::VecOps::RVec<Float_t> *Muon_pt,
        ROOT::VecOps::RVec<Float_t> *Muon_eta,
        ROOT::VecOps::RVec<Float_t> *Muon_phi,
        ROOT::VecOps::RVec<Float_t> *Muon_mass,
        ROOT::VecOps::RVec<UChar_t> *Muon_pfIsoId,
        ROOT::VecOps::RVec<Bool_t> *Muon_mediumId,
        ROOT::VecOps::RVec<Int_t> *Muon_charge
        ){
        Float_t deltaM = 1000;
        Int_t ind1, ind2;
        ind1 = -99;
        ind2 = -99;
        for(int i=0; i<nMuon; i++){
            if (Muon_pt->at(i) < 10) continue;
            if(fabs(Muon_eta->at(i) > 2.4)) continue;
            if(Muon_pfIsoId->at(i) < 2) continue;
            if(Muon_mediumId->at(i) == 0) continue;

            for(int j=i; j<nMuon; j++){
                if (Muon_pt->at(j) < 10) continue;
                if (Muon_pt->at(i) < 25 && Muon_pt->at(j) < 25) continue;
                if(fabs(Muon_eta->at(j) > 2.4)) continue;
                if(Muon_pfIsoId->at(j) < 2) continue;
                if(Muon_mediumId->at(j) == 0) continue;
                if(Muon_charge->at(i) * Muon_charge->at(j) > 0) continue;

                TLorentzVector mui, muj, Z;
                mui.SetPtEtaPhiM(
                    Muon_pt->at(i),
                    Muon_eta->at(i),
                    Muon_phi->at(i),
                    Muon_mass->at(i)
                );
                muj.SetPtEtaPhiM(
                    Muon_pt->at(j),
                    Muon_eta->at(j),
                    Muon_phi->at(j),
                    Muon_mass->at(j)
                );
                Z = mui + muj;
                if (fabs(Z.M() - 91.1876) < deltaM){
                    deltaM = fabs(Z.M() - 91.1876);
                    if (Muon_charge->at(i) < 0){
                        ind1 = i;
                        ind2 = j;
                    }
                    else {
                        ind1 = j;
                        ind2 = i;
                    }
                }
            }
        }

        ROOT::VecOps::RVec<Int_t> s(2);
        s[0] = ind1;
        s[1] = ind2;
        return s;
    }
""")



# function for trigger matching
ROOT.gInterpreter.Declare("""
    UInt_t trg_match_ind(
        Float_t eta,
        Float_t phi,
        Int_t nTrigObj,
        ROOT::VecOps::RVec<UShort_t> *TrigObj_id,
        ROOT::VecOps::RVec<Float_t> *TrigObj_eta,
        ROOT::VecOps::RVec<Float_t> *TrigObj_phi,
        Int_t match1
    ){
        Int_t index = -99;
        Float_t dRmin = 1000;
        Float_t dR, dEta, dPhi;
        for(int i=0; i<nTrigObj; i++){
            if (TrigObj_id->at(i) != 13) continue;
            if (TrigObj_id->at(i) == match1) continue;
            dEta = eta - TrigObj_eta->at(i);
            dPhi = (phi - TrigObj_phi->at(i));
            if (dPhi > 3.1415) dPhi = 2*3.1415 - dPhi;
            dR = sqrt(dEta*dEta + dPhi*dPhi);
            if (dR > 0.1) continue;
            if (index > -1){
                if(dR < dRmin) {
                    dRmin = dR;
                    index = i;
                }
                else continue;
            }
            else index = i;
        }
        return index;
    }
""")

# function in c++ code. Finds indices of muons matched to GEN
ROOT.gInterpreter.Declare(
"""
    Int_t muon_genmatch(
        Float_t eta,
        Float_t phi,
        Int_t charge,
        ROOT::VecOps::RVec<Int_t> *GenPart_status,
        ROOT::VecOps::RVec<UShort_t> *GenPart_statusFlags,
        ROOT::VecOps::RVec<Int_t> *GenPart_pdgId,
        ROOT::VecOps::RVec<Short_t> *GenPart_genPartIdxMother,
        ROOT::VecOps::RVec<Float_t> *GenPart_eta,
        ROOT::VecOps::RVec<Float_t> *GenPart_phi
    ){
        Int_t index=-99;
        
        Float_t deltaR=0.1;
        Float_t dEta, dPhi, dR;
        Int_t mother_idx, mother_id;

        for(int j=0; j<GenPart_eta->size(); j++){
            if (fabs(GenPart_pdgId->at(j)) != 13) continue;
            if (GenPart_pdgId->at(j) * charge > 0) continue;
            if (GenPart_status->at(j) != 1) continue;
            if ((GenPart_statusFlags->at(j) >> 13) % 2 == 0) continue;

            mother_idx = GenPart_genPartIdxMother->at(j);
            mother_id = GenPart_pdgId->at(mother_idx);
            
            while (mother_id != 23){
                mother_idx = GenPart_genPartIdxMother->at(mother_idx);
                if (mother_idx < 0) break;
                mother_id = GenPart_pdgId->at(mother_idx);
            }
            if (mother_id != 23) continue;

            dEta = eta - GenPart_eta->at(j);
            dPhi = abs(phi - GenPart_phi->at(j));
            if (dPhi > 3.1415) dPhi = 2*3.1415 - dPhi;
            dR = sqrt(dEta*dEta + dPhi*dPhi);
            
            if(dR < deltaR){
                deltaR = dR;
                index = j;
            }
        }
    return index;
    }
"""
)



ROOT.gInterpreter.Declare(
"""
    Int_t is_golden(
        Int_t luminosityBlock,
        Int_t run,
        std::vector<std::vector<std::vector<int>>> glumisecs,
        std::vector<int> runs
    ){
        int index = -1;
        auto runidx = std::find(runs.begin(), runs.end(), run);
        if (runidx != runs.end()){
            index = runidx - runs.begin();
        }

        if (index == -1) {
            return 0;
        }

        int n_glumsecs = glumisecs[index].size();
        for (int i=0; i<n_glumsecs; i++) {
            if (luminosityBlock >= glumisecs[index][i][0] && luminosityBlock <= glumisecs[index][i][1]){
                return 1;
            }
        }
        return 0;
    }
"""
)

def filter_lumi(rdf, golden_json):
    """
    function to get rdf with events filtered with golden json
    
    (ROOT.RDataFrame) rdf: dataframe
    (str) golden_json: path to golden json
    """

    # load content of golden json
    with open(golden_json) as cf:
        goldenruns = json.load(cf)

    # extract runs and lumi sections to lists
    runs = [r for r in goldenruns.keys()]
    lumlist = [goldenruns[r] for r in goldenruns.keys()]

    # make c++ vectors of runlist and lumilist for dataframe
    runstr = "{" + ",".join(runs) + "}"
    lumstr = str(lumlist).replace("[", "{").replace("]", "}")

# function to make z_pt distributions 
def hist_zpt(ntuples, pt_bins, hdir):
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)
    hists = []

    for typ in ntuples:
        for sample in ntuples[typ]:
            rdf = ROOT.RDataFrame("Events", ntuples[typ][sample])
            print(ntuples[typ][sample])

            h = rdf.Histo1D((f"h_Zboson_pt_{sample}", "", len(pt_bins)-1, array('f', pt_bins)), "pt_Z", 'genWeight')
            print(h.Integral())
            h.Scale(1./h.Integral())
            hists.append(h)

    tf = ROOT.TFile(f"{hdir}z_reweighting.root", "RECREATE")
    for h in hists:
        h.Write()
    tf.Close()

    return 


ROOT.gInterpreter.Declare("""
        float gaus(){
            return gRandom->Gaus(0,1);
        }
""")



# function to add z pt weight. TODO: also for GEN or use reco instead?
def weight_zpt(ntuples, hdir, eta_bins, phi_bins, SFs):
    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}z_reweighting.root", "READ");')
    ROOT.gROOT.ProcessLine('TH1D* h_dt = (TH1D*)tf->Get("h_Zboson_pt_DATA");')
    ROOT.gROOT.ProcessLine('TH1D* h_mc = (TH1D*)tf->Get("h_Zboson_pt_SIG");')
    ROOT.gROOT.ProcessLine('TH1D* h_ratio = (TH1D*)h_dt->Clone();')
    ROOT.gROOT.ProcessLine('h_ratio->Divide(h_mc)')
    trgpath = SFs['trg'][0]
    trgname = SFs['trg'][1]
    ROOT.gROOT.ProcessLine(f'TFile *trg_tf = TFile::Open("{trgpath}", "READ");')
    ROOT.gROOT.ProcessLine(f'TH2F *trg_mc = (TH2F*)trg_tf->Get("{trgname}_efficiencyMC");')
    ROOT.gROOT.ProcessLine(f'TH2F *trg_dt = (TH2F*)trg_tf->Get("{trgname}_efficiencyData");')


    for sf in SFs:
        if sf=='trg': continue
        with open(SFs[sf][0]) as _file:
            tmp_dicts = json.load(_file)["corrections"]
            for tmp_dict in tmp_dicts:
                if tmp_dict["name"] == SFs[sf][1]:
                    SFs[sf].append(tmp_dict)

    for typ in ntuples:
        for sample in ntuples[typ]:
            #if sample !='SIG' : continue
            print(sample)
            print(ntuples[typ][sample].replace('_*.root', '.yaml').replace(sample, f'ntuples/{sample}'))
            sampleyaml = yaml_loader(ntuples[typ][sample].replace('_*.root', '.yaml'))

            rdf = ROOT.RDataFrame("Events", ntuples[typ][sample])
            rdf = rdf.Define('sumwWeight', str(sampleyaml['genweight']))
            rdf = rdf.Define('xsec', str(sampleyaml['xsec']))
            if typ != "SIG" and typ != "GEN":
                rdf = rdf.Define("zPtWeight", "1")
            else:
                rdf = rdf.Define("zPtWeight", "h_ratio->GetBinContent(h_ratio->FindBin(pt_Z))")
        
            for sf in SFs:
                rdf = rdf.Define(f"sf_{sf}_1", "1")
                rdf = rdf.Define(f"sf_{sf}_2", "1")

            if typ != 'DATA':
                for sf in SFs:
                    if sf != 'trg':
                        etabins = SFs[sf][2]["data"]["edges"]
                        if SFs[sf][2]["inputs"][0]["name"] != "abseta":
                            print("The sf does not seem to be binned in abseta. Please adjust.")
                            break

                        for eta in range(len(etabins)-1):
                            tmp_dict = SFs[sf][2]["data"]["content"][eta]
                            ptbins = tmp_dict["edges"]

                            for pt in range(len(ptbins)-1):

                                systs = tmp_dict["content"][pt]["content"]
                                for syst in systs:
                                    if syst["key"] == 'nominal':
                                        sf_val = syst["value"]

                                rdf = rdf.Redefine(
                                    f"sf_{sf}_1",
                                    f'double sf;\
                                    if (abs(eta_1) > {etabins[eta]} && abs(eta_1) < {etabins[eta+1]} &&\
                                    pt_1 > {ptbins[pt]} && pt_1 < {ptbins[pt+1]}) sf = {sf_val};\
                                    else sf = sf_{sf}_1;\
                                    return sf;'
                                )
                                rdf = rdf.Redefine(
                                    f"sf_{sf}_2",
                                    f'double sf;\
                                    if (abs(eta_2) > {etabins[eta]} && abs(eta_2) < {etabins[eta+1]} &&\
                                    pt_2 > {ptbins[pt]} && pt_2 < {ptbins[pt+1]}) sf = {sf_val};\
                                    else sf = sf_{sf}_2;\
                                    return sf;'
                                )
                    else:
                        rdf = rdf.Redefine(
                            "sf_trg_1",
                            'trg_dt->GetBinContent(trg_dt->FindBin(abs(eta_1), pt_1)) /\
                             trg_mc->GetBinContent(trg_mc->FindBin(abs(eta_1), pt_1))'
                        )
                        rdf = rdf.Redefine(
                            "sf_trg_2",
                            'trg_dt->GetBinContent(trg_dt->FindBin(abs(eta_2), pt_2)) /\
                             trg_mc->GetBinContent(trg_mc->FindBin(abs(eta_2), pt_2))'
                        )


            for sf in SFs:
                if sf != 'trg' or typ=='DATA':
                    rdf = rdf.Define(f"sf_{sf}", f"sf_{sf}_1 * sf_{sf}_2")
                else:
                    rdf = rdf.Define(
                        "sf_trg",
                        '''
                        double sf = 1.0;  // Default scale factor
                        double eff_mc_1 = trg_mc->GetBinContent(trg_mc->FindBin(abs(eta_1), pt_1));
                        double eff_mc_2 = trg_mc->GetBinContent(trg_mc->FindBin(abs(eta_2), pt_2));
                        double eff_dt_1 = trg_dt->GetBinContent(trg_dt->FindBin(abs(eta_1), pt_1));
                        double eff_dt_2 = trg_dt->GetBinContent(trg_dt->FindBin(abs(eta_2), pt_2));

                        if(trg_match_1 == 0){
                            eff_mc_1 = 0;
                            eff_dt_1 = 0;
                        }
                        if(trg_match_2 == 0){
                            eff_mc_2 = 0;
                            eff_dt_2 = 0;
                        }

                        double eff_mc = 1.0 - (1.0 - eff_mc_1) * (1.0 - eff_mc_2);
                        double eff_dt = 1.0 - (1.0 - eff_dt_1) * (1.0 - eff_dt_2);

                        if (eff_mc == 0) {
                            std::cout << "mc efficiency is 0 in this bin: " << pt_1 << ", " << pt_2 << ", " << eta_1 << ", " << eta_2 << std::endl;
                            return 1.0;  // Return default scale factor if eff_mc is zero
                        } else {
                            sf = eff_dt / eff_mc;  // Calculate scale factor
                            return sf;
                        }
                        '''
                    )

            rdf.Snapshot("Events", ntuples[typ][sample].replace("*.root", "zPt.root").replace('ntuples/',''))

    return



# function which creates ntuple files from nanoaod
def make_ntuples(nanoAODs, datasets, datadir, golden_json):

    os.makedirs(datadir, exist_ok=True)

    nanoAODs = yaml_loader(nanoAODs)
    datasets = yaml_loader(datasets)

    for sample in nanoAODs:
        # if sample != 'DY': continue
        sum_genweights = 0

        start = time.time()
        print(f"Processing {sample} samples. Number of Files: {len(nanoAODs[sample])}")
        args = [(datasets, sample, _file, number, datadir, golden_json) for number, _file in enumerate(nanoAODs[sample])]
        nthreads = 32

        pool = Pool(nthreads, initargs=(RLock,), initializer=tqdm.set_lock)
        for genweight in tqdm(
            pool.imap_unordered(job_wrapper, args),
            total=len(args),
            desc="Total progess",
            dynamic_ncols=True,
            leave=True
        ): sum_genweights += genweight

        datasets[sample]['genweight'] = sum_genweights

        with open(datadir + sample + '.yaml', "w") as f:
            yaml.dump(datasets[sample], f)
        if sample == 'DY':
            with open(datadir + 'GEN.yaml', 'w') as f:
                yaml.dump(datasets[sample], f)

        end = time.time()
        print(f"Finished processing of {sample} samples in {round(end-start,1)} s.")


def job_wrapper(args):
    return process_ntuples(*args)


def process_ntuples(datasets, sample, _file, number, datadir, golden_json):
    quants = [
            "pt_Z", "mass_Z", "eta_Z", "phi_Z",
            "pt_1", "mass_1", "eta_1", "phi_1", "charge_1",
            "pt_2", "mass_2", "eta_2", "phi_2", "charge_2",
            "genWeight",
            "nTrkLayers_1", "nTrkLayers_2",
            "trg_match_1", "trg_match_2"
        ]
    # load nanoAOD
    rdf = ROOT.RDataFrame("Events", _file)

    if sample == 'DATA':
        rdf = rdf.Define("genWeight", "1")

    genweight = rdf.Sum("genWeight").GetValue()
    
    # only collect events w/ >1 muon and find muon pair closest to z mass. Muon1 is always charge -1 and muon2 always +1
    rdf = rdf.Filter("HLT_IsoMu24 == 1")
    rdf = rdf.Filter("Muon_pt.size() > 1")
    rdf = rdf.Define("ind", """ROOT::VecOps::RVec<Int_t> (get_indices(
                                nMuon,
                                &Muon_pt,
                                &Muon_eta,
                                &Muon_phi,
                                &Muon_mass,
                                &Muon_pfIsoId,
                                &Muon_mediumId,
                                &Muon_charge
                                ))""")
    rdf = rdf.Define("ind0", "ind[0]")
    rdf = rdf.Define("ind1", "ind[1]")
    rdf = rdf.Filter("ind0 + ind1 > 0")
    rdf = rdf.Define("pt_1", "Muon_pt[ind[0]]")
    rdf = rdf.Define("pt_2", "Muon_pt[ind[1]]")
    rdf = rdf.Define("mass_1", "Muon_mass[ind[0]]")
    rdf = rdf.Define("mass_2", "Muon_mass[ind[1]]")
    rdf = rdf.Define("eta_1", "Muon_eta[ind[0]]")
    rdf = rdf.Define("eta_2", "Muon_eta[ind[1]]")
    rdf = rdf.Define("phi_1", "Muon_phi[ind[0]]")
    rdf = rdf.Define("phi_2", "Muon_phi[ind[1]]")
    rdf = rdf.Define("charge_1", "Muon_charge[ind[0]]")
    rdf = rdf.Define("charge_2", "Muon_charge[ind[1]]")
    rdf = rdf.Define("nTrkLayers_1", "Muon_nTrackerLayers[ind[0]]")
    rdf = rdf.Define("nTrkLayers_2", "Muon_nTrackerLayers[ind[1]]")

    # check trigger matches
    rdf = rdf.Define("trg_ind0", "trg_match_ind(eta_1, phi_1, nTrigObj, &TrigObj_id, &TrigObj_eta, &TrigObj_phi, -99)")
    rdf = rdf.Define("trg_ind1", "trg_match_ind(eta_2, phi_2, nTrigObj, &TrigObj_id, &TrigObj_eta, &TrigObj_phi, trg_ind0)")
    rdf = rdf.Filter("trg_ind0 >= 0 || trg_ind1 >= 0")

    rdf = rdf.Define("trg_match_1", "int match; if(trg_ind0 >= 0) match = 1; else match = 0; return match;")
    rdf = rdf.Define("trg_match_2", "int match; if(trg_ind1 >= 0) match = 1; else match = 0; return match;")

    # Define quantities of Z boson and collect those events with 50 < m_Z < 130
    rdf = rdf.Define("p4_1", "ROOT::Math::PtEtaPhiMVector(pt_1, eta_1, phi_1, mass_1)")
    rdf = rdf.Define("p4_2", "ROOT::Math::PtEtaPhiMVector(pt_2, eta_2, phi_2, mass_2)")
    rdf = rdf.Define("p4_Z", "p4_1 + p4_2")
    rdf = rdf.Define("pt_Z", "p4_Z.Pt()")
    rdf = rdf.Define("eta_Z", "p4_Z.Eta()")
    rdf = rdf.Define("mass_Z", "p4_Z.M()")
    rdf = rdf.Define("phi_Z", "p4_Z.Phi()")
    
    rdf = rdf.Filter("mass_Z > 50 && mass_Z < 130")

    if sample=='DATA':
        # load content of golden json
        with open(golden_json) as cf:
            goldenruns = json.load(cf)

        # extract runs and lumi sections to lists
        runs = [r for r in goldenruns.keys()]
        lumlist = [goldenruns[r] for r in goldenruns.keys()]

        # make c++ vectors of runlist and lumilist for dataframe
        runstr = "{" + ",".join(runs) + "}"
        lumstr = str(lumlist).replace("[", "{").replace("]", "}")
        # print(runstr, lumstr)
        rdf = rdf.Define(
            "golden", 
            f"Int_t (is_golden(luminosityBlock, run, {lumstr}, {runstr}))"
        )
        n_all = rdf.Count().GetValue()
        rdf = rdf.Filter("golden > 0")
        n_golden = rdf.Count().GetValue()
        print(f"Total number of events: {n_all}, golden events: {n_golden}.")

    # make output with interesting data
    rdf.Snapshot("Events", f"{datadir+sample}_{number}.root", quants)
    
    if sample=="DY":
        start = time.time()
        sample = 'GEN'
        print(f"Calculation of Gen quantities for GEN samples.")
        # perform gen delta R matching and collect corresponding events and gen quantities
        rdf = rdf.Define("genind_1", """muon_genmatch(
                                            eta_1,
                                            phi_1,
                                            charge_1,
                                            &GenPart_status,
                                            &GenPart_statusFlags,
                                            &GenPart_pdgId,
                                            &GenPart_genPartIdxMother,
                                            &GenPart_eta,
                                            &GenPart_phi
                                            )""")
        rdf = rdf.Define("genind_2", """muon_genmatch(
                                            eta_2,
                                            phi_2,
                                            charge_2,
                                            &GenPart_status,
                                            &GenPart_statusFlags,
                                            &GenPart_pdgId,
                                            &GenPart_genPartIdxMother,
                                            &GenPart_eta,
                                            &GenPart_phi
                                            )""")
        # if gen matching successfull: save information of gen event
        rdf = rdf.Define("gen_mask", "genind_1 != -99 && genind_2 != -99 && genind_1 != genind_2")
        rdf = rdf.Filter("gen_mask")

        rdf = rdf.Define("genpt_1", "GenPart_pt[genind_1]")
        rdf = rdf.Define("genpt_2", "GenPart_pt[genind_2]")
        rdf = rdf.Define("geneta_1", "GenPart_eta[genind_1]")
        rdf = rdf.Define("geneta_2", "GenPart_eta[genind_2]")
        rdf = rdf.Define("genphi_1", "GenPart_phi[genind_1]")
        rdf = rdf.Define("genphi_2", "GenPart_phi[genind_2]")
        rdf = rdf.Define("genmass_1", "GenPart_mass[genind_1]")
        rdf = rdf.Define("genmass_2", "GenPart_mass[genind_2]")
        rdf = rdf.Define("gencharge_1", "- GenPart_pdgId[genind_1]/fabs(GenPart_pdgId[genind_1])")
        rdf = rdf.Define("gencharge_2", "- GenPart_pdgId[genind_2]/fabs(GenPart_pdgId[genind_2])")

        rdf = rdf.Define("genp4_1", "ROOT::Math::PtEtaPhiMVector(genpt_1, geneta_1, genphi_1, genmass_1)")
        rdf = rdf.Define("genp4_2", "ROOT::Math::PtEtaPhiMVector(genpt_2, geneta_2, genphi_2, genmass_2)")
        rdf = rdf.Define("genp4_Z", "genp4_1 + genp4_2")
        rdf = rdf.Define("genpt_Z", "genp4_Z.Pt()")
        rdf = rdf.Define("geneta_Z", "genp4_Z.Eta()")
        rdf = rdf.Define("genmass_Z", "genp4_Z.M()")
        rdf = rdf.Define("genphi_Z", "genp4_Z.Phi()")

        
        quants += [
            "genpt_1", "genmass_1", "geneta_1", "genphi_1", "gencharge_1",
            "genpt_2", "genmass_2", "geneta_2", "genphi_2", "gencharge_2",
            "genpt_Z", "genmass_Z", "geneta_Z", "genphi_Z"
        ]

        # make output with interesting data
        rdf.Snapshot("Events", f"{datadir+sample}_{number}.root", quants)
        end = time.time()
        print(f"Finished processing of GEN samples in {round(end-start,1)} s.")

    return genweight


def yaml_loader(fname):
    with open(fname, "r") as f:
        dsets = yaml.load(f, Loader=yaml.Loader)
    #print(dsets)
    return dsets



ROOT.gInterpreter.Declare("""
        int poisson(){
            return gRandom->Poisson(1);
        }
""")
