import ROOT
from array import array

# function in c++ code. Finds indices of muons closest to z boson mass
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<Int_t> get_indices(
        UInt_t nMuon,
        ROOT::VecOps::RVec<Float_t> *Muon_pt,
        ROOT::VecOps::RVec<Float_t> *Muon_eta,
        ROOT::VecOps::RVec<Float_t> *Muon_phi,
        ROOT::VecOps::RVec<Float_t> *Muon_mass,
        ROOT::VecOps::RVec<Float_t> *Muon_tkRelIso,
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
            if(Muon_tkRelIso->at(i) > 0.1) continue;
            if(Muon_mediumId->at(i)==0) continue;

            for(int j=i; j<nMuon; j++){
                if (Muon_pt->at(j) < 10) continue;
                if(fabs(Muon_eta->at(j) > 2.4)) continue;
                if(Muon_tkRelIso->at(j) > 0.1) continue;
                if(Muon_mediumId->at(j)==0) continue;
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

# function in c++ code. Finds indices of muons matched to HLT; TODO look up definition of hlt
#ROOT.gInterpreter.Declare(
"""
    ROOT::VecOps::RVec<Int_t> get_hlt_matches(
        UInt_t nMuon,
        ROOT::VecOps::RVec<Float_t> *Muon_pt,
        ROOT::VecOps::RVec<Float_t> *Muon_eta,
        ROOT::VecOps::RVec<Float_t> *Muon_phi,
        ROOT::VecOps::RVec<Float_t> *hlt_eta,
        ROOT::VecOps::RVec<Float_t> *hlt_phi,
        ROOT::VecOps::RVec<Float_t> *Muon_tkRelIso,
        ROOT::VecOps::RVec<Float_t> *Muon_mediumPromptId
        ){
        ROOT::VecOps::RVec<Int_t> muon_hlt_match;
        for(int i=0; i<nMuon; i++){
            if(muon_hlt_match.size()>1) continue;
            Float_t deltaR=1000;
            if (Muon_pt->at(i) < 10) continue;
            if(fabs(Muon_eta->at(i) > 2.4)) continue;
            if(Muon_tkRelIso->at(i) > 0.1) continue;
            if(Muon_mediumPromptId->at(i)) continue;
            for(int j=0; j<hlt_eta.size(); j++){
                dEta = Muon_eta->at(i) - hlt_eta->at(j);
                dPhi = Muon_phi->at(i) - hlt_phi->at(j);
                deltaR = sqrt(dEta*dEta + dPhi*dPhi);
                if(deltaR < 0.35) muon_hlt_match.push_back(i);
            }
        }
    return muon_hlt_matches;
    }
"""
#)


# function in c++ code. Finds indices of muons matched to GEN
ROOT.gInterpreter.Declare(
"""
    Int_t muon_genmatch(
        Float_t eta,
        Float_t phi,
        Int_t charge,
        ROOT::VecOps::RVec<Int_t> *GenPart_status,
        ROOT::VecOps::RVec<Int_t> *GenPart_pdgId,
        ROOT::VecOps::RVec<Int_t> *GenPart_genPartIdxMother,
        ROOT::VecOps::RVec<Float_t> *GenPart_eta,
        ROOT::VecOps::RVec<Float_t> *GenPart_phi
    ){
        Int_t index=-99;
        
        Float_t deltaR=0.1;
        Float_t dEta, dPhi, dR;
        
        for(int j=0; j<GenPart_eta->size(); j++){
            if (GenPart_status->at(j) != 1) continue;
            if (fabs(GenPart_pdgId->at(j)) != 13) continue;
            if (GenPart_pdgId->at(GenPart_genPartIdxMother->at(j)) != 23) continue;
            if (GenPart_pdgId->at(j) * charge > 0) continue;

            dEta = eta - GenPart_eta->at(j);
            dPhi = phi - GenPart_phi->at(j);
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


# function to make z_pt distributions 
def hist_zpt(ntuples, pt_bins, hists):
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)
    hists = []

    for s in ntuples:

        rdf = ROOT.RDataFrame("Events", ntuples[s])

        h = rdf.Histo1D((f"h_Zboson_pt_{s}", "", len(pt_bins)-1, array('f', pt_bins)), "pt_Z")
        h.Scale(1./h.Integral())
        hists.append(h)

    tf = ROOT.TFile(f"{hists}z_reweighting.root", "RECREATE")
    for h in hists:
        h.Write()
    tf.Close()

    return 



# function to add z pt weight. TODO: also for GEN or use reco instead?
def weight_zpt(ntuples):
    ROOT.gROOT.ProcessLine('TFile* tf = TFile::Open("z_reweighting.root", "READ");')
    ROOT.gROOT.ProcessLine('TH1D* h_dt = (TH1D*)tf->Get("h_Zboson_pt_DATA");')
    ROOT.gROOT.ProcessLine('TH1D* h_mc = (TH1D*)tf->Get("h_Zboson_pt_MC");')
    ROOT.gROOT.ProcessLine('TH1D* h_ratio = (TH1D*)h_dt->Clone();')
    ROOT.gROOT.ProcessLine('h_ratio->Divide(h_mc)')
   
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s])
        if s== "DATA":
            rdf = rdf.Define("zPtWeight", "1")
        else:
            rdf = rdf.Define("zPtWeight", "h_ratio->GetBinContent(h_ratio->FindBin(pt_Z))")
            
        quants = list(rdf.GetColumnNames())
        rdf.Snapshot("Events", ntuples[s], quants)

    return



def make_ntuples(nanoAODs, ntuples, pt_bins):
    for s in nanoAODs:
        quants = [
            "pt_Z", "mass_Z", "eta_Z", "phi_Z",
            "pt_1", "mass_1", "eta_1", "phi_1", "charge_1",
            "pt_2", "mass_2", "eta_2", "phi_2", "charge_2",
        ]
        # load nanoAOD
        rdf = ROOT.RDataFrame("Events", nanoAODs[s])
        
        # only collect events w/ >1 muon and find muon pair closest to z mass. Muon1 is always charge -1 and muon2 always +1
        rdf = rdf.Filter("Muon_pt.size() > 1")
        rdf = rdf.Define("ind", """ROOT::VecOps::RVec<Int_t> (get_indices(
                         nMuon,
                         &Muon_pt,
                         &Muon_eta,
                         &Muon_phi,
                         &Muon_mass,
                         &Muon_tkRelIso,
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

        # Define quantities of Z boson and collect those events with 50 < m_Z < 130
        rdf = rdf.Define("p4_1", "ROOT::Math::PtEtaPhiMVector(pt_1, eta_1, phi_1, mass_1)")
        rdf = rdf.Define("p4_2", "ROOT::Math::PtEtaPhiMVector(pt_2, eta_2, phi_2, mass_2)")
        rdf = rdf.Define("p4_Z", "p4_1 + p4_2")
        rdf = rdf.Define("pt_Z", "p4_Z.Pt()")
        rdf = rdf.Define("eta_Z", "p4_Z.Eta()")
        rdf = rdf.Define("mass_Z", "p4_Z.M()")
        rdf = rdf.Define("phi_Z", "p4_Z.Phi()")
        
        rdf = rdf.Filter("mass_Z > 50 && mass_Z < 130")
        
        if s=="MC":
            # perform gen delta R matching and collect corresponding events and gen quantities
            rdf = rdf.Define("genind_1", """muon_genmatch(
                                                eta_1,
                                                phi_1,
                                                charge_1,
                                                &GenPart_status,
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
                                                &GenPart_pdgId,
                                                &GenPart_genPartIdxMother,
                                                &GenPart_eta,
                                                &GenPart_phi
                                                )""")
            rdf = rdf.Filter("genind_1 != -99 && genind_2 != -99 && genind_1 != genind_2")
            rdf = rdf.Define("genpt_1", "GenPart_pt[genind_1]")
            rdf = rdf.Define("genpt_2", "GenPart_pt[genind_2]")
            rdf = rdf.Define("geneta_1", "GenPart_eta[genind_1]")
            rdf = rdf.Define("geneta_2", "GenPart_eta[genind_2]")
            rdf = rdf.Define("genphi_1", "GenPart_phi[genind_1]")
            rdf = rdf.Define("genphi_2", "GenPart_phi[genind_2]")
            rdf = rdf.Define("genmass_1", "GenPart_mass[genind_1]")
            rdf = rdf.Define("genmass_2", "GenPart_mass[genind_2]")


            
            quants += [
                "genpt_1", "genmass_1", "geneta_1", "genphi_1",
                "genpt_2", "genmass_2", "geneta_2", "genphi_2"
            ]

        # make output with interesting data
        rdf.Snapshot("Events", ntuples[s], quants)

    return



if __name__=='__main__':
    pt_bins = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 140, 200]
    nanoAODs = {
        'MC': f"{datadir}MC.root",
        'DATA': f"{datadir}DATA.root",
    }
    ntuples = {
        'MC': f"{datadir}MC_ntuples.root",
        'DATA': f"{datadir}DATA_ntuples.root",
    }
    hists = 'hists/'

    ntuple.make_ntuples(nanoAODs, ntuples, pt_bins, hists)
    ntuple.hist_zpt(ntuples, pt_bins)
    ntuple.weight_zpt(ntuples)