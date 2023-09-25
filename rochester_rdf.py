import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array

# function in c++ code. Finds indices of muons closest to z boson mass
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<Int_t> get_indices(UInt_t nMuon, ROOT::VecOps::RVec<Float_t> *Muon_pt, ROOT::VecOps::RVec<Float_t> *Muon_eta, ROOT::VecOps::RVec<Float_t> *Muon_phi, ROOT::VecOps::RVec<Float_t> *Muon_mass){
        Float_t deltaM = 1000;
        Int_t ind1, ind2;
        ind1 = -99;
        ind2 = -99;
        for(int i=0; i<nMuon; i++){
            for(int j=i; j<nMuon; j++){
                TLorentzVector mui, muj, Z;
                mui.SetPtEtaPhiM(Muon_pt->at(i), Muon_eta->at(i), Muon_phi->at(i), Muon_mass->at(i));
                muj.SetPtEtaPhiM(Muon_pt->at(j), Muon_eta->at(j), Muon_phi->at(j), Muon_mass->at(j));
                Z = mui + muj;
                if (fabs(Z.M() - 91.1876) < deltaM){
                    deltaM = fabs(Z.M() - 91.1876);
                    ind1 = i;
                    ind2 = j;
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
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<Int_t> get_hlt_matches(UInt_t nMuon, ROOT::VecOps::RVec<Float_t> *Muon_eta, ROOT::VecOps::RVec<Float_t> *Muon_phi, ROOT::VecOps::RVec<Float_t> *hlt_eta, ROOT::VecOps::RVec<Float_t> *hlt_phi, ROOT::VecOps::RVec<Float_t> *Muon_tkRelIso, ROOT::VecOps::RVec<Float_t> *Muon_mediumPromptId){
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
""")


# function in c++ code. Finds indices of muons matched to GEN
ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<Int_t> get_gen_matches(UInt_t nMuon, ROOT::VecOps::RVec<Float_t> *Muon_eta, ROOT::VecOps::RVec<Float_t> *Muon_phi, ROOT::VecOps::RVec<Float_t> *GenPart_statusFlags, ROOT::VecOps::RVec<Float_t> *GenPart_pdgId, ROOT::VecOps::RVec<Float_t> *GenPart_genPartIdxMother, ROOT::VecOps::RVec<Float_t> *GenPart_eta, ROOT::VecOps::RVec<Float_t> *GenPart_phi, , ROOT::VecOps::RVec<Float_t> *Muon_tkRelIso, ROOT::VecOps::RVec<Float_t> *Muon_mediumPromptId){
        ROOT::VecOps::RVec<Int_t> muon_gen_match;
        for(int i=0; i<nMuon; i++){
            Float_t deltaR=1000;
            if (Muon_pt->at(i) < 10) continue;
            if(fabs(Muon_eta->at(i) > 2.4)) continue;
            if(Muon_tkRelIso->at(i) > 0.1) continue;
            if(Muon_mediumPromptId->at(i)) continue;
            for(int j=0; j<GenPart_eta.size(); j++){
                if (GenPart_statusFlags->at(j) != 1) continue;
                if (fabs(GenPart_pdgId->at(j)) != 13) continue;
                if (GenPart_pdgId->at(GenPart_genPartIdxMother->at(j)) != 23) continue;

                dEta = Muon_eta->at(i) - GenPart_eta->at(j);
                dPhi = Muon_phi->at(i) - GenPart_phi->at(j);
                deltaR = sqrt(dEta*dEta + dPhi*dPhi);
                if(deltaR < 0.1) muon_gen_match.push_back(i);
            }
        }
    return muon_gen_match;
    }
""")



# start of function to extract z pt distribution
def hist_zpt(samples, pt_bins):
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)
    hists = []

    for s in samples[:2]:
        rdf = ROOT.RDataFrame("Events", f"/ceph/jdriesch/rochester/{s}.root")

        rdf = rdf.Filter("Muon_pt.size() > 1")
        rdf = rdf.Define("ind", "ROOT::VecOps::RVec<Int_t> (get_indices(nMuon, &Muon_pt, &Muon_eta, &Muon_phi, &Muon_mass))")
        rdf = rdf.Define("ind0", "ind[0]").Define("ind1", "ind[1]")
        rdf = rdf.Filter("ind0 + ind1 > 0")
        rdf = rdf.Define("pt_1", "Muon_pt[ind[0]]").Define("pt_2", "Muon_pt[ind[1]]")
        rdf = rdf.Define("mass_1", "Muon_mass[ind[0]]").Define("mass_2", "Muon_mass[ind[1]]")
        rdf = rdf.Define("eta_1", "Muon_eta[ind[0]]").Define("eta_2", "Muon_eta[ind[1]]")
        rdf = rdf.Define("phi_1", "Muon_phi[ind[0]]").Define("phi_2", "Muon_phi[ind[1]]")

        rdf = rdf.Define("p4_1", "ROOT::Math::PtEtaPhiMVector(pt_1, eta_1, phi_1, mass_1)")
        rdf = rdf.Define("p4_2", "ROOT::Math::PtEtaPhiMVector(pt_2, eta_2, phi_2, mass_2)")
        rdf = rdf.Define("p4_Z", "p4_1 + p4_2")
        rdf = rdf.Define("pt_Z", "p4_Z.Pt()")

        h = rdf.Histo1D((f"h_Zboson_pt_{s}", "", len(pt_bins)-1, array('f', pt_bins)), "pt_Z")
        hists.append(h)

    tf = ROOT.TFile("z_reweighting.root", "RECREATE")
    for h in hists:
        h.Write()
    tf.Close()

if __name__=='__main__':

    pt_bins = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 140, 200]
    eta_bins = [-2.4, -0.8, 0.8, 2.4]
    phi_bins = [-3.2, 0, 3.2]

    lumi = 31906.78
    xsec = 5600

    samples = ['MC', 'DATA', 'GEN']

    hist_zpt(samples, pt_bins)
