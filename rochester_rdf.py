import ROOT
import numpy as np
import matplotlib.pyplot as plt


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

# start of function to extract z pt distribution
def hist_zpt(samples):
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)

    rdf = ROOT.RDataFrame("Events", "/ceph/jdriesch/rochester/MC.root")

    rdf = rdf.Filter("Muon_pt.size() > 1")
    rdf = rdf.Define("ind", "ROOT::VecOps::RVec<Int_t> (get_indices(nMuon, &Muon_pt, &Muon_eta, &Muon_phi, &Muon_mass))")
    rdf = rdf.Define("ind0", "ind[0]").Define("ind1", "ind[1]")
    rdf = rdf.Filter("ind0 + ind1 > 0")
    rdf = rdf.Define("pt_1", "Muon_pt[ind[0]]").Define("pt_2", "Muon_pt[ind[1]]")




if __name__=='__main__':

    pt_bins = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 140, 200]
    eta_bins = [-2.4, -0.8, 0.8, 2.4]
    phi_bins = [-3.2, 0, 3.2]

    lumi = 31906.78
    xsec = 5600

    samples = ['MC', 'DATA', 'GEN']

    hist_zpt(samples)
