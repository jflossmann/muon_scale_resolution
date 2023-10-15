import ROOT
from array import array
from tqdm import tqdm
from multiprocessing import Pool, RLock
import os

#define different statistics to base corrections on
modes = ['mean', 'median']    

def hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='')->None:
    """create histograms for inverse transversal momentum from ntuple data"""
    
    hists = {}
    hists["tosave"] = []
    negpos = ["neg", "pos"]

    #iterate over data sources
    for typ in ntuples:
        hists[typ+'0'] = []
        hists[typ+'1'] = []

        for sample in ntuples[typ]:
            gen = ""
            roccor = corr
            if typ == "GEN":
                gen = "gen"
                roccor = ''
            #create RDataframe to acess data
            rdf = ROOT.RDataFrame("Events", ntuples[typ][sample])
            print(ntuples[typ][sample])
            rdf = rdf.Define("weight", "zPtWeight*genWeight*sumwWeight*xsec")

            for np in range(2):
                #define new column for 1/pt
                rdf = rdf.Define(f"oneOverPt_{np+1}", f"1./{gen}pt_{np+1}{roccor}")
                #create 3D histogram
                h_3d = rdf.Histo3D(
                    (
                        f"h_oneOverPt_{sample}_{negpos[np]}", "", 
                        len(eta_bins)-1, array('f', eta_bins), 
                        len(phi_bins)-1, array('f', phi_bins), 
                        len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
                    ),
                    f"{gen}eta_{np+1}",
                    f"{gen}phi_{np+1}",
                    f"oneOverPt_{np+1}",
                    "weight"
                )
                hists[typ+str(np)] += [h_3d]

        for np in range(2):
            h_sum = hists[typ+str(np)][0].Clone(f"h_oneOverPt_{typ}_{negpos[np]}")
            for i in range(len(hists[typ+str(np)])-1):
                h_tmp = hists[typ+str(np)][i+1].Clone("h_tmp")
                h_sum.Add(h_tmp)

            #create histogram for median
            h_median = ROOT.TH2D(
                f"h_oneOverPt_{typ}_median_{negpos[np]}", "",
                len(eta_bins)-1, array('d', eta_bins), 
                len(phi_bins)-1, array('d', phi_bins)
            )

            quantile = array('d', [0.5])
            median = array('d', [0])

            for e in range(len(eta_bins)-1):
                for p in range(len(phi_bins)-1):
                    h_ep = h_sum.ProjectionZ(f"h_tmp_{e}_{p}", e+1, e+1, p+1, p+1)

                    h_ep.GetQuantiles(1, median, quantile)
                    h_median.SetBinContent(e+1, p+1, median[0])

            hists["tosave"] += [h_sum, h_sum.Project3DProfile("yx"), h_median]            
        
    #save
    tf = ROOT.TFile(f"{hdir}oneOverPt{corr}.root","RECREATE")
    for h in hists["tosave"]:
        h.Write()
    tf.Close()



def get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir)->None: 
    """extract scale corrections from ntuple data"""

    #charges
    negpos = ["neg", "pos"]
    
    # get 3D histograms from TFile
    tf = ROOT.TFile(f"{hdir}oneOverPt.root", "READ")
    oneOverPt_hists = {}

    #iterate over data, mc and gen
    for typ in ntuples:
        for np in negpos:
            #read mean-histogram
            oneOverPt_hists[f"{typ}_mean_{np}"] = tf.Get(f'h_oneOverPt_{typ}_{np}_pyx')
            oneOverPt_hists[f"{typ}_mean_{np}"].SetDirectory(ROOT.nullptr)
            #read median histogram
            oneOverPt_hists[f"{typ}_median_{np}"] = tf.Get(f'h_oneOverPt_{typ}_median_{np}')
            oneOverPt_hists[f"{typ}_median_{np}"].SetDirectory(ROOT.nullptr)
    tf.Close()

    # define correction factor C from paper as root 3d histogram
    C, Dm, Da, M, A = {}, {}, {}, {}, {}
    for mode in modes:
        for typ in list(ntuples.keys())[:-1]:
            #print(s)
            #make 3D-histograms for C, Dm, Da, M and A
            C[mode+typ] = ROOT.TH3D(
                f"C_{mode}_{typ}", "",
                len(charge_bins)-1, array('f', charge_bins),
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            Dm[mode+typ] = ROOT.TH2D(
                f"Dm_{mode}_{typ}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            Da[mode+typ] = ROOT.TH2D(
                f"Da_{mode}_{typ}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            M[mode+typ] = ROOT.TH2D(
                f"M_{mode}_{typ}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            A[mode+typ] = ROOT.TH2D(
                f"A_{mode}_{typ}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            #iterate over eta, phi and charge
            for eta in range(len(eta_bins)-1):
                for phi in range(len(phi_bins)-1):
                    for charge in range(len(charge_bins)-1):
                        mean_gen = oneOverPt_hists[f"GEN_{mode}_{negpos[charge]}"].GetBinContent(eta+1, phi+1)
                        mean = oneOverPt_hists[f"{typ}_{mode}_{negpos[charge]}"].GetBinContent(eta+1, phi+1)

                        #set value of C[mode+s] for current bin
                        C[mode+typ].SetBinContent(
                            charge+1,
                            eta+1,
                            phi+1,
                            mean_gen - mean
                        )
                        #print(mean_gen, mean)
                    #calculate Dm, Da, M, A for current bin
                    Dm[mode+typ].SetBinContent(
                        eta+1,
                        phi+1,
                        (C[mode+typ].GetBinContent(2, eta+1, phi+1) + C[mode+typ].GetBinContent(1, eta+1, phi+1)) / 2.
                    )
                    Da[mode+typ].SetBinContent(
                        eta+1,
                        phi+1,
                        (C[mode+typ].GetBinContent(2, eta+1, phi+1) - C[mode+typ].GetBinContent(1, eta+1, phi+1)) / 2.
                    )
                    M[mode+typ].SetBinContent(
                        eta+1,
                        phi+1,
                        1 + 2*Dm[mode+typ].GetBinContent(eta+1, phi+1) / (
                            oneOverPt_hists[f"{typ}_{mode}_{negpos[0]}"].GetBinContent(eta+1, phi+1) + 
                            oneOverPt_hists[f"{typ}_{mode}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                        )
                    )

                    A[mode+typ].SetBinContent(
                        eta+1,
                        phi+1,
                        Da[mode+typ].GetBinContent(eta+1, phi+1) - (
                            Dm[mode+typ].GetBinContent(eta+1, phi+1)*(
                            oneOverPt_hists[f"{typ}_{mode}_{negpos[0]}"].GetBinContent(eta+1, phi+1) -
                            oneOverPt_hists[f"{typ}_{mode}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                        ) / (
                            oneOverPt_hists[f"{typ}_{mode}_{negpos[0]}"].GetBinContent(eta+1, phi+1) +
                            oneOverPt_hists[f"{typ}_{mode}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                        )
                        )
                    )
    #safe histograms as Tfile
    tf = ROOT.TFile(f"{hdir}C.root", "RECREATE")
    for mode in modes:
        for typ in list(ntuples.keys())[:-1]:
            C[mode+typ].Write()
            Dm[mode+typ].Write()
            Da[mode+typ].Write()
            M[mode+typ].Write()
            A[mode+typ].Write()
    tf.Close()



def apply_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir):
    #open Tfile with scale corrections
    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}C.root", "READ");')
    for mode in modes:
        #read histograms for curent mode
        ROOT.gROOT.ProcessLine(f'TH2D* M_{mode}_DATA = (TH2D*)tf->Get("M_{mode}_DATA");')
        ROOT.gROOT.ProcessLine(f'TH2D* M_{mode}_MC = (TH2D*)tf->Get("M_{mode}_MC");')
        ROOT.gROOT.ProcessLine(f'TH2D* A_{mode}_DATA = (TH2D*)tf->Get("A_{mode}_DATA");')
        ROOT.gROOT.ProcessLine(f'TH2D* A_{mode}_MC = (TH2D*)tf->Get("A_{mode}_MC");')   
    
    for typ in list(ntuples.keys())[:-1]:
        for sample in ntuples[typ]:
            #make RDataFrame object for current dataset
            rdf = ROOT.RDataFrame("Events", ntuples[typ][sample])
            quants = list(rdf.GetColumnNames())
            #calculate reciprocal of momentum for both muons
            rdf = rdf.Define("oneOverPt_1", "1./pt_1")
            rdf = rdf.Define("oneOverPt_2", "1./pt_2")

            
            for mode in modes:
                #calculate the corrected reciprocal pt
                rdf = rdf.Define(
                    f"oneOverPt_1_{mode}_roccor",
                    f"""oneOverPt_1 * M_{mode}_{typ}->GetBinContent(M_{mode}_{typ}->FindBin(eta_1, phi_1)) + 
                    charge_1 * A_{mode}_{typ}->GetBinContent(A_{mode}_{typ}->FindBin(eta_1, phi_1))
                    """
                )
                rdf = rdf.Define(f"pt_1_{mode}_roccor", f"1./oneOverPt_1_{mode}_roccor")
                rdf = rdf.Define(
                    f"oneOverPt_2_{mode}_roccor",
                    f"""oneOverPt_2 * M_{mode}_{typ}->GetBinContent(M_{mode}_{typ}->FindBin(eta_2, phi_2)) + 
                    charge_2 * A_{mode}_{typ}->GetBinContent(A_{mode}_{typ}->FindBin(eta_2, phi_2))
                    """
                )
                rdf = rdf.Define(f"pt_2_{mode}_roccor", f"1./oneOverPt_2_{mode}_roccor")

                #calculate corrected 4-momenta and corrected Z-mass
                rdf = rdf.Define(f"p4_1_{mode}", f"ROOT::Math::PtEtaPhiMVector(pt_1_{mode}_roccor, eta_1, phi_1, mass_1)")
                rdf = rdf.Define(f"p4_2_{mode}", f"ROOT::Math::PtEtaPhiMVector(pt_2_{mode}_roccor, eta_2, phi_2, mass_2)")
                rdf = rdf.Define(f"p4_Z_{mode}", f"p4_1_{mode} + p4_2_{mode}")
                rdf = rdf.Define(f"mass_Z_{mode}_roccor", f"p4_Z_{mode}.M()")

                #add calculated quantities to quantities
                quants += [f"pt_1_{mode}_roccor", f"pt_2_{mode}_roccor", f"mass_Z_{mode}_roccor"]
                
            #save everything in new ntuple
            rdf.Snapshot("Events", ntuples[typ][sample].replace('.root', '_corr.root'), quants)
