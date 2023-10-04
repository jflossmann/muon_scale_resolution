import ROOT
from array import array

def hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir):
    # ROOT.gROOT.ProcessLine('tf->Close()')
    hists = []
    ntuples_tmp = dict(ntuples)
    ntuples_tmp["GEN"] = ntuples_tmp["MC"]
    for s in ntuples_tmp:
        gen = ""
        if s == "GEN":
            gen == "gen"

        rdf = ROOT.RDataFrame("Events", ntuples_tmp[s])
        rdf = rdf.Define("oneOverPt_1", f"1./{gen}pt_1")
        rdf = rdf.Define("oneOverPt_2", f"1./{gen}pt_2")
        h_neg = rdf.Histo3D(
            (
                f"h_oneOverPt_{s}_neg", "", 
                len(eta_bins)-1, array('f', eta_bins), 
                len(phi_bins)-1, array('f', phi_bins), 
                len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
            ),
            f"{gen}eta_1",
            f"{gen}phi_1",
            "oneOverPt_1",
            "zPtWeight" # TODO: improve method. averaging over bins not precise enough
        )
        hists += [h_neg, h_neg.Project3DProfile("xy")]
        h_pos = rdf.Histo3D(
            (
                f"h_oneOverPt_{s}_pos", "",
                len(eta_bins)-1, array('f', eta_bins), 
                len(phi_bins)-1, array('f', phi_bins), 
                len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
            ),
            f"{gen}eta_2",
            f"{gen}phi_2",
            "oneOverPt_2",
            "zPtWeight"
        )
        hists += [h_pos, h_pos.Project3DProfile("xy")]
    tf = ROOT.TFile(f"{hdir}oneOverPt.root","RECREATE")
    for h in hists:
        h.Write()
    tf.Close()



def get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir):
    negpos = ["neg", "pos"]

    # get 3D histograms from TFile
    tf = ROOT.TFile(f"{hdir}oneOverPt.root", "READ")
    oneOverPt_hists = {}
    for s in list(ntuples.keys())+["GEN"]:
        for np in negpos:
            oneOverPt_hists[f"{s}_mean_{np}"] = tf.Get(f'h_oneOverPt_{s}_{np}_pxy')
            oneOverPt_hists[f"{s}_mean_{np}"].SetDirectory(ROOT.nullptr)
    tf.Close()

    # define correction factor C from paper as root 3d histogram
    C, Dm, Da, M, A = {}, {}, {}, {}, {}
    for s in ntuples:
        #print(s)
        C[s] = ROOT.TH3D(
            f"C_{s}", "",
            len(charge_bins)-1, array('f', charge_bins),
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        Dm[s] = ROOT.TH2D(
            f"Dm_{s}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        Da[s] = ROOT.TH2D(
            f"Da_{s}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        M[s] = ROOT.TH2D(
            f"M_{s}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        A[s] = ROOT.TH2D(
            f"A_{s}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )

        for eta in range(len(eta_bins)-1):
            for phi in range(len(phi_bins)-1):
                for charge in range(len(charge_bins)-1):
                    mean_gen = oneOverPt_hists[f"GEN_mean_{negpos[charge]}"].GetBinContent(phi+1, eta+1)
                    mean = oneOverPt_hists[f"{s}_mean_{negpos[charge]}"].GetBinContent(phi+1, eta+1)
                    C[s].SetBinContent(
                        charge+1,
                        eta+1,
                        phi+1,
                        mean_gen - mean
                    )
                    #print(mean_gen, mean)
                Dm[s].SetBinContent(
                    eta+1,
                    phi+1,
                    (C[s].GetBinContent(2, eta+1, phi+1) + C[s].GetBinContent(1, eta+1, phi+1)) / 2.
                )
                Da[s].SetBinContent(
                    eta+1,
                    phi+1,
                    (C[s].GetBinContent(2, eta+1, phi+1) - C[s].GetBinContent(1, eta+1, phi+1)) / 2.
                )
                M[s].SetBinContent(
                    eta+1,
                    phi+1,
                    1 + 2*Dm[s].GetBinContent(eta+1, phi+1) / (
                        oneOverPt_hists[f"{s}_mean_{negpos[0]}"].GetBinContent(phi+1, eta+1) + 
                        oneOverPt_hists[f"{s}_mean_{negpos[1]}"].GetBinContent(phi+1, eta+1)
                    )
                )

                A[s].SetBinContent(
                    eta+1,
                    phi+1,
                    Da[s].GetBinContent(eta+1, phi+1) - (
                        Dm[s].GetBinContent(eta+1, phi+1)*(
                        oneOverPt_hists[f"{s}_mean_{negpos[0]}"].GetBinContent(phi+1, eta+1) -
                        oneOverPt_hists[f"{s}_mean_{negpos[1]}"].GetBinContent(phi+1, eta+1)
                    ) / (
                        oneOverPt_hists[f"{s}_mean_{negpos[0]}"].GetBinContent(phi+1, eta+1) +
                        oneOverPt_hists[f"{s}_mean_{negpos[1]}"].GetBinContent(phi+1, eta+1)
                    )
                    )
                )

    tf = ROOT.TFile(f"{hdir}C.root", "RECREATE")
    for s in ntuples:
        C[s].Write()
        Dm[s].Write()
        Da[s].Write()
        M[s].Write()
        A[s].Write()
    tf.Close()



def apply_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir):
    
    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}C.root", "READ");')
    ROOT.gROOT.ProcessLine('TH2D* M_DATA = (TH2D*)tf->Get("M_DATA");')
    ROOT.gROOT.ProcessLine('TH2D* M_MC = (TH2D*)tf->Get("M_MC");')
    ROOT.gROOT.ProcessLine('TH2D* A_DATA = (TH2D*)tf->Get("A_DATA");')
    ROOT.gROOT.ProcessLine('TH2D* A_MC = (TH2D*)tf->Get("A_MC");')   
    
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s])
        rdf = rdf.Define("oneOverPt_1", "1./pt_1")
        rdf = rdf.Define(
            "oneOverPt_1_roccor",
            f"""oneOverPt_1 * M_{s}->GetBinContent(M_{s}->FindBin(eta_1, phi_1)) + 
            charge_1 * A_{s}->GetBinContent(A_{s}->FindBin(eta_1, phi_1))
            """
        )
        rdf = rdf.Define("pt_1_roccor", "1./oneOverPt_1_roccor")
        rdf = rdf.Define("oneOverPt_2", "1./pt_2")
        rdf = rdf.Define(
            "oneOverPt_2_roccor",
            f"""oneOverPt_2 * M_{s}->GetBinContent(M_{s}->FindBin(eta_2, phi_2)) + 
            charge_2 * A_{s}->GetBinContent(A_{s}->FindBin(eta_2, phi_2))
            """
        )
        rdf = rdf.Define("pt_2_roccor", "1./oneOverPt_2_roccor")

        quants = list(rdf.GetColumnNames())
        
        rdf = rdf.Define("p4_1", "ROOT::Math::PtEtaPhiMVector(pt_1_roccor, eta_1, phi_1, mass_1)")
        rdf = rdf.Define("p4_2", "ROOT::Math::PtEtaPhiMVector(pt_2_roccor, eta_2, phi_2, mass_2)")
        rdf = rdf.Define("p4_Z", "p4_1 + p4_2")
        rdf = rdf.Define("mass_Z_roccor", "p4_Z.M()")
        
        rdf.Snapshot("Events", ntuples[s].replace('.root', '_corr.root'), quants + ["mass_Z_roccor"])




def scale_closure(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir):
    hists = []
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s].replace('.root', '_corr.root'))
        h_neg = rdf.Histo3D(
            (
                f"h_oneOverPt_{s}_neg_roccor", "", 
                len(eta_bins)-1, array('f', eta_bins), 
                len(phi_bins)-1, array('f', phi_bins), 
                len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
            ),
            f"eta_1",
            f"phi_1",
            "oneOverPt_1_roccor",
            "zPtWeight" # TODO: improve method. averaging over bins not precise enough
        )

        hists += [h_neg, h_neg.Project3DProfile("xy")]

        h_pos = rdf.Histo3D(
            (
                f"h_oneOverPt_{s}_pos_roccor", "",
                len(eta_bins)-1, array('f', eta_bins), 
                len(phi_bins)-1, array('f', phi_bins), 
                len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
            ),
            f"eta_2",
            f"phi_2",
            "oneOverPt_2_roccor",
            "zPtWeight"
        )

        hists += [h_pos, h_pos.Project3DProfile("xy")]

    tf = ROOT.TFile(f"{hdir}oneOverPt_roccor.root","RECREATE")
    for h in hists:
        h.Write()
    tf.Close()