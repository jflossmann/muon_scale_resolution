import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array
import python.ntuple as ntuple
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(
        description="""
            Rochester correction calculator
            *** add more information ***
        """
    )
    parser.add_argument(
        '-N',
        '--ntuples',
        action=store_true,
        default=False,
        description="produce ntuples from NanoAOD"
    )
    parser.add_argument(
        '-S',
        '--scale',
        action=store_true,
        default=False,
        description="Make scale correction"
    )

def hist_oneOverpT(ntuples, oneOverpT_bins, eta_bins, phi_bins):
    # ROOT.gROOT.ProcessLine('tf->Close()')
    hists = []
    ntuples["GEN"] = ntuples["MC"]
    for s in ntuples:
        gen = ""
        if s == "GEN":
            gen == "gen"

        rdf = ROOT.RDataFrame("Events", ntuples[s])
        rdf = rdf.Define("oneOverpT_1", f"1./{gen}pt_1")
        rdf = rdf.Define("oneOverpT_2", f"1./{gen}pt_2")
        h_neg = rdf.Histo3D(
            (
                f"h_oneOverpT_{s}_neg", "", 
                len(eta_bins)-1, array('f', eta_bins), 
                len(phi_bins)-1, array('f', phi_bins), 
                len(oneOverpT_bins)-1, array('f', oneOverpT_bins)
            ),
            f"{gen}eta_1",
            f"{gen}phi_1",
            "oneOverpT_1",
            "z_pT_weight" # TODO: improve method. averaging over bins not precise enough
        )
        hists.append(h_neg)
        h_pos = rdf.Histo3D(
            (
                f"h_oneOverpT_{s}_pos", "",
                len(eta_bins)-1, array('f', eta_bins), 
                len(phi_bins)-1, array('f', phi_bins), 
                len(oneOverpT_bins)-1, array('f', oneOverpT_bins)
            ),
            f"{gen}eta_1",
            f"{gen}phi_1",
            "oneOverpT_2",
            "z_pT_weight"
        )
        hists.append(h_pos)
    tf = ROOT.TFile("oneOverpT.root","RECREATE")
    for h in hists:
        h.Write()
    tf.Close()



def get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins):
    negpos = ["neg", "pos"]

    # get 3D histograms from TFile
    tf = ROOT.TFile("oneOverpT.root", "READ")
    oneOverPt_hists = {}
    for s in list(ntuples.keys())+["GEN"]:
        for np in negpos:
            h_tmp = tf.Get(f"h_oneOverpT_{s}_{np}")
            oneOverPt_hists[f"{s}_mean_{np}"] = h_tmp.Project3DProfile("xy")
            oneOverPt_hists[f"{s}_mean_{np}"].SetDirectory(ROOT.nullptr)
    tf.Close()
    #print(samples)
    # define correction factor C from paper as root 3d histogram
    C = {}
    for s in ntuples:
        print(s)
        C[s] = ROOT.TH3D(
            f"C_{s}", "",
            len(charge_bins)-1, array('f', charge_bins),
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )

        for eta in range(len(eta_bins)-1):
            for phi in range(len(phi_bins)-1):
                for charge in range(len(charge_bins)-1):
                    mean_gen = oneOverPt_hists[f"GEN_mean_{negpos[charge]}"].GetBinContent(eta+1, phi+1)
                    mean = oneOverPt_hists[f"{s}_mean_{negpos[charge]}"].GetBinContent(eta+1, phi+1)
                    C[s].SetBinContent(
                        charge+1, eta+1, phi+1,
                        mean_gen - mean
                    )
                    print(mean_gen, mean)
    tf = ROOT.TFile("C.root", "RECREATE")
    for s in ntuples:
        C[s].Write()
    tf.Close()


def hist_ntuples(ntuples, var, nbins, low, up):  
    #make histogram of any variable in ntuple (e.g. "mass_Z" )
    hists=[]
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s])
        h_info=(var+"_"+s, var+" "+s, nbins, low, up)
        hists.append(rdf.Histo1D(h_info, var))
    
    tf = ROOT.TFile("mass_z_hist.root", "RECREATE")
    for h in hists:
        h.Write()
    tf.Close()



if __name__=='__main__':

    pt_bins = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 140, 200]
    oneOverpT_bins = np.linspace(1/200,1/5, 100000)#[i/2000000. for i in range(200000)]
    eta_bins = [-2.4, -0.8, 0.8, 2.4]
    phi_bins = [-3.2, 0, 3.2]
    charge_bins = [-2,0,2]

    lumi = 31906.78
    xsec = 5600
    
    datadir = "/ceph/jdriesch/rochester/"
    nanoAODs = {
        'DATA': f"{datadir}DATA.root",
        'MC': f"{datadir}MC.root",
    }
    ntuples = {
        'DATA': f"{datadir}DATA_ntuples.root",
        'MC': f"{datadir}MC_ntuples.root",
    }

    args = parse_args()

    if args.ntuples:
        ntuple.make_ntuples(nanoAODs, ntuples, pt_bins)
        ntuple.hist_zpt(ntuples, pt_bins)
        ntuple.weight_zpt(ntuples)

    if args.scale:
        hist_oneOverpT(ntuples, oneOverpT_bins, eta_bins, phi_bins)
        get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins)
        hist_ntuples(ntuples, "mass_Z", 90, 60, 120)
