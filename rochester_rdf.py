import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array
from argparse import ArgumentParser
import os

# import modules
import python.ntuple as ntuple
import python.scale_corr as corr



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
        action='store_true',
        default=False,
        help="produce ntuples from NanoAOD"
    )
    parser.add_argument(
        '-S',
        '--scale',
        action='store_true',
        default=False,
        help="Make scale correction"
    )
    parser.add_argument(
        '-P',
        '--plot',
        default=False,
        action='store_true',
        help='Make plots'
    )
    args = parser.parse_args()
    return args




def hist_ntuples(ntuples, var, nbins, low, up, hdir):  
    #make histogram of any variable in ntuple (e.g. "mass_Z" )
    hists=[]
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s])
        h_info=(var+"_"+s, var+" "+s, nbins, low, up)
        hists.append(rdf.Histo1D(h_info, var))
    
    tf = ROOT.TFile(f"{hdir}mass_z_hist.root", "RECREATE")
    for h in hists:
        h.Write()
    tf.Close()




if __name__=='__main__':

    pt_bins = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 140, 200]
    oneOverpT_bins = np.linspace(1/200,1/5, 1000000)#[i/2000000. for i in range(200000)]
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
    ntuples_corr = {
        'DATA': f"{datadir}DATA_ntuples_corr.root",
        'MC': f"{datadir}MC_ntuples_corr.root",
    }
    hdir = 'hists/'
    pdir = 'plots/'

    args = parse_args()

    if args.ntuples:
        ntuple.make_ntuples(nanoAODs, ntuples, pt_bins)
        os.makedirs(hdir, exist_ok=True)
        ntuple.hist_zpt(ntuples, pt_bins, hdir)
        ntuple.weight_zpt(ntuples)

    if args.scale:
        os.makedirs(hdir, exist_ok=True)
        corr.hist_oneOverpT(ntuples, oneOverpT_bins, eta_bins, phi_bins, hdir)
        corr.get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir)
        corr.apply_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir)

    if args.plot:
        os.makedirs(pdir, exist_ok=True)
        corr.hist_ntuples(ntuples, "mass_Z", 90, 60, 120, hdir)
