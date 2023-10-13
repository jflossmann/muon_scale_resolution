import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array
from argparse import ArgumentParser
import os

# import modules
import python.ntuple as ntuple
import python.scale_corr as corr
import python.plot as plot
import python.zmass as zmass

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
    parser.add_argument(
        '-Z',
        '--zmass',
        default=False,
        action='store_true',
        help='Check Z mass'
    )
    args = parser.parse_args()
    return args



if __name__=='__main__':

    pt_bins = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 200]
    oneOverPt_bins = np.linspace(1/200,1/20, 200)#[i/2000000. for i in range(200000)]
    eta_bins = [-2.4, -2.1, -1.85, -0.4, 0, 0.4, 1.85, 2.1, 2.4]
    phi_bins = np.linspace(-3.2, 3.2, 17) #[-3.2, -2.4, -1.6, -.8, 0, .8, 1.6, 2.4, 3.2]
    charge_bins = [-2,0,2]
    mass_bins = np.linspace(75, 105, 61)

    lumi = 31906.78
    xsec = 5600
    
    datadir = "/ceph/jdriesch/rochester/"
    nanoAODs = ntuple.yaml_loader('configs/nanoAODs.yaml')
    datasets = ntuple.yaml_loader('configs/datasets.yaml')
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
        ntuple.make_ntuples(nanoAODs, datasets, ntuples, pt_bins)
        os.makedirs(hdir, exist_ok=True)
        ntuple.hist_zpt(ntuples, pt_bins, hdir)
        ntuple.weight_zpt(ntuples, hdir)

    if args.scale:
        os.makedirs(hdir, exist_ok=True)
        corr.hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir)
        corr.get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir)
        corr.apply_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir)
        corr.hist_oneOverpT(ntuples_corr, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='_peak_roccor')

    if args.plot:
        os.makedirs(pdir, exist_ok=True)
        plot.hist_ntuples(ntuples_corr, "mass_Z", len(mass_bins)-1, mass_bins[0], mass_bins[-1], hdir, "mass_z", "RECREATE")
        for corr in ['_mean_roccor', '_median_roccor', '_peak_roccor']:
            plot.hist_ntuples(ntuples_corr, "mass_Z"+corr, len(mass_bins)-1, mass_bins[0], mass_bins[-1], hdir, "mass_z")
        plot.plot_stuff(pdir, eta_bins, phi_bins)
        
    if args.zmass:
        zmass.hist_zmass(ntuples_corr, eta_bins, phi_bins, mass_bins, hdir)
        zmass.fit_zmass(eta_bins, phi_bins, hdir)
        zmass.plot_zmass(eta_bins, phi_bins, hdir)