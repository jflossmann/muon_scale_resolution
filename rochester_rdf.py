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
import python.res_corr as rc
import python.iterative as it

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
        '-G',
        '--gen_cor',
        default=False,
        action='store_true',
        help='Correct GEN-ntuple-reco values with MC-SIG kappa and lambda from 1/pt correction'
    parser.add_argument(
        '-R',
        '--res',
        default=False,
        action='store_true',
        help='Make Resolution calculation'
    )
    parser.add_argument(
    '-I',
    '--iterative',
    default=False,
    action='store_true',
    help='Make iterative correction'
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
    diffpt_bins = np.linspace(-0.01, 0.01, 200)
    pt_bins = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 200]
    oneOverPt_bins = np.linspace(0.019,0.027, 50)#[i/2000000. for i in range(200000)]
    eta_bins = [-2.4, -2.1, -1.85, -0.4, 0, 0.4, 1.85, 2.1, 2.4]
    phi_bins = np.linspace(-3.2, 3.2, 17) #[-3.2, -2.4, -1.6, -.8, 0, .8, 1.6, 2.4, 3.2]
    charge_bins = [-2,0,2]
    mass_bins = np.linspace(75, 105, 61)
    
    datadir = "/ceph/jdriesch/rochester/"
    nanoAODs = ntuple.yaml_loader('data/nanoAODs.yaml')
    datasets = ntuple.yaml_loader('data/datasets.yaml')
    samples_to_plot = ["DATA", "SIG", "GEN"]
    ntuples = {
        'DATA': {'DATA': f"{datadir}DATA_ntuples.root"},
        'SIG': {'SIG': f"{datadir}DY_ntuples.root"},
        'BKG': {            
            'WW': f"{datadir}WW_ntuples.root",
            'WZ': f"{datadir}WZ_ntuples.root",
            'ZZ': f"{datadir}ZZ_ntuples.root",
            'TT': f"{datadir}TT_ntuples.root",
        },
        'GEN': {'GEN': f"{datadir}GEN_ntuples.root"}
    }
    ntuples_zPt = {
        'DATA': {'DATA': f"{datadir}DATA_ntuples_zPt.root"},
        'SIG': {'SIG': f"{datadir}DY_ntuples_zPt.root"},
        'BKG': {
            'WW': f"{datadir}WW_ntuples_zPt.root",
            'WZ': f"{datadir}WZ_ntuples_zPt.root",
            'ZZ': f"{datadir}ZZ_ntuples_zPt.root",
            'TT': f"{datadir}TT_ntuples_zPt.root",
        },
        'GEN': {'GEN': f"{datadir}GEN_ntuples_zPt.root"}
    }
    ntuples_corr = {
        'DATA': {'DATA': f"{datadir}DATA_ntuples_zPt_corr.root"},
        'SIG': {'SIG': f"{datadir}DY_ntuples_zPt_corr.root"},
        'BKG': {
            'WW': f"{datadir}WW_ntuples_zPt_corr.root",
            'WZ': f"{datadir}WZ_ntuples_zPt_corr.root",
            'ZZ': f"{datadir}ZZ_ntuples_zPt_corr.root",
            'TT': f"{datadir}TT_ntuples_zPt_corr.root",
        },
        'GEN': {'GEN': f"{datadir}GEN_ntuples_zPt.root"}
    }

    sf_path = 'data/scaleFactors/Run2/UL/2018/2018_Z/Efficiencies_muon_generalTracks_Z_Run2018_UL_'
    SFs = {
        'ID': ntuple.load_hist(sf_path+'ID.root', 'NUM_MediumID_DEN_TrackerMuons_abseta_pt'),
        'ISO': ntuple.load_hist(sf_path+'ISO.root', 'NUM_TightRelIso_DEN_MediumID_abseta_pt')
    }
    hdir = 'hists/'
    pdir = 'plots/'

    args = parse_args()

    if args.ntuples:
        # ntuple.make_ntuples(nanoAODs, datasets, datadir)
        os.makedirs(hdir, exist_ok=True)
        ntuple.hist_zpt(ntuples, pt_bins, hdir)
        ntuple.weight_zpt(ntuples, hdir, diffpt_bins, eta_bins, phi_bins) # weight also backgrounds and GEN?

    if args.scale:
        os.makedirs(hdir, exist_ok=True)
        corr.hist_oneOverpT(ntuples_zPt, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir)
        corr.get_scale_corrections(samples_to_plot, eta_bins, phi_bins, charge_bins, hdir)
        corr.apply_scale_corrections(ntuples_zPt, hdir)
        corr.hist_oneOverpT(ntuples_corr, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='_mean_roccor')
        corr.hist_oneOverpT(ntuples_corr, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='_median_roccor')

    if args.plot:
        os.makedirs(pdir, exist_ok=True)
        plot.hist_ntuples(ntuples_corr, "mass_Z", len(mass_bins)-1, mass_bins[0], mass_bins[-1], hdir, "mass_z", "RECREATE")
        for corr in ['_mean_roccor', '_median_roccor']:
            plot.hist_ntuples(ntuples_corr, "mass_Z", len(mass_bins)-1, mass_bins[0], mass_bins[-1], hdir, "mass_z", "update", corr)
        plot.plot_stuff(pdir, eta_bins, phi_bins)
        
    if args.zmass:
        zmass.hist_zmass(ntuples_corr, eta_bins, phi_bins, mass_bins, hdir)
        zmass.fit_zmass(eta_bins, phi_bins, hdir)
        zmass.plot_zmass(eta_bins, phi_bins, hdir)

    if args.gen_cor:
        rc.cor_gen(ntuples_corr["GEN"]["GEN"],ntuples_corr['SIG']['SIG'], eta_bins=eta_bins, phi_bins=phi_bins)
        ntuples_corr["GEN"]["GEN"]=f"{datadir}GEN_zPt_corr.root"
        
    if args.res:
        pull_bins=np.linspace(-5,5,100)
        abseta_bins=np.linspace(0, 2.4, 13)
        nl_bins=[6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5]
        pt_bins=[25,30,35,40,50,60,80,110,150,200]
        rc.get_res_correction(ntuples_corr["GEN"]["GEN"], pull_bins, abseta_bins, nl_bins, pt_bins, pdir, hdir, do_plot=True)
        rc.apply_res_corr(ntuples_corr["GEN"]["GEN"], hdir, pdir, do_plot=True)

    if args.iterative:
        ntuples_corr["GEN"]["GEN"]=f"{datadir}GEN_ntuples_zPt_corr.root"
        it.iterative_correction(samples=ntuples_corr, eta_bins=eta_bins, phi_bins=phi_bins, hdir=hdir, pdir=pdir)
