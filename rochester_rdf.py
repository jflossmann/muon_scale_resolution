import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array
from argparse import ArgumentParser
import os

# import modules
import python.ntuple as ntuple
import python.step1 as step1
import python.plot as plot
import python.zmass as zmass
import python.step2 as step2
import python.step3 as step3
import python.step4 as step4

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
        '-1',
        '--scale',
        action='store_true',
        default=False,
        help="Make scale correction"
    )
    parser.add_argument(
        '-2',
        '--resolution',
        default=False,
        action='store_true',
        help='Make Resolution calculation'
    )
    parser.add_argument(
        '-3',
        '--iterative',
        default=False,
        action='store_true',
        help='Make iterative correction'
    )
    parser.add_argument(
        '-4',
        '--residual',
        default=False,
        action='store_true',
        help='Residual_correction'
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
    
    datadir = "/ceph/jdriesch/rochester/ntuples/"
    nanoAODs = ntuple.yaml_loader('data/nanoAODs.yaml')
    datasets = ntuple.yaml_loader('data/datasets.yaml')
    samples_to_plot = ["DATA", "SIG", "GEN"]
    ntuples = {
        'DATA': {'DATA': f"{datadir}DATA_*.root"},
        'SIG': {'SIG': f"{datadir}DY_*.root"},
        'BKG': {            
            'WW': f"{datadir}WW_*.root",
            'WZ': f"{datadir}WZ_*.root",
            'ZZ': f"{datadir}ZZ_*.root",
            'TT': f"{datadir}TT_*.root",
        },
        'GEN': {'GEN': f"{datadir}GEN_*.root"}
    }
    ntuples_zPt = {
        'DATA': {'DATA': f"{datadir}DATA_zPt.root"},
        'SIG': {'SIG': f"{datadir}DY_zPt.root"},
        'BKG': {
            'WW': f"{datadir}WW_zPt.root",
            'WZ': f"{datadir}WZ_zPt.root",
            'ZZ': f"{datadir}ZZ_zPt.root",
            'TT': f"{datadir}TT_zPt.root",
        },
        'GEN': {'GEN': f"{datadir}GEN_zPt.root"}
    }
    ntuples_corr = {
        'DATA': {'DATA': f"{datadir}DATA_zPt_corr.root"},
        'SIG': {'SIG': f"{datadir}DY_zPt_corr.root"},
        'BKG': {
            'WW': f"{datadir}WW_zPt_corr.root",
            'WZ': f"{datadir}WZ_zPt_corr.root",
            'ZZ': f"{datadir}ZZ_zPt_corr.root",
            'TT': f"{datadir}TT_zPt_corr.root",
        },
        'GEN': {'GEN': f"{datadir}GEN_zPt_corr.root"}
    }

    sf_path = 'data/scaleFactors/Run2/UL/2018/2018_Z/Efficiencies_muon_generalTracks_Z_Run2018_UL_'
    SFs = {
        'ID': ntuple.load_hist(sf_path+'ID.root', 'NUM_MediumID_DEN_TrackerMuons_abseta_pt'),
        'ISO': ntuple.load_hist(sf_path+'ISO.root', 'NUM_TightRelIso_DEN_MediumID_abseta_pt')
    }
    hdir = 'hists/'
    pdir = 'plots/'

    args = parse_args()
    ROOT.gROOT.SetBatch()

    if args.ntuples:
        # ntuple.make_ntuples(nanoAODs, datasets, datadir)
        os.makedirs(hdir, exist_ok=True)
        ntuple.hist_zpt(ntuples, pt_bins, hdir)
        ntuple.weight_zpt(ntuples, hdir, diffpt_bins, eta_bins, phi_bins) # weight also backgrounds and GEN?

    if args.scale:
        os.makedirs(hdir, exist_ok=True)
        step1.hist_oneOverpT(ntuples_zPt, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir)
        step1.get_scale_corrections(samples_to_plot, eta_bins, phi_bins, charge_bins, hdir)
        # step1.apply_scale_corrections(ntuples_zPt, hdir)
        # step1.hist_oneOverpT(ntuples_corr, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='_roccor')
        
    # if args.zmass:
    #     zmass.hist_zmass(ntuples_corr, eta_bins, phi_bins, mass_bins, hdir)
    #     zmass.fit_zmass(eta_bins, phi_bins, hdir)
    #     zmass.plot_zmass(eta_bins, phi_bins, hdir)

    if args.resolution:
        pull_bins=np.linspace(-5,5,100)
        r_bins = np.linspace(0.5, 1.5, 1000)
        abseta_bins=np.linspace(0, 2.4, 13)
        nl_bins=[6.5,8.5,9.5,10.5,11.5,12.5,13.5,17.5]
        pt_bins=[25,30,35,40,50,60,80,110,150,200]
        step2.get_res_correction(ntuples_zPt["GEN"]["GEN"], pull_bins, r_bins, abseta_bins, nl_bins, pt_bins, pdir, hdir, do_plot=True)
        step2.plot_closure(ntuples_zPtples["GEN"]["GEN"], hdir, pdir)

    if args.iterative:
        step3.iterative_correction(samples=ntuples_zPt, eta_bins=eta_bins, phi_bins=phi_bins, hdir=hdir, pdir=pdir)
        step3.plot_closure(samples=ntuples_zPt, hdir=hdir, pdir=pdir, eta_bins=eta_bins, phi_bins=phi_bins, iterationsteps=20)

    if args.residual:
        abseta_bins=np.linspace(0, 2.4, 13)

        step4.residual_correction(samples=ntuples_zPt, abseta_bins=abseta_bins, hdir=hdir, pdir=pdir)
        step4.perform_fits(ntuples_zPt, abseta_bins, hdir, pdir)
        step4.plot_closure(ntuples_zPt, hdir, pdir)