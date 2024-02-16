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
    ROOT.gROOT.SetBatch()

    args = parse_args()

    # adjustable paths
    datadir = '/ceph/jdriesch/rochester/ntuples/'


    # definition of paths
    nanoAODs = 'data/nanoAODs.yaml'
    datasets = 'data/datasets.yaml'
    sf_path = 'data/scaleFactors/Run2/UL/2018/2018_Z/\
        Efficiencies_muon_generalTracks_Z_Run2018_UL_'
    SFs = {
        'id': [
            sf_path+'ID.root',
            'NUM_MediumID_DEN_TrackerMuons_abseta_pt'
        ],
        'iso': [
            sf_path+'ISO.root',
            'NUM_TightRelIso_DEN_MediumID_abseta_pt'
        ]
    }
    
    ntuples_raw = {
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

    ntuples = {
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

    hdir = 'hists/'
    pdir = 'plots/'
    os.makedirs(hdir, exist_ok=True)
    os.makedirs(pdir, exist_ok=True)


    # bin defintion
    pt_bins = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 200]
    oneOverPt_bins = np.linspace(0.019,0.027, 50)
    r_bins = np.linspace(0.5, 1.5, 1000)
    pull_bins=np.linspace(-5,5,100)

    eta_bins = [-2.4, -2.1, -1.85, -0.4, 0, 0.4, 1.85, 2.1, 2.4]
    abseta_bins = np.linspace(0, 2.4, 13)

    phi_bins = np.linspace(-3.2, 3.2, 17)

    charge_bins = [-2, 0, 2]

    mass_bins = np.linspace(75, 105, 61)

    nl_bins=[6.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 17.5]


    if args.ntuples:
        ntuple.make_ntuples(nanoAODs, datasets, datadir)
        ntuple.hist_zpt(ntuples_raw, pt_bins, hdir)
        ntuple.weight_zpt(ntuples_raw, hdir, eta_bins, phi_bins, SFs) # weight also backgrounds and GEN?

    if args.scale:
        step1.hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir)
        step1.get_scale_corrections(["DATA", "SIG", "GEN"], eta_bins, phi_bins, charge_bins, hdir)

    if args.resolution:
        step2.get_res_correction(ntuples["GEN"]["GEN"], pull_bins, r_bins, abseta_bins, nl_bins, pt_bins, pdir, hdir, do_plot=True)
        step2.plot_closure(ntuples["GEN"]["GEN"], hdir, pdir)

    if args.iterative:
        step3.iterative_correction(samples=ntuples, eta_bins=eta_bins, phi_bins=phi_bins, hdir=hdir, pdir=pdir)
        step3.plot_closure(samples=ntuples, hdir=hdir, pdir=pdir, eta_bins=eta_bins, phi_bins=phi_bins, iterationsteps=20)

    if args.residual:
        step4.residual_correction(samples=ntuples, abseta_bins=abseta_bins, hdir=hdir, pdir=pdir)
        step4.perform_fits(ntuples, abseta_bins, hdir, pdir)
        step4.plot_closure(ntuples, hdir, pdir)