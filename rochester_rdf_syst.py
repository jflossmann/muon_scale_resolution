import ROOT
import numpy as np
import matplotlib.pyplot as plt
from array import array
from argparse import ArgumentParser
import os
from time import time

# import modules
import python.ntuple as ntuple
import python.step1 as step1
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
    parser.add_argument(
        '-B',
        '--bootstrap',
        default=False,
        action='store_true',
        help='Activate bootstrapping'
    )

    parser.add_argument(
        "--process", 
        type=int,
        default=-1,
        help="Number of condor process"
    )
    parser.add_argument(
        "--syst", 
        type=str,
        default='',
        help="Evaluate systematic. Only one at a time!"
    )

    parser.add_argument(
        '-Y',
        '--year',
        type=str,
        default='',
        help='specific run'
    )
    parser.add_argument(
        '-P',
        '--plot',
        default=False,
        action='store_true',
        help='Save plots'
    )

    args = parser.parse_args()
    return args


if __name__=='__main__':
    t0 = time()
    ROOT.gROOT.SetBatch()

    args = parse_args()

    year = args.year

    # adjustable paths
    datadir = f'/ceph/jdriesch/rochester/{year}/'

    postEE = ''
    if year in ['2022E', '2022F', '2022G']:
        postEE = '_EE'

    # definition of paths
    nanoAODs = f'data/{year}/nanoAODs.yaml'
    datasets = f'data/{year}/datasets.yaml'
    sf_path = f'data/scaleFactors/Run3/2022{postEE}/2022_Z/'\
        f'ScaleFactors_Muon_Z_ID_ISO_2022{postEE}_schemaV2.json'
    SFs = {
        'id': [
            sf_path,
            'NUM_MediumID_DEN_TrackerMuons'
        ],
        'iso': [
            sf_path,
            'NUM_LoosePFIso_DEN_MediumID'
        ]
    }
    
    ntuples_raw = {
        'DATA': {'DATA': f"{datadir}ntuples/DATA_*.root"},
        'SIG': {'SIG': f"{datadir}ntuples/DY_*.root"},
        'BKG': {            
            'WW': f"{datadir}ntuples/WW_*.root",
            'WZ': f"{datadir}ntuples/WZ_*.root",
            'ZZ': f"{datadir}ntuples/ZZ_*.root",
            'TT': f"{datadir}ntuples/TT_*.root",
        },
        'GEN': {'GEN': f"{datadir}ntuples/GEN_*.root"}
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

    hdir = f'hists/{year}/'
    pdir = f'plots/{year}/'

    golden_json = 'data/jsons/Run3_2022_2023_Golden.json'

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
    m_bins_3 = np.linspace(86, 96, 100)
    m_bins_4 = np.linspace(80, 102, 44)
    m_bins_4 = np.linspace(85, 97, 44)

    nl_bins=[6.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 17.5]

    weight = "zPtWeight*genWeight/sumwWeight*xsec*sf_id*sf_iso"

    ROOT.gROOT.ProcessLine(f'gRandom->SetSeed({args.process});')

    # make sure that hists and plots are saved in different places for syst evaluation
    if args.process == -1 and args.syst=='':
        # this is the nominal case
        print("Nominal case.")
        # ROOT.EnableImplicitMT(5)

    elif '1' in args.syst and args.process > -1:
        oneOverPt_bins = np.linspace(
            0.022 - 0.001 * (args.process+1),
            0.024 - 0.001 * (args.process+1),
            40 + (args.process+1) * 5
        )
        hdir += f'condor/bin_step1/{args.process}/'
        pdir += f'condor/bin_step1/{args.process}/'

    elif '3' in args.syst and args.process > -1:
        m_bins_3 = np.linspace(
            91 - 2.5 - 0.5*args.process,
            91 + 2.5 + 0.5*args.process,
            50 + 10 * args.process
        )
        hdir += f'condor/bin_step3/{args.process}/'
        pdir += f'condor/bin_step3/{args.process}/'

    elif '4' in args.syst and args.process > -1:
        m_bins_4 = np.linspace(
            91 - 5 - args.process,
            91 + 5 + args.process,
            20 + 4 * args.process
        )
        hdir += f'condor/bin_step4/{args.process}/'
        pdir += f'condor/bin_step4/{args.process}/'

    elif 'zPt' in args.syst:
        weight = "genWeight/sumwWeight*xsec*sf_id*sf_iso"
        hdir += 'condor/no_zPt/'
        pdir += 'condor/no_zPt/'

    elif 'stat' in args.syst and args.process > -1:
        weight = "zPtWeight*genWeight/sumwWeight*xsec*sf_id*sf_iso*bs_weight"
        hdir += f'condor/{args.process}/'
        pdir += f'condor/{args.process}/'

    else:
        print("somethings broken. Please check input!")

    os.makedirs(hdir, exist_ok=True)
    os.makedirs(pdir, exist_ok=True)

 
    if args.ntuples:
        ntuple.make_ntuples(nanoAODs, datasets, datadir+'ntuples/', golden_json)
        ntuple.hist_zpt(ntuples_raw, pt_bins, hdir)
        ntuple.weight_zpt(ntuples_raw, hdir, eta_bins, phi_bins, SFs)

    if args.scale:
        step1.hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, weight)
        step1.get_scale_corrections(["DATA", "SIG", "GEN"], eta_bins, phi_bins, charge_bins, hdir)

    if args.resolution:
        step2.get_res_correction(ntuples["GEN"]["GEN"], pull_bins, r_bins, abseta_bins, nl_bins, pt_bins, pdir, hdir, True, weight)
        if args.plot:
            step2.plot_closure(ntuples["GEN"]["GEN"], hdir, pdir, weight)

    if args.iterative:
        itsteps = 25
        step3.iterative_correction(samples=ntuples, eta_bins=eta_bins, phi_bins=phi_bins, mass_bins=m_bins_3, hdir=hdir, pdir=pdir, iterationsteps=itsteps, weight=weight)
        if args.plot:
            step3.plot_closure(samples=ntuples, hdir=hdir, pdir=pdir, eta_bins=eta_bins, phi_bins=phi_bins, iterationsteps=itsteps)

    if args.residual:
        step4.residual_correction(samples=ntuples, abseta_bins=abseta_bins, hdir=hdir, pdir=pdir, weight=weight)
        step4.perform_fits(ntuples, abseta_bins, hdir, pdir, m_bins_4)
        if args.plot:
            step4.plot_closure(ntuples, hdir, pdir, weight, m_bins_4)

    print(f"Done in {round(time()-t0, 1)}s.")
