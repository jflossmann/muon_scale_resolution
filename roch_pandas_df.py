import numpy as np
import os
import time

from python.pd_oopt import oopt
from python.pd_resolution import resolution
from python.pd_iterative import iterative
from python.pd_residual import residual
from python.pd_plots import *


pt_bins = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 200]
eta_bins = [-2.40, -2.10, -1.85, -1.60, -1.20, -0.80, -0.40,  0.00,  0.40,  0.80,  1.20,  1.60,  1.85,  2.10,  2.40]
phi_bins = np.linspace(-3.2, 3.2, 17)
abseta_bins = np.linspace(0, 2.4, 13)
nT_bins = [6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5]
pt_bins = [20,25,30,35,37.5,40,42.5,45,47.5,50,52.5,55,60,65,75,90,120,150,200]
pull_bins = np.linspace(-4,4,80)

datadir = "/home/dguthmann/Documents/MA/Rochester/Daten/"
hdir = 'hists/'
os.makedirs(hdir, exist_ok=True)

pdir = 'plots/'

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


if __name__=='__main__':
    
    t0=time.time()

    oopt(ntuples = ntuples_zPt, 
         eta_bins = eta_bins, 
         phi_bins = phi_bins, 
         datadir = datadir, 
         hdir = hdir, 
         n_cores=16)
    
    t1=time.time()
    print(t1-t0)

    resolution(abseta_bins, 
               nT_bins = nT_bins, 
               pt_bins = pt_bins, 
               pull_bins = pull_bins, 
               datadir = datadir, 
               hdir = hdir, 
               pdir = pdir, 
               plot=True, 
               n_cores=16)
    
    t2=time.time()
    print(t2-t1)

    iterative(ntuples = ntuples_zPt, 
              eta_bins = eta_bins, 
              phi_bins = phi_bins, 
              mass_range = [86,96], 
              n_iterations = 15, 
              datadir = datadir, 
              hdir = hdir, 
              plotdir=pdir,
              n_cores=2)
    
    t3=time.time()
    print(t3-t2)

    residual(ntuples=ntuples_zPt,
             abseta_bins = abseta_bins,
             mass_range=[86,96],
             datadir=datadir,
             hdir=hdir,
             pdir=pdir,
             fit_bins=40,
             n_cores=8)

    t4=time.time()
    print(t4-t3)
   
    plot_lambda_fixed_eta(eta_bins, phi_bins, hdir, pdir)
    plot_lambda_fixed_phi(eta_bins, phi_bins, hdir, pdir)
    plot_kappa_fixed_eta(eta_bins, phi_bins, hdir, pdir)
    plot_kappa_fixed_phi(eta_bins, phi_bins, hdir, pdir)
