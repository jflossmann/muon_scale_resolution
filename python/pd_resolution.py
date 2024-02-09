import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
pd.options.mode.chained_assignment = None
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
from python.utils import load_from_root, CrystalBall, Poly2, Inv_mass
import os

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)

def resolution(abseta_bins, nT_bins, pt_bins, pull_bins, datadir, hdir, pdir, plot=True, n_cores=8):
    
    df_gen = pd.read_hdf(f"{datadir}df_gen.hdf5", key="df")

    abseta_nT_bins=[]
    for abseta_index in range(len(abseta_bins)-1):
        for nT_index in range(len(nT_bins)-1):
            abseta_nT_bins.append([abseta_index, nT_index])

    #calc abseta, R
    df_gen["absgeneta_1"] = np.abs(df_gen.geneta_1)
    df_gen["absgeneta_2"] = np.abs(df_gen.geneta_2)
    df_gen["R_1"]         = df_gen["genpt_1"]/df_gen["pt_1_cor"]
    df_gen["R_2"]         = df_gen["genpt_2"]/df_gen["pt_2_cor"]
    df_gen["R_1_std"]     = np.std(df_gen["R_1"])
    df_gen["R_2_std"]     = np.std(df_gen["R_2"])

    def Fits_bin(bin):
        """nested function, performs poly and CB-fit for given bin"""
        i, j = bin
        
        if plot:
            pdir_CB = pdir+'resolution/CB/'
            os.makedirs(pdir_CB, exist_ok=True)
            pdir_Poly = pdir+'resolution/Poly/'
            os.makedirs(pdir_Poly, exist_ok=True)
       
        #go to correct bin in dataframes, calc 1/pt
        abseta_filter_gen_1 = (abseta_bins[i] < df_gen[f"absgeneta_1"]) & (df_gen[f"absgeneta_1"] <= abseta_bins[i+1])
        nT_filter_gen_1     = (nT_bins[j] < df_gen[f"nTrkLayers_1"]) & (df_gen[f"nTrkLayers_1"] <= nT_bins[j+1])
        df_gen_bin_1        = df_gen[abseta_filter_gen_1 & nT_filter_gen_1]

        abseta_filter_gen_2 = (abseta_bins[i] < df_gen[f"absgeneta_2"]) & (df_gen[f"absgeneta_2"] <= abseta_bins[i+1])
        nT_filter_gen_2     = (nT_bins[j] < df_gen[f"nTrkLayers_2"]) & (df_gen[f"nTrkLayers_2"] <= nT_bins[j+1])
        df_gen_bin_2        = df_gen[abseta_filter_gen_2 & nT_filter_gen_2]

        R=pd.concat([df_gen_bin_1["R_1"], df_gen_bin_2["R_2"]])

        R_std_binned=[]
        for k in range(len(pt_bins)-1):
            pt_filter_gen_1 = (pt_bins[k] < df_gen_bin_1[f"genpt_1"]) & (df_gen_bin_1[f"genpt_1"] <= pt_bins[k+1])
            df_gen_bin_1_pt = df_gen_bin_1[pt_filter_gen_1]

            pt_filter_gen_2 = (pt_bins[k] < df_gen_bin_2[f"genpt_2"]) & (df_gen_bin_2[f"genpt_2"] <= pt_bins[k+1])
            df_gen_bin_2_pt =  df_gen_bin_2[pt_filter_gen_2] 

            R_std_binned.append(np.std(pd.concat([df_gen_bin_1_pt["R_1"], df_gen_bin_2_pt["R_2"]])))

        #calc needed quantities for all muons in this bin
        R_std = np.array(R_std_binned)
        pull  = (R-np.mean(R))/np.std(R)

        pt_bincenters=[]
        for l in range(len(pt_bins)-1):
            pt_bincenters.append(pt_bins[l]+(pt_bins[l+1]-pt_bins[l])/2)

        pt_bincenters=np.array(pt_bincenters)
        
        #find and replace nan values in fit-arrays
        ind_nan_poly=np.argwhere(np.isnan(R_std))
        R_std=np.delete(R_std, ind_nan_poly)
        pt_bincenters_left=np.delete(pt_bincenters, ind_nan_poly)

        if len(R_std)<3:
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        
        popt, pcov = curve_fit(Poly2, pt_bincenters_left, R_std)
        poly_a, poly_b, poly_c = popt
        
        if plot:
            plt.plot(pt_bincenters_left, R_std, "_", c="k",label="Data")
            x=np.linspace(min(pt_bins),max(pt_bins),200)
            plt.plot(x,Poly2(x,*popt), c="b", label="Fit")
            plt.ylabel("std(R)")
            plt.xlabel("genpt (GeV)")
            plt.title(f"η_{i}, nT_{j}")
            plt.legend()
            plt.savefig(f"{pdir_Poly}abseta{i}_nT{j}.png")
            plt.clf()

        #make pull histogram and fit
        pull_hist       = np.histogram(pull, bins=pull_bins, density=True)[0]
        pull_bincenters = pull_bins[:-1] + (pull_bins[1]-pull_bins[0])/2

        popt,pcov = curve_fit(CrystalBall, pull_bincenters, pull_hist, p0=np.array([0, 0.7, 1, 4]))
        CB_mu, CB_sig, CB_alpha, CB_n = popt
        
        if plot:
            #make CB-fit-plot
            plt.plot(pull_bincenters,pull_hist,"_", c="k",label="Hist")
            x=np.linspace(-4,4,200)
            plt.plot(x,CrystalBall(x,*popt), c="b", label="Fit")
            plt.ylabel("normalized event count")
            plt.xlabel("Pull(R)")
            plt.title(f"η_{i}, nT_{j}")
            plt.legend()
            plt.savefig(f"{pdir_CB}abseta{i}_nT{j}.png")
            plt.clf()

        return poly_a, poly_b, poly_c, CB_mu, CB_sig, CB_alpha, CB_n
    
    print("perform fits")
    
    FitPar=np.array(Parallel(n_jobs=n_cores)(delayed(Fits_bin)(bin) for bin in tqdm(abseta_nT_bins))).T
        
    poly_a   = np.reshape(FitPar[0], [len(abseta_bins)-1, len(nT_bins)-1])
    poly_b   = np.reshape(FitPar[1], [len(abseta_bins)-1, len(nT_bins)-1])
    poly_c   = np.reshape(FitPar[2], [len(abseta_bins)-1, len(nT_bins)-1])
    CB_mu    = np.reshape(FitPar[3], [len(abseta_bins)-1, len(nT_bins)-1])
    CB_sig   = np.reshape(FitPar[4], [len(abseta_bins)-1, len(nT_bins)-1])
    CB_alpha = np.reshape(FitPar[5], [len(abseta_bins)-1, len(nT_bins)-1])
    CB_n     = np.reshape(FitPar[6], [len(abseta_bins)-1, len(nT_bins)-1])

    def apply_smearing_bin_k(bin, k):
        """nested function, can be executed in paralell, returns smeared gen-pt of df_gen for given bin"""
        i, j = bin
        #go to correct bin in dataframes, calc 1/pt
        abseta_filter_gen = (abseta_bins[i] < df_gen[f"absgeneta_{k}"]) & (df_gen[f"absgeneta_{k}"] <= abseta_bins[i+1])
        nT_filter_gen = (nT_bins[j] < df_gen[f"nTrkLayers_{k}"]) & (df_gen[f"nTrkLayers_{k}"] <= nT_bins[j+1])
        df_gen_bin     = df_gen[abseta_filter_gen & nT_filter_gen].copy()


        pt=np.array(df_gen_bin[f"genpt_{k}"])
        std_poly = pt**2*poly_a[i][j] + pt*poly_b[i][j] + poly_c[i][j]
        std_poly[std_poly<0] = 0
        std_CB   = CB_sig[i][j]

        df_gen_bin[f"genpt_{k}_smeared"] = df_gen_bin[f"genpt_{k}"]*np.random.normal(1, std_poly*std_CB)
        
        return df_gen_bin
    
    print("smear muon 1")
    df_gen = pd.concat(Parallel(n_jobs=n_cores)(delayed(apply_smearing_bin_k)(bin, k=1) for bin in tqdm(abseta_nT_bins)))
    print("smear muon 2")
    df_gen = pd.concat(Parallel(n_jobs=n_cores)(delayed(apply_smearing_bin_k)(bin, k=2) for bin in tqdm(abseta_nT_bins)))

    df_gen["mass_Z_smeared"] = Inv_mass(df_gen["genpt_1_smeared"], 
                                                df_gen["genpt_2_smeared"],
                                                df_gen["geneta_1"],
                                                df_gen["geneta_2"],
                                                df_gen["genphi_1"],
                                                df_gen["genphi_2"])
    
    print(f"save {datadir}df_gen.hdf5")
    df_gen.to_hdf(f"{datadir}df_gen.hdf5", key="df")
    np.savetxt(f"{hdir}poly_a.csv", poly_a, delimiter=",")
    np.savetxt(f"{hdir}poly_b.csv", poly_b, delimiter=",")
    np.savetxt(f"{hdir}poly_c.csv", poly_c, delimiter=",")
    np.savetxt(f"{hdir}CB_mu.csv", CB_mu, delimiter=",")
    np.savetxt(f"{hdir}CB_sig.csv", CB_sig, delimiter=",")
    np.savetxt(f"{hdir}CB_alpha.csv", CB_alpha, delimiter=",")
    np.savetxt(f"{hdir}CB_n.csv", CB_n, delimiter=",")
    return
