import pandas as pd
import numpy as np
import uproot
import matplotlib.pyplot as plt
from python.utils import load_from_root, Inv_mass, Minfinder, datahist_minus_bkghist
from joblib import Parallel, delayed
import os
from tqdm import tqdm
import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)

def chi2(k, z1, z2, df_gen, reco_hist, bins, rang):

    #reconstruct Z-Mass with smearing
    MZ=Inv_mass(pt1  = (1+k*z1)*df_gen["genpt_1_smeared"],
                pt2  = (1+k*z2)*df_gen["genpt_2_smeared"],
                eta1 = df_gen["geneta_1"],
                eta2 = df_gen["geneta_2"],
                phi1 = df_gen["genphi_1"],
                phi2 = df_gen["genphi_2"])
    
    #make histogram
    gen_hist=np.histogram(MZ, bins=bins, range=rang, density=True)
    #calc chi2
    x2=np.sum((gen_hist[0]-reco_hist[0])**2/reco_hist[0])
    return x2

def residual(ntuples, abseta_bins, mass_range, datadir, hdir, pdir, fit_bins=30, n_cores=8):

    pdir = pdir+'residual/k/'
    os.makedirs(pdir, exist_ok=True)

    #load gen, cut on mZ
    df_gen = pd.read_hdf(f"{datadir}df_gen.hdf5", key="df")
    mass_filter_gen = (mass_range[0]-2<df_gen["mass_Z_smeared"]) & (df_gen["mass_Z_smeared"]<mass_range[1]+2)
    df_gen_m = df_gen[mass_filter_gen]

    #load reco, cut on mZ, calc abseta
    df_reco = pd.read_hdf(f"{datadir}df_reco.hdf5", key="df")
    df_reco["abseta_1"] = np.abs(df_reco["eta_1"])
    df_reco["abseta_2"] = np.abs(df_reco["eta_2"])
    df_reco["weight"] = df_reco["zPtWeight"]*df_reco["genWeight"]/df_reco["sumwWeight"]*df_reco["xsec"]*df_reco["sf_id"]*df_reco["sf_iso"]
    mass_filter_reco = (mass_range[0]-2<df_reco["mass_Z_cor"]) & (df_reco["mass_Z_cor"]<mass_range[1]+2)
    df_reco_m=df_reco[mass_filter_reco]

    #load data, cut on mZ, calc abseta
    df_data = pd.read_hdf(f"{datadir}df_data.hdf5", key="df")
    df_data["abseta_1"] = np.abs(df_data["eta_1"])
    df_data["abseta_2"] = np.abs(df_data["eta_2"])
    mass_filter_data = (mass_range[0]-2<df_data["mass_Z_cor"]) & (df_data["mass_Z_cor"]<mass_range[1]+2)
    df_data_m = df_data[mass_filter_data]
    
    #load bkgs, cut on mZ, calc abseta
    dfs_bkg_m = []
    for subtyp in ntuples["BKG"]:
        df_bkg = load_from_root(ntuples["BKG"][subtyp])
        df_bkg["abseta_1"] = np.abs(df_bkg["eta_1"])
        df_bkg["abseta_2"] = np.abs(df_bkg["eta_2"])
        df_bkg["weight"] = df_bkg["zPtWeight"]*df_bkg["genWeight"]/df_bkg["sumwWeight"]*df_bkg["xsec"]*df_bkg["sf_id"]*df_bkg["sf_iso"]
        mass_filter_bkg = (mass_range[0]-2<df_bkg["mass_Z"]) & (df_bkg["mass_Z"]<mass_range[1]+2)
        dfs_bkg_m.append(df_bkg[mass_filter_bkg])

    def residual_bin(bin):
        i = bin

        abseta_filter_gen_1 = (abseta_bins[i] < df_gen_m[f"absgeneta_1"]) & (df_gen_m[f"absgeneta_1"] <= abseta_bins[i+1])
        abseta_filter_gen_2 = (abseta_bins[i] < df_gen_m[f"absgeneta_2"]) & (df_gen_m[f"absgeneta_2"] <= abseta_bins[i+1])
        df_gen_me = df_gen_m[abseta_filter_gen_1 & abseta_filter_gen_2]

        abseta_filter_reco_1 = (abseta_bins[i] < df_reco_m[f"abseta_1"]) & (df_reco_m[f"abseta_1"] <= abseta_bins[i+1])
        abseta_filter_reco_2 = (abseta_bins[i] < df_reco_m[f"abseta_2"]) & (df_reco_m[f"abseta_2"] <= abseta_bins[i+1])
        df_reco_me = df_reco_m[abseta_filter_reco_1 & abseta_filter_reco_2]

        abseta_filter_data_1 = (abseta_bins[i] < df_data_m[f"abseta_1"]) & (df_data_m[f"abseta_1"] <= abseta_bins[i+1])
        abseta_filter_data_2 = (abseta_bins[i] < df_data_m[f"abseta_2"]) & (df_data_m[f"abseta_2"] <= abseta_bins[i+1])
        df_data_me = df_data_m[abseta_filter_data_1 & abseta_filter_data_2]
        
        dfs_bkg_me=[]
        for df_bkg_m in dfs_bkg_m:
            abseta_filter_bkg_1 = (abseta_bins[i] < df_bkg_m[f"abseta_1"]) & (df_bkg_m[f"abseta_1"] <= abseta_bins[i+1])
            abseta_filter_bkg_2 = (abseta_bins[i] < df_bkg_m[f"abseta_2"]) & (df_bkg_m[f"abseta_2"] <= abseta_bins[i+1])
            dfs_bkg_me.append(df_bkg_m[abseta_filter_bkg_1 & abseta_filter_bkg_2])

        #perform fit for mc, plot result
            
        fig, (ax0,ax1) = plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
        hist_reco = list(ax0.hist(df_reco_me["mass_Z_cor"], bins=fit_bins, range=mass_range, density=True, color="k", alpha=0.5, label="reco"))
        n_gen=len(df_gen_me)
        #draw random numbers for smearing
        z1 = np.random.normal(0, 1, n_gen)
        z2 = np.random.normal(0, 1, n_gen)

        #calc initial k
        k_mc = Minfinder(fun=chi2, args=(z1, z2, df_gen_me, hist_reco, fit_bins, mass_range), bounds=[0,2], n=10, N=8)

        
        MZ = Inv_mass(pt1  = (1+k_mc*z1)*df_gen_me["genpt_1_smeared"],
                      pt2  = (1+k_mc*z2)*df_gen_me["genpt_2_smeared"],
                      eta1 = df_gen_me["geneta_1"],
                      eta2 = df_gen_me["geneta_2"],
                      phi1 = df_gen_me["genphi_1"],
                      phi2 = df_gen_me["genphi_2"])
        
        
        
        gen_hist_old = ax0.hist(df_gen_me.mass_Z_smeared, bins=fit_bins, range=mass_range, density=True, histtype="step", color="r", label="old")
        gen_hist_new = ax0.hist(MZ,                       bins=fit_bins, range=mass_range, density=True, histtype="step", color="b", label="new")

        x2_old=np.sum((gen_hist_old[0]-hist_reco[0])**2/hist_reco[0])
        x2_new=np.sum((gen_hist_new[0]-hist_reco[0])**2/hist_reco[0])

        ax0.text(mass_range[0]+0.1, 0.001, s=f"χ²old:{round(x2_old,4)}, χ²new:{round(x2_new,4)}")
                        
        ax0.set_ylabel("normalized event count")
        ax0.set_xlim(left=mass_range[0], right=mass_range[1])
        ax0.set_xticks([])
        ax0.legend()

        ax1.set_xlim(left=mass_range[0],right=mass_range[1])
        #ax1.set_ylim(bottom=0.9,top=1.1)
        x=np.linspace(mass_range[0],mass_range[1], fit_bins)
        ax1.plot(x ,gen_hist_old[0]/hist_reco[0], ".", c="r", label="old/RECO")
        ax1.plot(x ,gen_hist_new[0]/hist_reco[0], "_", c="b", label="new/RECO")
        ax1.set_xlabel("M_µµ (GeV)")
        ax1.set_ylabel("ratio")
        ax1.grid(True)
        ax1.legend()
            
        plt.savefig(f"{pdir}eta_{i}_mc_residual.png")
        plt.clf()

        #perform fit for data, plot result

        fig, (ax0,ax1) = plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
        hist_data=datahist_minus_bkghist(variable="mass_Z_cor", 
                                            df_data=df_data_me, 
                                            df_mc=df_reco_me, 
                                            variable_bkg="mass_Z", 
                                            dfs_bkg=dfs_bkg_me, 
                                            bins=fit_bins,
                                            Range = mass_range,
                                            normalized=True)
        
        binwidth=hist_data[1][1]-hist_data[1][0]
        ax0.bar(hist_data[1], hist_data[0], width=binwidth, color="k", alpha=0.5, label="DATA")
        k_data=Minfinder(fun=chi2, args=(z1, z2, df_gen_me, hist_data, fit_bins, mass_range), bounds=[0,2], n=10, N=8)
        
        MZ = Inv_mass(pt1  = (1+k_data*z1)*df_gen_me["genpt_1_smeared"],
                      pt2  = (1+k_data*z2)*df_gen_me["genpt_2_smeared"],
                      eta1 = df_gen_me["geneta_1"],
                      eta2 = df_gen_me["geneta_2"],
                      phi1 = df_gen_me["genphi_1"],
                      phi2 = df_gen_me["genphi_2"])
        
        x2_old=np.sum((gen_hist_old[0]-hist_data[0])**2/hist_data[0])
        x2_new=np.sum((gen_hist_new[0]-hist_data[0])**2/hist_data[0])

        ax0.text(mass_range[0]+0.1, 0.001, s=f"χ²old:{round(x2_old,4)}, χ²new:{round(x2_new,4)}")
            
        ax0.set_ylabel("normalized event count")
        ax0.set_xlim(left=mass_range[0],right=mass_range[1])
        ax0.set_xticks([])
        ax0.legend()

        ax1.set_xlim(left=mass_range[0], right=mass_range[1])
        ax1.set_ylim(bottom=0.8,top=1.2)
        ax1.plot(x, gen_hist_old[0]/hist_data[0],".", c="r", label="old/DATA")
        ax1.plot(x, gen_hist_new[0]/hist_data[0],"_", c="b", label="new/DATA")
        ax1.set_xlabel("M_µµ (GeV)")
        ax1.set_ylabel("ratio")
        ax1.grid(True)
        ax1.legend()
            
        plt.savefig(f"{pdir}eta_{i}_data_residual.png")
        plt.clf()
        
        return np.array([k_mc, k_data])

    Ks = np.array(Parallel(n_jobs=n_cores)(delayed(residual_bin)(bin) for bin in tqdm(range(len(abseta_bins)-1)))).T

    #for bin in range(len(abseta_bins)-1):
    #    residual_bin(bin)
   
    print(Ks)
    k_mc = Ks[0]
    k_data = Ks[1]
    np.savetxt(f"{hdir}k_residual_mc.csv", k_mc, delimiter=",")
    np.savetxt(f"{hdir}k_residual_data.csv", k_data, delimiter=",")

    return
