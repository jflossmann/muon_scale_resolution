import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
from python.utils import load_from_root



def oopt(ntuples, eta_bins, phi_bins, datadir, hdir, n_cores=8):
    
    df_gen = load_from_root(ntuples["GEN"]["GEN"])
    df_reco = load_from_root(ntuples["SIG"]["SIG"])

    eta_phi_bins=[]
    for eta_index in range(len(eta_bins)-1):
        for phi_index in range(len(phi_bins)-1):
            eta_phi_bins.append([eta_index, phi_index])

    #calc one over pt for gen and reco
    df_gen[f"oogenpt_1"] = 1/df_gen[f"genpt_1"]
    df_gen[f"oogenpt_2"] = 1/df_gen[f"genpt_2"]
    df_reco[f"oopt_1"]   = 1/df_reco[f"pt_1"]
    df_reco[f"oopt_2"]   = 1/df_reco[f"pt_2"]

    def one_over_pt_bin(bin):
        i,j = bin

        for k in [1,2]:
            #go to correct bin in dataframes, calc 1/pt
            eta_filter_gen = (eta_bins[i] < df_gen[f"geneta_{k}"]) & (df_gen[f"geneta_{k}"] <= eta_bins[i+1])
            phi_filter_gen = (phi_bins[j] < df_gen[f"genphi_{k}"]) & (df_gen[f"genphi_{k}"] <= phi_bins[j+1])
            df_gen_bin     = df_gen[eta_filter_gen & phi_filter_gen]

            eta_filter_reco = (eta_bins[i] < df_reco[f"eta_{k}"]) & (df_reco[f"eta_{k}"] <= eta_bins[i+1])
            phi_filter_reco = (phi_bins[j] < df_reco[f"phi_{k}"]) & (df_reco[f"phi_{k}"] <= phi_bins[j+1])
            df_reco_bin     = df_reco[eta_filter_reco & phi_filter_reco]

            if k==1:
                gm = np.mean(df_gen_bin[f"oogenpt_{k}"])
                xm = np.mean(df_reco_bin[f"oopt_{k}"])
            elif k==2:
                gp = np.mean(df_gen_bin[f"oogenpt_{k}"])
                xp = np.mean(df_reco_bin[f"oopt_{k}"])

            #calc kappa and lambda
                    
        kap = (gp+gm)/(xp+xm)
        lam = (gp*xm-gm*xp)/(xp+xm)

        return np.array([kap, lam])
    
    KapLam=np.array(Parallel(n_jobs=n_cores)(delayed(one_over_pt_bin)(bin) for bin in tqdm(eta_phi_bins))).T
            
    kap = np.reshape(KapLam[0],[len(eta_bins)-1,len(phi_bins)-1])
    lam = np.reshape(KapLam[1],[len(eta_bins)-1,len(phi_bins)-1])

    def apply_correction_bin_k(bin, k):
        i, j = bin
            
        #go to correct bin in dataframes, calc 1/pt
        eta_filter_gen = (eta_bins[i] < df_gen[f"geneta_{k}"]) & (df_gen[f"geneta_{k}"] <= eta_bins[i+1])
        phi_filter_gen = (phi_bins[j] < df_gen[f"genphi_{k}"]) & (df_gen[f"genphi_{k}"] <= phi_bins[j+1])
        df_gen_bin     = df_gen[eta_filter_gen & phi_filter_gen].copy()

        if k==1:
            df_gen_bin[f"pt_{k}_cor"] = 1 / (kap[i][j] / df_gen_bin[f"pt_{k}"] + lam[i][j])
        elif k==2:
            df_gen_bin[f"pt_{k}_cor"] = 1 / (kap[i][j] / df_gen_bin[f"pt_{k}"] - lam[i][j])
        
        return df_gen_bin
    
    df_gen = pd.concat(Parallel(n_jobs=n_cores)(delayed(apply_correction_bin_k)(bin, k=1) for bin in tqdm(eta_phi_bins)))  
    df_gen = pd.concat(Parallel(n_jobs=n_cores)(delayed(apply_correction_bin_k)(bin, k=2) for bin in tqdm(eta_phi_bins)))

    #safe gen-file as hdf
    print("save")
    df_gen.to_hdf(f"{datadir}df_gen.hdf5", key="df")
    np.savetxt(f"{hdir}kappa_step1.csv", kap, delimiter=",")
    np.savetxt(f"{hdir}lambda_step1.csv", lam, delimiter=",")

    return 
