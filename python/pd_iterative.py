import numpy as np 
import pandas as pd
from python.utils import load_from_root, get_background_scale_factor, get_bkg_average_product, get_bkg_average_single, Inv_mass
from joblib import Parallel, delayed
from tqdm import tqdm
import os
import matplotlib.pyplot as plt

def iterative(ntuples, eta_bins, phi_bins, mass_range, n_iterations, datadir, hdir, plotdir, n_cores=2):
    #make plotdir
    pdir = plotdir+'iterative/V/'
    os.makedirs(pdir, exist_ok=True)
    #load gen
    df_gen = pd.read_hdf(f"{datadir}df_gen.hdf5", key="df")

    #cut gen on Mµµ
    mass_filter_gen = (mass_range[0]<df_gen["genmass_Z"]) & (df_gen["genmass_Z"]<mass_range[1])
    df_gen_m = df_gen[mass_filter_gen]

    #fill 2d list with eta-phi-binned gen-dfs
    df_gen_mep1=[]
    df_gen_mep2=[]
    for i in range(len(eta_bins)-1):
        df_gen_mep1.append([])
        df_gen_mep2.append([])

        eta_filter_gen_1 = (eta_bins[i]<df_gen_m["geneta_1"]) & (df_gen_m["geneta_1"]<eta_bins[i+1])
        eta_filter_gen_2 = (eta_bins[i]<df_gen_m["geneta_2"]) & (df_gen_m["geneta_2"]<eta_bins[i+1])

        df_gen_me1 = df_gen_m[eta_filter_gen_1]
        df_gen_me2 = df_gen_m[eta_filter_gen_2]

        for j in range(len(phi_bins)-1):
            phi_filter_gen_1 = (phi_bins[j]<df_gen_me1["genphi_1"]) & (df_gen_me1["genphi_1"]<phi_bins[j+1])
            phi_filter_gen_2 = (phi_bins[j]<df_gen_me2["genphi_2"]) & (df_gen_me2["genphi_2"]<phi_bins[j+1])
            df_gen_mep1[i].append(df_gen_me1[phi_filter_gen_1])
            df_gen_mep2[i].append(df_gen_me2[phi_filter_gen_2])


    def iterative_type(typ):

        
        

        #arrays for updated correction-parameters
        kap = np.ones([len(eta_bins)-1, len(phi_bins)-1])
        lam = np.zeros([len(eta_bins)-1, len(phi_bins)-1])

        #lists to monitor and plot optimized statistic
        V1=[]
        V2=[]

        if typ=="mc":

            #load reco, set inital values for the corrected quantities
            df_reco = load_from_root(ntuples["SIG"]["SIG"])
            df_reco["pt_1_cor"] = df_reco["pt_1"]
            df_reco["pt_2_cor"] = df_reco["pt_2"]
            
            #iterate for n_iteration times
            for n in tqdm(range(n_iterations)):
                V1.append([])
                V2.append([])
                #reconstruct and cut on Mµµ
                df_reco["mass_Z_cor"] = Inv_mass(pt1=df_reco["pt_1_cor"],
                                                 pt2=df_reco["pt_2_cor"],
                                                 eta1=df_reco["eta_1"],
                                                 eta2=df_reco["eta_2"],
                                                 phi1=df_reco["phi_1"],
                                                 phi2=df_reco["phi_2"])
                
                mass_filter_reco = (mass_range[0]<df_reco["mass_Z_cor"]) & (df_reco["mass_Z_cor"]<mass_range[1])
                df_reco_m = df_reco[mass_filter_reco]

                #bin reco in eta-phi
                for i in range(len(eta_bins)-1):
                    V1[n].append([])
                    V2[n].append([])
                    eta_filter_reco_1 = (eta_bins[i]<df_reco_m["eta_1"]) & (df_reco_m["eta_1"]<eta_bins[i+1])
                    eta_filter_reco_2 = (eta_bins[i]<df_reco_m["eta_2"]) & (df_reco_m["eta_2"]<eta_bins[i+1])
                    
                    df_reco_me1 = df_reco_m[eta_filter_reco_1]
                    df_reco_me2 = df_reco_m[eta_filter_reco_2]
                    
                    for j in range(len(phi_bins)-1):
                        
                        phi_filter_reco_1 = (phi_bins[j]<df_reco_me1["phi_1"]) & (df_reco_me1["phi_1"]<phi_bins[j+1])
                        phi_filter_reco_2 = (phi_bins[j]<df_reco_me2["phi_2"]) & (df_reco_me2["phi_2"]<phi_bins[j+1])
                        
                        df_reco_mep1 = df_reco_me1[phi_filter_reco_1]
                        df_reco_mep2 = df_reco_me2[phi_filter_reco_2]

                        #get statistics

                        mean_MZ_reco_1=np.mean(df_reco_mep1["mass_Z_cor"])
                        mean_MZ_reco_2=np.mean(df_reco_mep2["mass_Z_cor"])
                        
                        Vm = -2*(np.mean(df_gen_mep1[i][j]["genmass_Z"]) - mean_MZ_reco_1)
                        Mm = 2*mean_MZ_reco_1
                        Km = -np.mean(df_reco_mep1["pt_1_cor"]*df_reco_mep1["mass_Z_cor"])

                        Vp = -2*(np.mean(df_gen_mep2[i][j]["genmass_Z"]) - mean_MZ_reco_2)
                        Mp = 2*mean_MZ_reco_2
                        Kp = np.mean(df_reco_mep2["pt_2_cor"]*df_reco_mep2["mass_Z_cor"])
                        
                        V1[n][i].append(Vm)
                        V2[n][i].append(Vp)

                        #solve equation system

                        mu = (Vp/Kp-Vm/Km)/(Mp/Kp-Mm/Km)
                        kappa_new = 1+mu
                        lambda_new = ((Vp-Mp*mu)/Kp+(Vm-Mm*mu)/Km)/2

                        #update kap, lam
                        kap[i][j] = kap[i][j]*kappa_new
                        lam[i][j] = lam[i][j]*kappa_new+lambda_new

                        #update pt_cor in eta-phi-bin for complete df_reco 
                        eta_filter_reco_1 = (eta_bins[i]<df_reco["eta_1"]) & (df_reco["eta_1"]<eta_bins[i+1])
                        phi_filter_reco_1 = (phi_bins[j]<df_reco["phi_1"]) & (df_reco["phi_1"]<phi_bins[j+1]) 
                        df_reco.loc[eta_filter_reco_1 & phi_filter_reco_1, "pt_1_cor"] = 1/(kap[i][j]/df_reco.loc[eta_filter_reco_1 & phi_filter_reco_1, "pt_1"] - lam[i][j])

                        eta_filter_reco_2 = (eta_bins[i]<df_reco["eta_2"]) & (df_reco["eta_2"]<eta_bins[i+1])
                        phi_filter_reco_2 = (phi_bins[j]<df_reco["phi_2"]) & (df_reco["phi_2"]<phi_bins[j+1])
                        df_reco.loc[eta_filter_reco_2 & phi_filter_reco_2, "pt_2_cor"] = 1/(kap[i][j]/df_reco.loc[eta_filter_reco_2 & phi_filter_reco_2, "pt_2"] + lam[i][j])
                    
            #final m_Z
            df_reco["mass_Z_cor"] = Inv_mass(pt1=df_reco["pt_1_cor"],
                                             pt2=df_reco["pt_2_cor"],
                                             eta1=df_reco["eta_1"],
                                             eta2=df_reco["eta_2"],
                                             phi1=df_reco["phi_1"],
                                             phi2=df_reco["phi_2"])
            
            print(f"save {datadir}df_reco.hdf5")
            df_reco.to_hdf(f"{datadir}df_reco.hdf5", key="df")

        elif typ=="data":
            #load data, set inital values for the corrected quantities
            df_data = load_from_root(ntuples["DATA"]["DATA"])
            df_data["pt_1_cor"] = df_data["pt_1"]
            df_data["pt_2_cor"] = df_data["pt_2"]

            #load mc, go to mass-range
            df_reco = load_from_root(ntuples["SIG"]["SIG"])
            mass_filter_reco = (mass_range[0]<df_reco["mass_Z"]) & (df_reco["mass_Z"]<mass_range[1])
            df_reco["weight"]=df_reco["zPtWeight"]*df_reco["genWeight"]/df_reco["sumwWeight"]*df_reco["xsec"]*df_reco["sf_id"]*df_reco["sf_iso"]
            df_reco_m = df_reco[mass_filter_reco]

            #load background, go to mass-range
            dfs_bkg_m=[]
            for subtyp in ntuples["BKG"]:
                df_bkg = load_from_root(ntuples["BKG"][subtyp])
                df_bkg["weight"]=df_bkg["zPtWeight"]*df_bkg["genWeight"]/df_bkg["sumwWeight"]*df_bkg["xsec"]*df_bkg["sf_id"]*df_bkg["sf_iso"]
                mass_filter_bkg = (mass_range[0]<df_bkg["mass_Z"]) & (df_bkg["mass_Z"]<mass_range[1])
                dfs_bkg_m.append(df_bkg[mass_filter_bkg])

            #objects to store bkg information
                
            #event counts
            n_bkg1 = np.ones([len(eta_bins)-1, len(phi_bins)-1])
            n_bkg2 = np.ones([len(eta_bins)-1, len(phi_bins)-1])

            #means
            mean_MZ_1_bkg   = np.ones([len(eta_bins)-1, len(phi_bins)-1])
            mean_MZ_2_bkg   = np.ones([len(eta_bins)-1, len(phi_bins)-1])
            mean_ptMZ_1_bkg = np.ones([len(eta_bins)-1, len(phi_bins)-1])
            mean_ptMZ_2_bkg = np.ones([len(eta_bins)-1, len(phi_bins)-1])
            
            #iterate for n_iteration times
            for n in tqdm(range(n_iterations)):
                V1.append([])
                V2.append([])
                #reconstruct and cut on Mµµ
                df_data["mass_Z_cor"] = Inv_mass(pt1=df_data["pt_1_cor"],
                                                 pt2=df_data["pt_2_cor"],
                                                 eta1=df_data["eta_1"],
                                                 eta2=df_data["eta_2"],
                                                 phi1=df_data["phi_1"],
                                                 phi2=df_data["phi_2"])
                
                mass_filter_data = (mass_range[0]<df_data["mass_Z_cor"]) & (df_data["mass_Z_cor"]<mass_range[1])
                df_data_m = df_data[mass_filter_data]

                #bin data in eta-phi
                for i in range(len(eta_bins)-1):
                    V1[n].append([])
                    V2[n].append([])
                    eta_filter_data_1 = (eta_bins[i]<df_data_m["eta_1"]) & (df_data_m["eta_1"]<eta_bins[i+1])
                    eta_filter_data_2 = (eta_bins[i]<df_data_m["eta_2"]) & (df_data_m["eta_2"]<eta_bins[i+1])
                    
                    df_data_me1 = df_data_m[eta_filter_data_1]
                    df_data_me2 = df_data_m[eta_filter_data_2]

                    if n==0:
                        eta_filter_reco_1 = (eta_bins[i] < df_reco_m["eta_1"]) & (df_reco_m["eta_1"] <= eta_bins[i+1])
                        eta_filter_reco_2 = (eta_bins[i] < df_reco_m["eta_2"]) & (df_reco_m["eta_2"] <= eta_bins[i+1])
                        df_reco_me1 = df_reco_m[eta_filter_reco_1]
                        df_reco_me2 = df_reco_m[eta_filter_reco_2]

                        dfs_bkg_me1=[]
                        dfs_bkg_me2=[]
                        for df_bkg_m in dfs_bkg_m:
                            eta_filter_bkg_1 = (eta_bins[i] < df_bkg_m["eta_1"]) & (df_bkg_m["eta_1"] <= eta_bins[i+1])
                            eta_filter_bkg_2 = (eta_bins[i] < df_bkg_m["eta_2"]) & (df_bkg_m["eta_2"] <= eta_bins[i+1])
                            dfs_bkg_me1.append(df_bkg_m[eta_filter_bkg_1])
                            dfs_bkg_me2.append(df_bkg_m[eta_filter_bkg_2])

                    for j in range(len(phi_bins)-1):
                        
                        phi_filter_data_1 = (phi_bins[j]<df_data_me1["phi_1"]) & (df_data_me1["phi_1"]<phi_bins[j+1])
                        phi_filter_data_2 = (phi_bins[j]<df_data_me2["phi_2"]) & (df_data_me2["phi_2"]<phi_bins[j+1])
                        
                        df_data_mep1 = df_data_me1[phi_filter_data_1]
                        df_data_mep2 = df_data_me2[phi_filter_data_2]

                        if n==0:
                            #get background scale-factors, calc background-means
                            phi_filter_reco_1 = (phi_bins[j] < df_reco_me1["phi_1"]) & (df_reco_me1["phi_1"] <= phi_bins[j+1])
                            phi_filter_reco_2 = (phi_bins[j] < df_reco_me2["phi_2"]) & (df_reco_me2["phi_2"] <= phi_bins[j+1])
                            df_reco_mep1 = df_reco_me1[phi_filter_reco_1]
                            df_reco_mep2 = df_reco_me2[phi_filter_reco_2]

                            dfs_bkg_mep1 = []
                            dfs_bkg_mep2 = []
                            for df_bkg_me1 in dfs_bkg_me1:
                                phi_filter_bkg_1 = (phi_bins[j] < df_bkg_me1["phi_1"]) & (df_bkg_me1["phi_1"] <= phi_bins[j+1])
                                dfs_bkg_mep1.append(df_bkg_me1[phi_filter_bkg_1])
                            for df_bkg_me2 in dfs_bkg_me2:
                                phi_filter_bkg_2 = (phi_bins[j] < df_bkg_me2["phi_2"]) & (df_bkg_me2["phi_2"] <= phi_bins[j+1])
                                dfs_bkg_mep2.append(df_bkg_me2[phi_filter_bkg_2])

                            #calc statistics and event counts for bkg
                            sf1 = get_background_scale_factor(df_data_mep1, df_reco_mep1, dfs_bkg_mep1, 100, mass_range)
                            sf2 = get_background_scale_factor(df_data_mep1, df_reco_mep1, dfs_bkg_mep1, 100, mass_range)
                            mean_MZ_1_bkg[i][j], n_bkg1[i][j] = get_bkg_average_single(dfs_bkg_mep1, "mass_Z", 100, mass_range, sf1)
                            mean_MZ_2_bkg[i][j], n_bkg2[i][j] = get_bkg_average_single(dfs_bkg_mep2, "mass_Z", 100, mass_range, sf2)
                            mean_ptMZ_1_bkg[i][j]             = get_bkg_average_product(dfs_bkg_mep1, "mass_Z", "pt_1", 1000, [mass_range[0]*20, mass_range[1]*200])
                            mean_ptMZ_2_bkg[i][j]             = get_bkg_average_product(dfs_bkg_mep2, "mass_Z", "pt_2", 1000, [mass_range[0]*20, mass_range[1]*200])

                        #get statistics
                        
                        #data event counts
                        n_data1 = len(df_data_mep1)
                        n_data2 = len(df_data_mep2)

                        #scaled means
                        mean_MZ_data_sig_1   = (
                            np.mean(df_data_mep1["mass_Z_cor"])*n_data1 - mean_MZ_1_bkg[i][j]*n_bkg1[i][j])/(n_data1-n_bkg1[i][j])
                        mean_MZ_data_sig_2   = (
                            np.mean(df_data_mep2["mass_Z_cor"])*n_data2 - mean_MZ_2_bkg[i][j]*n_bkg2[i][j])/(n_data2-n_bkg2[i][j])
                        mean_MZpt_data_sig_1 = (
                            np.mean(df_data_mep1["mass_Z_cor"]*df_data_mep1["pt_1_cor"])*n_data1 - mean_ptMZ_1_bkg[i][j]*n_bkg1[i][j])/(n_data1-n_bkg1[i][j])
                        mean_MZpt_data_sig_2 = (
                            np.mean(df_data_mep2["mass_Z_cor"]*df_data_mep2["pt_2_cor"])*n_data2 - mean_ptMZ_2_bkg[i][j]*n_bkg2[i][j])/(n_data2-n_bkg2[i][j])
                        
                        Vm = -2*(np.mean(df_gen_mep1[i][j]["genmass_Z"]) - mean_MZ_data_sig_1)
                        Mm = 2*mean_MZ_data_sig_1
                        Km = -mean_MZpt_data_sig_1

                        Vp = -2*(np.mean(df_gen_mep2[i][j]["genmass_Z"]) - mean_MZ_data_sig_2)
                        Mp = 2*mean_MZ_data_sig_2
                        Kp = mean_MZpt_data_sig_2
                        
                        V1[n][i].append(Vm)
                        V2[n][i].append(Vp)

                        #solve equation system

                        mu = (Vp/Kp-Vm/Km)/(Mp/Kp-Mm/Km)
                        kappa_new = 1+mu
                        lambda_new = ((Vp-Mp*mu)/Kp+(Vm-Mm*mu)/Km)/2

                        #update kap, lam
                        kap[i][j] = kap[i][j]*kappa_new
                        lam[i][j] = lam[i][j]*kappa_new+lambda_new

                        #update pt_cor in eta-phi-bin for complete df_data 
                        eta_filter_data_1 = (eta_bins[i]<df_data["eta_1"]) & (df_data["eta_1"]<eta_bins[i+1])
                        phi_filter_data_1 = (phi_bins[j]<df_data["phi_1"]) & (df_data["phi_1"]<phi_bins[j+1]) 
                        df_data.loc[eta_filter_data_1 & phi_filter_data_1, "pt_1_cor"] = 1/(
                            kap[i][j]/df_data.loc[eta_filter_data_1 & phi_filter_data_1, "pt_1"] - lam[i][j])

                        eta_filter_data_2 = (eta_bins[i]<df_data["eta_2"]) & (df_data["eta_2"]<eta_bins[i+1])
                        phi_filter_data_2 = (phi_bins[j]<df_data["phi_2"]) & (df_data["phi_2"]<phi_bins[j+1])
                        df_data.loc[eta_filter_data_2 & phi_filter_data_2, "pt_2_cor"] = 1/(
                            kap[i][j]/df_data.loc[eta_filter_data_2 & phi_filter_data_2, "pt_2"] + lam[i][j])
                    
            #final m_Z
            df_data["mass_Z_cor"] = Inv_mass(pt1=df_data["pt_1_cor"],
                                             pt2=df_data["pt_2_cor"],
                                             eta1=df_data["eta_1"],
                                             eta2=df_data["eta_2"],
                                             phi1=df_data["phi_1"],
                                             phi2=df_data["phi_2"])
            
            print(f"save {datadir}df_data.hdf5")
            df_data.to_hdf(f"{datadir}df_data.hdf5", key="df")
        

        #plot the develpoment of V in every bin
        for i in range(len(eta_bins)-1):
            for j in range(len(phi_bins)-1):
                    for n in range(n_iterations):
                        plt.plot(n, V1[n][i][j],"o",c="r")
                        plt.plot(n, V2[n][i][j],"o",c="b")
                    plt.savefig(f"{pdir}{typ}eta{i}_phi{j}_V_iterations")
                    plt.clf()

        return kap, lam

    KapLam = Parallel(n_jobs=n_cores)(delayed(iterative_type)(typ) for typ in ["mc", "data"])
    
    np.savetxt(f"{hdir}kappa_mc.csv", KapLam[0][0], delimiter=",")
    np.savetxt(f"{hdir}lambda_mc.csv", KapLam[0][1], delimiter=",")
    np.savetxt(f"{hdir}kappa_data.csv", KapLam[1][0], delimiter=",")
    np.savetxt(f"{hdir}lambda_data.csv", KapLam[1][1], delimiter=",")
