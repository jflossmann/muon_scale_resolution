import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import uproot
from tqdm import tqdm
import json
import os

def minfinder(fun, args, bounds, n=10, N=10):
    i=0
    while i<N:
        X2=[]
        K=[]
        width=(bounds[1]-bounds[0])/n
        for k in np.linspace(bounds[0],bounds[1],n+1):
            X2.append(fun(k,*args))
            K.append(k)
        m=np.argmin(X2)

        if m==0:
            bounds=[bounds[0],K[m]+width]
        elif m==n:
            bounds=[K[m]-width,bounds[1]]
        else:
            bounds=[K[m]-width,K[m]+width]
        
        i+=1
        print(m, K[m], X2[m])
    return K[m]

def chi2(k, z1, z2, df_gen, reco_hist, bins=50, rang=[81,101]):

    #reconstruct Z-Mass with smearing
    MZ=np.sqrt( 2*(1+k*z1)*(1+k*z2)*df_gen.genpt_1_smeared*df_gen.genpt_2_smeared*(np.cosh(df_gen.geneta_1-df_gen.geneta_2)-np.cos(df_gen.genphi_1-df_gen.genphi_2)) )
    
    #make histogram
    gen_hist=np.histogram(MZ, bins=bins, range=rang, density=True)

    #calc chi2/me/mse
    #x2
    x2=np.sum((gen_hist[0]-reco_hist[0])**2/reco_hist[0])
    #me
    #me=np.sum(np.abs(gen_hist[0]-reco_hist[0]))/bins
    #mse
    #me=np.sum((gen_hist[0]-reco_hist[0])**2)/bins
    
    return x2

def residual_correction(samples, abseta_bins, hdir, pdir):
    pdir = pdir+'residual/resid/'
    os.makedirs(pdir, exist_ok=True)

    #dict for results
    k_results={}
    bins=40
    rang = [86,96]
    
    n=8
    #get gen data
    file=uproot.open(samples["GEN"]["GEN"])
    tree=file["Events"]
    variables=tree.keys()
    df_GEN=tree.arrays(variables, library="pd")
    df_GEN["genabseta_1"]=np.abs(df_GEN.geneta_1)
    df_GEN["genabseta_2"]=np.abs(df_GEN.geneta_2)

    #cut events with 81<M_z_smeared<101
    mask_M_Z=(df_GEN.genmass_Z>rang[0]-5) & (df_GEN.genmass_Z<rang[1]+5)
    df_GEN_m=df_GEN[mask_M_Z]

    for typ in samples:
        
        if typ == "SIG":
            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")
                
                #read data into dataframe
                file=uproot.open(samples[typ][subtyp])
                tree=file["Events"]
                variables=tree.keys()
                df_RECO=tree.arrays(variables, library="pd")
                df_RECO["abseta_1"]=np.abs(df_RECO.eta_1)
                df_RECO["abseta_2"]=np.abs(df_RECO.eta_2)

                mass_Z_filter_reco=(rang[0]<df_RECO["MASS_Z_COR"]) & (df_RECO["MASS_Z_COR"]<rang[1])

                df_RECO_m=df_RECO[mass_Z_filter_reco]

                #list for k(eta)
                K=[]
                for i in tqdm(range(len(abseta_bins)-1)):
                    
                    for j in range(i+1):
                        print(i,j)
                        if i==j:
                            #bin reco
                            abseta_1_filter_reco=(abseta_bins[i]<df_RECO["abseta_1"]) & (df_RECO["abseta_1"]<=abseta_bins[i+1])
                            abseta_2_filter_reco=(abseta_bins[i]<df_RECO["abseta_2"]) & (df_RECO["abseta_2"]<=abseta_bins[i+1])
                            df_RECO_me=df_RECO_m[abseta_1_filter_reco & abseta_2_filter_reco]

                            #bin gen
                            abseta_1_filter_gen=(abseta_bins[i]<df_GEN_m["genabseta_1"]) & (df_GEN_m["genabseta_1"]<=abseta_bins[i+1])
                            abseta_2_filter_gen=(abseta_bins[i]<df_GEN_m["genabseta_2"]) & (df_GEN_m["genabseta_2"]<=abseta_bins[i+1])
                            df_GEN_me=df_GEN_m[abseta_1_filter_gen & abseta_2_filter_gen]
                        else:
                            abseta_1_filter_reco_i=(abseta_bins[i]<df_RECO["abseta_1"]) & (df_RECO["abseta_1"]<=abseta_bins[i+1])
                            abseta_2_filter_reco_j=(abseta_bins[j]<df_RECO["abseta_2"]) & (df_RECO["abseta_2"]<=abseta_bins[j+1])
                            df_RECO_me_ij=df_RECO_m[abseta_1_filter_reco_i & abseta_2_filter_reco_j]
                            abseta_1_filter_reco_j=(abseta_bins[j]<df_RECO["abseta_1"]) & (df_RECO["abseta_1"]<=abseta_bins[j+1])
                            abseta_2_filter_reco_i=(abseta_bins[i]<df_RECO["abseta_2"]) & (df_RECO["abseta_2"]<=abseta_bins[i+1])
                            df_RECO_me_ji=df_RECO_m[abseta_1_filter_reco_j & abseta_2_filter_reco_i]

                            df_RECO_me=pd.concat([df_RECO_me_ij,df_RECO_me_ji])

                            #bin gen
                            abseta_1_filter_gen_i=(abseta_bins[i]<df_GEN_m["genabseta_1"]) & (df_GEN_m["genabseta_1"]<=abseta_bins[i+1])
                            abseta_2_filter_gen_j=(abseta_bins[j]<df_GEN_m["genabseta_2"]) & (df_GEN_m["genabseta_2"]<=abseta_bins[j+1])
                            df_GEN_me_ij=df_GEN_m[abseta_1_filter_gen_i & abseta_2_filter_gen_j]
                            abseta_1_filter_gen_j=(abseta_bins[j]<df_GEN_m["genabseta_1"]) & (df_GEN_m["genabseta_1"]<=abseta_bins[j+1])
                            abseta_2_filter_gen_i=(abseta_bins[i]<df_GEN_m["genabseta_2"]) & (df_GEN_m["genabseta_2"]<=abseta_bins[i+1])
                            df_GEN_me_ji=df_GEN_m[abseta_1_filter_gen_j & abseta_2_filter_gen_i]

                            df_GEN_me=pd.concat([df_GEN_me_ij,df_GEN_me_ji])

                        fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
                        hist_reco=list(ax0.hist(df_RECO_me["MASS_Z_COR"],bins=bins, density=True,  range=rang, color="k", alpha=0.5, label="reco"))

                        z1 = np.random.normal(0,1,len(df_GEN_me))
                        z2 = np.random.normal(0,1,len(df_GEN_me))

                        k=minfinder(fun=chi2, args=(z1, z2, df_GEN_me, hist_reco, bins, rang), bounds=[0,0.2], n=10, N=n)
                        print(k)
                        K.append(k)
                        MZ = np.sqrt( 2*(1+k*z1)*(1+k*z2)*df_GEN_me.genpt_1_smeared*df_GEN_me.genpt_2_smeared*(np.cosh(df_GEN_me.geneta_1-df_GEN_me.geneta_2)-np.cos(df_GEN_me.genphi_1-df_GEN_me.genphi_2)) )


                        gen_hist_old = ax0.hist(df_GEN_me.mass_Z_smeared, bins=bins, range=rang, density=True, histtype="step", color="r")
                        gen_hist_new = ax0.hist(MZ,                       bins=bins, range=rang, density=True, histtype="step", color="b")
                            
                        x2_old=np.sum((gen_hist_old[0]-hist_reco[0])**2/hist_reco[0])
                        x2_new=np.sum((gen_hist_new[0]-hist_reco[0])**2/hist_reco[0])
    
                        ax0.text(85, 0.001, s=f"χ²old:{round(x2_old,4)}, χ²new:{round(x2_new,4)}")
                            
                        ax0.set_ylabel("normalized event count")
                        ax0.set_xlim(left=rang[0],right=rang[1])
                        ax0.set_xticks([])
                        ax0.legend()

                        ax1.set_xlim(left=rang[0],right=rang[1])
                        ax1.plot(np.linspace(rang[0],rang[1],bins),hist_reco[0]/gen_hist_old[0],".", c="r", label="RECO/old")
                        ax1.plot(np.linspace(rang[0],rang[1],bins),hist_reco[0]/gen_hist_new[0],"_", c="b", label="RECO/new")
                        ax1.set_xlabel("M_µµ (GeV)")
                        ax1.set_ylabel("ratio")
                        ax1.grid(True)
                        ax1.legend()
                            
                        plt.savefig(f"{pdir}eta_{i}_{j}_{subtyp}_residual.png")
                        plt.clf()

                        #write k to all muons in this eta bin
                        if i==j:
                            
                            abseta_1_filter_gen=(abseta_bins[i]<df_GEN["genabseta_1"]) & (df_GEN["genabseta_1"]<=abseta_bins[i+1])
                            abseta_2_filter_gen=(abseta_bins[i]<df_GEN["genabseta_2"]) & (df_GEN["genabseta_2"]<=abseta_bins[i+1])

                            df_GEN.loc[abseta_1_filter_gen, "k1"]=k
                            df_GEN.loc[abseta_2_filter_gen, "k2"]=k

                        else:
                            abseta_1_filter_gen_i=(abseta_bins[i]<df_GEN["genabseta_1"]) & (df_GEN["genabseta_1"]<=abseta_bins[i+1])
                            abseta_2_filter_gen_j=(abseta_bins[j]<df_GEN["genabseta_2"]) & (df_GEN["genabseta_2"]<=abseta_bins[j+1])
                            
                            abseta_1_filter_gen_j=(abseta_bins[j]<df_GEN["genabseta_1"]) & (df_GEN["genabseta_1"]<=abseta_bins[j+1])
                            abseta_2_filter_gen_i=(abseta_bins[i]<df_GEN["genabseta_2"]) & (df_GEN["genabseta_2"]<=abseta_bins[i+1])
                        

                            df_GEN.loc[abseta_1_filter_gen_i & abseta_2_filter_gen_j, "k"]=k
                            df_GEN.loc[abseta_2_filter_gen_j & abseta_2_filter_gen_i, "k"]=k

                #safe ks in k_results-dict
                k_results[typ]=K
                rang=[81,101]
                bins=50
                #plot MZ for all eta bins combined
                fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
                
                #smear all pt according to respective k
                
                #calc random numbers
                z1=np.random.normal(0,1,len(df_GEN))
                z2=np.random.normal(0,1,len(df_GEN))

                MZ=np.sqrt( 2*(1+df_GEN.k*z1)*df_GEN.genpt_1_smeared*
                           (1+df_GEN.k*z2)*df_GEN.genpt_2_smeared*
                           (np.cosh(df_GEN.geneta_1-df_GEN.geneta_2)-np.cos(df_GEN.genphi_1-df_GEN.genphi_2)) )

                fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })

                hist_reco   = ax0.hist(df_RECO["MASS_Z_COR"], bins=bins, range=rang, density=True, color="k", alpha=0.5, label="reco")    
                gen_hist_old= ax0.hist(df_GEN.mass_Z_smeared, bins=bins, range=rang, density=True, histtype="step", color="r")
                gen_hist_new= ax0.hist(MZ,                    bins=bins, range=rang, density=True, histtype="step", color="b")
                
                x2_old=np.sum((gen_hist_old[0]-hist_reco[0])**2/hist_reco[0])
                x2_new=np.sum((gen_hist_new[0]-hist_reco[0])**2/hist_reco[0])

                ax0.text(85, 0.001, s=f"χ²old:{round(x2_old,4)}, χ²new:{round(x2_new,4)}")  
                ax0.set_ylabel("normalized event count")
                ax0.set_xlim(left=rang[0],right=rang[1])
                ax0.set_xticks([])
                ax0.legend()

                ax1.set_xlim(left=rang[0],right=rang[1])
                ax1.plot(np.linspace(rang[0],rang[1],bins),hist_reco[0]/gen_hist_old[0],".", c="r", label="RECO/old")
                ax1.plot(np.linspace(rang[0],rang[1],bins),hist_reco[0]/gen_hist_new[0],"_", c="b", label="RECO/new")
                ax1.set_xlabel("M_µµ (GeV)")
                ax1.set_ylabel("ratio")
                ax1.grid(True)
                ax1.legend()
                    
                plt.savefig(f"{pdir}{subtyp}_residual.png")
                plt.clf()


    json_object=json.dumps(k_results, indent=4)
    with open(f"{hdir}/residual_results.json", "w") as outfile:
        outfile.write(json_object)
