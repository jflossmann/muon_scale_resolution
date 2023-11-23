import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ROOT

def  iterative_correction(samples, eta_bins, phi_bins,  hdir, pdir):
    n_iterationsteps=40
    #write gen in dataframe
    
    file=uproot.open(samples["GEN"]["GEN"])
    tree=file["Events"]
    variables=tree.keys()
    df_GEN=tree.arrays(variables, library="pd")

    #cut events with 86<M_z_smeared<96
    mask_M_Z=(df_GEN.mass_Z_smeared>86) & (df_GEN.mass_Z_smeared<96)
    df_GEN_cut=df_GEN[mask_M_Z]
    
    #loop over samples
    
    for typ in samples:
        if typ=="DATA":#not typ=="GEN":
            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")
                #read data into dataframe
                file=uproot.open(samples[typ][subtyp])
                tree=file["Events"]
                variables=tree.keys()
                df_RECO=tree.arrays(variables, library="pd")

                ####################!!!!!!!!!!!!!!!!####################
                df_RECO["kappa_1"]=1
                df_RECO["lambda_1"]=0
                df_RECO["kappa_2"]=1
                df_RECO["lambda_2"]=0

                df_RECO["MASS_Z_COR"]=df_RECO["mass_Z_mean_roccor"]
                df_RECO["PT_1_COR"]=df_RECO["pt_1_mean_roccor"]
                df_RECO["PT_2_COR"]=df_RECO["pt_2_mean_roccor"]

                ########################################

                #table for kappa and lambda
                #... import from first step
                kappa_table=np.ones([len(eta_bins)-1,len(phi_bins)-1])
                lambd_table=np.zeros([len(eta_bins)-1,len(phi_bins)-1])

                #append empty df for current type to dfs_corr
                


                #iterate over eta, phi bins
                mass_means=[]
                iterations=[]
                for n in range(n_iterationsteps):
                    print("iteration:",n+1)

                    M_Z_filter=(86<df_RECO["MASS_Z_COR"]) & (df_RECO["MASS_Z_COR"]<96)

                    for i in tqdm(range(len(eta_bins)-1)):
                        #filter for eta bin
                        eta_1_filter=(eta_bins[i]<df_RECO["eta_1"]) & (df_RECO["eta_1"]<=eta_bins[i+1])
                        eta_2_filter=(eta_bins[i]<df_RECO["eta_2"]) & (df_RECO["eta_2"]<=eta_bins[i+1])
                        eta_1_filter_GEN=(eta_bins[i]<df_GEN_cut["geneta_1"]) & (df_GEN_cut["geneta_1"]<=eta_bins[i+1])
                        eta_2_filter_GEN=(eta_bins[i]<df_GEN_cut["geneta_2"]) & (df_GEN_cut["geneta_2"]<=eta_bins[i+1])

                        for j in range(len(phi_bins)-1):

                            #filter for phi bin
                            phi_1_filter=(phi_bins[j]<df_RECO["phi_1"]) & (df_RECO["phi_1"]<=phi_bins[j+1])
                            phi_2_filter=(phi_bins[j]<df_RECO["phi_2"]) & (df_RECO["phi_2"]<=phi_bins[j+1])
                            phi_1_filter_GEN=(phi_bins[j]<df_GEN_cut["genphi_1"]) & (df_GEN_cut["genphi_1"]<=phi_bins[j+1])
                            phi_2_filter_GEN=(phi_bins[j]<df_GEN_cut["genphi_2"]) & (df_GEN_cut["genphi_2"]<=phi_bins[j+1])

                            #select eta bin, phi bin and mass region
                            df_epm1=df_RECO.loc[M_Z_filter & eta_1_filter & phi_1_filter]
                            df_epm2=df_RECO.loc[M_Z_filter & eta_2_filter & phi_2_filter]

                            #select eta bin, phi bin for gen
                            df_ep1_GEN=df_GEN_cut.loc[eta_1_filter_GEN & phi_1_filter_GEN]
                            df_ep2_GEN=df_GEN_cut.loc[eta_2_filter_GEN & phi_2_filter_GEN]

                            #calculate kappa and lambda from equation-system
                            kappa, lambd = iterationsfunktion(df_gen_bin_p =  df_ep2_GEN,
                                                              df_reco_bin_p = df_epm2,
                                                              df_gen_bin_m =  df_ep1_GEN,
                                                              df_reco_bin_m = df_epm1)
                            

                            #update kappa and lambda for given bin

                            kappa_table[i][j]=kappa_table[i][j]*kappa
                            lambd_table[i][j]=kappa*lambd_table[i][j]+lambd

                            #print(f"l_gen1: {len(df_ep1_GEN)} l_gen2: {len(df_ep2_GEN)} m_Z_gen: {np.mean(df_ep1_GEN.mass_Z_smeared)} l1: {len(df_epm1)} l2: {len(df_epm2)} kap: {kappa_table[i][j]} lam: {lambd_table[i][j]}")

                            df_RECO.loc[eta_1_filter & phi_1_filter, "kappa"]=kappa_table[i][j]
                            df_RECO.loc[eta_1_filter & phi_1_filter, "lambda"]=lambd_table[i][j]

                            df_RECO.loc[eta_2_filter & phi_2_filter, "kappa"]=kappa_table[i][j]
                            df_RECO.loc[eta_2_filter & phi_2_filter, "lambda"]=lambd_table[i][j]
                            


                    
                    #calculate new corrected reconstructed pT for all events in df_RECO
                    df_RECO["PT_1_COR"]=1/(df_RECO["kappa"]/df_RECO["pt_1_mean_roccor"] - df_RECO["lambda"])
                    df_RECO["PT_2_COR"]=1/(df_RECO["kappa"]/df_RECO["pt_2_mean_roccor"] + df_RECO["lambda"])
                    #calculate new corrected reconstructed M_Z for all events in df_RECO
                    
                    df_RECO["MASS_Z_COR"]=np.sqrt( 2*df_RECO.PT_1_COR*df_RECO.PT_2_COR*(np.cosh(df_RECO.eta_1-df_RECO.eta_2)-np.cos(df_RECO.phi_1-df_RECO.phi_2)))
                    
                    #make lists of mass means for plot
                    mass_means.append(np.mean(df_RECO.loc[M_Z_filter,"MASS_Z_COR"]))
                    iterations.append(n+1)
                    
                    #define bins for plot
                    bins=300
                    rang=[60,120]
                        
                    if n==0:
                        #plt.hist(df_GEN["mass_Z_smeared"],bins=bins,range=rang,histtype="step",label="smeared gen",color='k',density=True)
                        plt.hist(df_RECO["mass_Z_mean_roccor"],bins=bins,range=rang,histtype="step",label="wo. it. corr",color='r',density=True)
                        plt.title(subtyp)
                        plt.xlabel("M_µµ (GeV)")
                            
                    k=int(n_iterationsteps/4)
                    if k<1:
                        k=1
                    if (n+1)%k==0:
                        plt.hist(df_RECO["MASS_Z_COR"],bins=bins,range=rang,histtype="step",label=f"iteration {n+1}",density=True)
                    if n==n_iterationsteps-1:
                        plt.legend()
                        plt.savefig(f"{pdir}iterative/{subtyp}_{n}.png")
                        plt.savefig(f"{pdir}iterative/{subtyp}_{n}.pdf")

            plt.clf()
            plt.plot(iterations,mass_means,"x",label="Corection")
            plt.plot(iterations, np.ones(len(iterations))*np.mean(df_GEN_cut.genmass_Z),label="Gen")
            plt.xlabel("iteration")
            plt.ylabel("mean(M_µµ)")
            plt.legend()
            plt.savefig(f"{pdir}iterative/iteration_mean_{n}.png")
            plt.savefig(f"{pdir}iterative/iteration_mean_{n}.pdf")

                #save data
                #print(f"saving corrected dataframe to /home/dguthmann/Documents/MA/Rochester/Daten/{subtyp}_iterative.ROOT")
                #data={key: df_RECO[key].values for key in df_RECO.columns}
                #rdf = ROOT.RDF.MakeNumpyDataFrame(data)
                #rdf.Snapshot("Events", f"/home/dguthmann/Documents/MA/Rochester/Daten/{subtyp}_iterative.ROOT")
                #print("done")                    


def iterationsfunktion(df_gen_bin_p, df_reco_bin_p, df_gen_bin_m, df_reco_bin_m):
    #define parameters of equation system
    Vp=-2*(np.mean(df_gen_bin_p.mass_Z_smeared)-np.mean(df_reco_bin_p["MASS_Z_COR"]))
    Vm=-2*(np.mean(df_gen_bin_m.mass_Z_smeared)-np.mean(df_reco_bin_m["MASS_Z_COR"]))

    Mp=2*np.mean(df_reco_bin_p["MASS_Z_COR"])
    Mm=2*np.mean(df_reco_bin_m["MASS_Z_COR"])

    Kp=np.mean(df_reco_bin_p["PT_2_COR"]*df_reco_bin_p["MASS_Z_COR"])
    Km=-np.mean(df_reco_bin_m["PT_1_COR"]*df_reco_bin_m["MASS_Z_COR"])

    #calc correction
    
    mu = (Vp/Kp-Vm/Km)/(Mp/Kp-Mm/Km)
    kappa=1+mu
    lambd = Vp/Kp-Mp/Kp*mu
    return kappa, lambd
