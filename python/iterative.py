import uproot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ROOT

def iterative_correction(samples, eta_bins, phi_bins, hdir, pdir):
    iterationsteps=20

    #get gen data
    file=uproot.open(samples["GEN"]["GEN"])
    print(samples["GEN"]["GEN"])
    tree=file["Events"]
    variables=tree.keys()
    df_GEN=tree.arrays(variables, library="pd")

    #make columns for corrected reconstructed gen-values
    df_GEN["kappa_1"]=0
    df_GEN["lambda_1"]=0
    df_GEN["kappa_2"]=0
    df_GEN["lambda_2"]=0

    df_GEN["MASS_Z_COR"]=df_GEN["mass_Z_mean_roccor"]
    df_GEN["PT_1_COR"]=df_GEN["pt_1_mean_roccor"]
    df_GEN["PT_2_COR"]=df_GEN["pt_2_mean_roccor"]

    #cut events with 86<M_z_smeared<96
    mask_M_Z=(df_GEN.genmass_Z>86) & (df_GEN.genmass_Z<96)
    df_GEN_m=df_GEN[mask_M_Z]

    for typ in samples:
        if typ == "SIG" or typ=="DATA":
            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")

                #read data into dataframe
                file=uproot.open(samples[typ][subtyp])
                tree=file["Events"]
                variables=tree.keys()
                df_RECO=tree.arrays(variables, library="pd")

                df_RECO["kappa_1"]=0
                df_RECO["lambda_1"]=0
                df_RECO["kappa_2"]=0
                df_RECO["lambda_2"]=0

                df_RECO["MASS_Z_COR"]=df_RECO["mass_Z_mean_roccor"]
                df_RECO["PT_1_COR"]=df_RECO["pt_1_mean_roccor"]
                df_RECO["PT_2_COR"]=df_RECO["pt_2_mean_roccor"]

                #table for kappa and lambda
                #... import from first step
                kappa_table, lambda_table=get_k_l(df_RECO=df_RECO, eta_bins=eta_bins, phi_bins=phi_bins)

                iterations=[]
                gen_means=[]
                reco_means=[]
                
                VM=[]
                VP=[]
                for n in tqdm(range(iterationsteps)):
                    iterations.append(n+1)
                    print("iteration ",n+1)
                    print("calc correction")

                    #mass cut on working region
                    M_Z_filter_RECO=(86<df_RECO["MASS_Z_COR"]) & (df_RECO["MASS_Z_COR"]<96)
                    df_RECO_m=df_RECO[M_Z_filter_RECO]

                    #loops to get correction
                    VM_n=[]
                    VP_n=[]

                    reco_means.append([])
                    gen_means.append([])
                    for i in range(len(eta_bins)-1):
                        #filter for eta bin
                        eta_1_filter_RECO=(eta_bins[i]<df_RECO_m["eta_1"]) & (df_RECO_m["eta_1"]<=eta_bins[i+1])
                        eta_2_filter_RECO=(eta_bins[i]<df_RECO_m["eta_2"]) & (df_RECO_m["eta_2"]<=eta_bins[i+1])
                        df_RECO_me1=df_RECO_m[eta_1_filter_RECO]
                        df_RECO_me2=df_RECO_m[eta_2_filter_RECO]

                        eta_1_filter_GEN=(eta_bins[i]<df_GEN_m["geneta_1"]) & (df_GEN_m["geneta_1"]<=eta_bins[i+1])
                        eta_2_filter_GEN=(eta_bins[i]<df_GEN_m["geneta_2"]) & (df_GEN_m["geneta_2"]<=eta_bins[i+1]) 
                        df_GEN_me1=df_GEN_m[eta_1_filter_GEN]
                        df_GEN_me2=df_GEN_m[eta_2_filter_GEN]

                        gen_means[n].append([])
                        reco_means[n].append([])
                        for j in range(len(phi_bins)-1):
                            phi_1_filter_RECO=(phi_bins[j]<df_RECO_me1["phi_1"]) & (df_RECO_me1["phi_1"]<=phi_bins[j+1])
                            phi_2_filter_RECO=(phi_bins[j]<df_RECO_me2["phi_2"]) & (df_RECO_me2["phi_2"]<=phi_bins[j+1])
                            df_RECO_mep1=df_RECO_me1[phi_1_filter_RECO]
                            df_RECO_mep2=df_RECO_me2[phi_2_filter_RECO]

                            phi_1_filter_GEN=(phi_bins[j]<df_GEN_me1["genphi_1"]) & (df_GEN_me1["genphi_1"]<=phi_bins[j+1])
                            phi_2_filter_GEN=(phi_bins[j]<df_GEN_me2["genphi_2"]) & (df_GEN_me2["genphi_2"]<=phi_bins[j+1])
                            df_GEN_mep1=df_GEN_me1[phi_1_filter_GEN]
                            df_GEN_mep2=df_GEN_me2[phi_2_filter_GEN]
  
                            gen_means[n][i].append(np.mean(df_GEN_mep1["genmass_Z"]))
                            reco_means[n][i].append(np.mean(df_RECO_mep1["MASS_Z_COR"]))

                            #calculate correction
                            Vp=-2*(np.mean(df_GEN_mep2["genmass_Z"])-np.mean(df_RECO_mep2["MASS_Z_COR"]))
                            VP_n.append(Vp)
                            Vm=-2*(np.mean(df_GEN_mep1["genmass_Z"])-np.mean(df_RECO_mep1["MASS_Z_COR"]))
                            VM_n.append(Vm)
                            Mp=2*np.mean(df_RECO_mep2["MASS_Z_COR"])
                            Mm=2*np.mean(df_RECO_mep1["MASS_Z_COR"])
                            Kp=np.mean(df_RECO_mep2["PT_2_COR"]*df_RECO_mep2["MASS_Z_COR"])
                            Km=-np.mean(df_RECO_mep1["PT_1_COR"]*df_RECO_mep1["MASS_Z_COR"])
                            kappa, lambd = iterationsfunktion(Mp=Mp, Mm=Mm, Vp=Vp, Vm=Vm, Kp=Kp, Km=Km)
                            
                            #update kappa and lambda table for given bin
                            kappa_table[i][j]=kappa_table[i][j]*kappa
                            lambda_table[i][j]=kappa*lambda_table[i][j]+lambd

                    VP.append(np.mean(VP_n))
                    VM.append(np.mean(VM_n))

                    #apply corrections
                    print("apply correction")

                    #correct gen-reco values
                    if subtyp=="SIG":
                        for i in range(len(eta_bins)-1):
                            #filter for eta bin
                            eta_1_filter=(eta_bins[i]<df_GEN["eta_1"]) & (df_GEN["eta_1"]<=eta_bins[i+1])
                            eta_2_filter=(eta_bins[i]<df_GEN["eta_2"]) & (df_GEN["eta_2"]<=eta_bins[i+1])
                            
                            for j in range(len(phi_bins)-1):
                                phi_1_filter=(phi_bins[j]<df_GEN["phi_1"]) & (df_GEN["phi_1"]<=phi_bins[j+1])
                                phi_2_filter=(phi_bins[j]<df_GEN["phi_2"]) & (df_GEN["phi_2"]<=phi_bins[j+1])
                                
                                df_GEN.loc[eta_1_filter & phi_1_filter, "kappa_1"]=kappa_table[i][j]
                                df_GEN.loc[eta_1_filter & phi_1_filter, "lambda_1"]=lambda_table[i][j]

                                df_GEN.loc[eta_2_filter & phi_2_filter, "kappa_2"]=kappa_table[i][j]
                                df_GEN.loc[eta_2_filter & phi_2_filter, "lambda_2"]=lambda_table[i][j]

                        #calc correctet muon pt and z_mass
                        df_GEN["PT_1_COR"]=1/(df_GEN["kappa_1"]/df_GEN["pt_1"] - df_GEN["lambda_1"])
                        df_GEN["PT_2_COR"]=1/(df_GEN["kappa_2"]/df_GEN["pt_2"] + df_GEN["lambda_2"])
                        df_GEN["MASS_Z_COR"]=np.sqrt( 2*df_GEN.PT_1_COR*df_GEN.PT_2_COR*(np.cosh(df_GEN.eta_1-df_GEN.eta_2)-np.cos(df_GEN.phi_1-df_GEN.phi_2)))

                        #save corrected gen results
                        print(f"saving corrected gen to {samples['GEN']['GEN']}")
                        data={key: df_GEN[key].values for key in df_GEN.columns}
                        rdf = ROOT.RDF.MakeNumpyDataFrame(data)
                        rdf.Snapshot("Events", samples["GEN"]["GEN"])
                        print("done")

                    #correct reco values
                    for i in range(len(eta_bins)-1):
                        #filter for eta bin
                        eta_1_filter=(eta_bins[i]<df_RECO["eta_1"]) & (df_RECO["eta_1"]<=eta_bins[i+1])
                        eta_2_filter=(eta_bins[i]<df_RECO["eta_2"]) & (df_RECO["eta_2"]<=eta_bins[i+1])
                        
                        for j in range(len(phi_bins)-1):
                            phi_1_filter=(phi_bins[j]<df_RECO["phi_1"]) & (df_RECO["phi_1"]<=phi_bins[j+1])
                            phi_2_filter=(phi_bins[j]<df_RECO["phi_2"]) & (df_RECO["phi_2"]<=phi_bins[j+1])
                            
                            df_RECO.loc[eta_1_filter & phi_1_filter, "kappa_1"]=kappa_table[i][j]
                            df_RECO.loc[eta_1_filter & phi_1_filter, "lambda_1"]=lambda_table[i][j]

                            df_RECO.loc[eta_2_filter & phi_2_filter, "kappa_2"]=kappa_table[i][j]
                            df_RECO.loc[eta_2_filter & phi_2_filter, "lambda_2"]=lambda_table[i][j]

                    #calc correctet muon pt and z_mass
                    df_RECO["PT_1_COR"]=1/(df_RECO["kappa_1"]/df_RECO["pt_1"] - df_RECO["lambda_1"])
                    df_RECO["PT_2_COR"]=1/(df_RECO["kappa_2"]/df_RECO["pt_2"] + df_RECO["lambda_2"])
                    df_RECO["MASS_Z_COR"]=np.sqrt( 2*df_RECO.PT_1_COR*df_RECO.PT_2_COR*(np.cosh(df_RECO.eta_1-df_RECO.eta_2)-np.cos(df_RECO.phi_1-df_RECO.phi_2)))

            
            
                   
            #define bins for plot
            bins=80
            rang=[86,96]        
            fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
            plt.subplots_adjust(wspace=0, hspace=0)
            ax0.set_xlim(left=rang[0],right=rang[1])
            ax0.set_xticks([])
            h_gen=ax0.hist(df_GEN["mass_Z_smeared"],bins=bins,range=rang,histtype="step",label="smeared gen",color='r',density=True)
            h_reco=ax0.hist(df_RECO["MASS_Z_COR"],bins=bins,range=rang,histtype="step",label=f"iteration {n+1}",density=True)
            h_reco0=ax0.hist(df_RECO["mass_Z_mean_roccor"],bins=bins,range=rang,histtype="step",label=f"iteration {0}",density=True)
            ax0.legend()

            ax1.set_ylim(bottom=0.85,top=1.15)
            ax1.set_xlim(left=rang[0],right=rang[1])
            ax1.plot(np.linspace(rang[0],rang[1],bins),h_reco[0]/h_gen[0],'.',label="reco_15it/gen")
            ax1.plot(np.linspace(rang[0],rang[1],bins),h_reco0[0]/h_gen[0],'.',label="reco_0it/gen")
            ax1.set_xlabel("M_µµ (GeV)")
            ax1.set_ylabel("ratio")
            ax1.legend()
            ax1.grid(True)
            plt.savefig(f"{pdir}iterative/{subtyp}_z_mass_cut_{n+1}.png")
            plt.clf()

            bins=100
            rang=[60,120]        
            fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
            plt.subplots_adjust(wspace=0, hspace=0)
            ax0.set_xlim(left=rang[0],right=rang[1])
            ax0.set_xticks([])
            h_gen=ax0.hist(df_GEN["mass_Z_smeared"],bins=bins,range=rang,histtype="step",label="smeared gen",color='r',density=True)
            h_reco=ax0.hist(df_RECO["MASS_Z_COR"],bins=bins,range=rang,histtype="step",label=f"iteration {n+1}",density=True)
            h_reco0=ax0.hist(df_RECO["mass_Z_mean_roccor"],bins=bins,range=rang,histtype="step",label=f"iteration {0}",density=True)
            ax0.legend()

            ax1.set_ylim(bottom=0.85,top=1.15)
            ax1.set_xlim(left=rang[0],right=rang[1])
            ax1.plot(np.linspace(rang[0],rang[1],bins),h_reco[0]/h_gen[0],'.',label="reco_15it/gen")
            ax1.plot(np.linspace(rang[0],rang[1],bins),h_reco0[0]/h_gen[0],'.',label="reco_0it/gen")
            ax1.set_xlabel("M_µµ (GeV)")
            ax1.set_ylabel("ratio")
            ax1.legend()
            ax1.grid(True)
            plt.savefig(f"{pdir}iterative/{subtyp}_z_mass_{n+1}.png")
            plt.clf()
                           

            #make binwise plot
            for i in range(len(eta_bins)-1):
                for j in range(len(phi_bins)-1):
                    reco=[]
                    gen=[]
                    for n in range(len(reco_means)):
                        reco.append(reco_means[n][i][j])
                        gen.append(gen_means[n][i][j])

                    plt.clf()
                    plt.plot(iterations, reco,"x",c="blue",label="RECO")
                    plt.plot(iterations, gen,"--",c="darkblue",label="GEN")
                    plt.xlabel("iteration")
                    plt.ylabel("mean(M_µµ)")
                    plt.title(f"eta[{eta_bins[i]}, {eta_bins[i+1]}], phi[{round(phi_bins[i],1)}, {round(phi_bins[i+1],1)}]")
                    plt.legend()
                    plt.savefig(f"{pdir}iterative/binwise/{subtyp}/{subtyp}iteration_mean_eta{i}_phi{j}.png")
                    plt.clf()
            
            plt.plot(iterations, VP)
            plt.plot(iterations,VM)
            plt.savefig(f"{pdir}iterative/{subtyp}iteration_min_var.png")
            plt.clf()

            #save results
            print(f"saving corrected dfs to {samples[typ][subtyp]}")
            data={key: df_RECO[key].values for key in df_RECO.columns}
            rdf = ROOT.RDF.MakeNumpyDataFrame(data)
            rdf.Snapshot("Events", samples[typ][subtyp])
            print("done")  

def iterationsfunktion(Mp,Mm,Vp,Vm,Kp,Km):

    #calc correction
    
    mu = (Vp/Kp-Vm/Km)/(Mp/Kp-Mm/Km)
    kappa=1+mu
    lambd = ((Vp-Mp*mu)/Kp+(Vm-Mm*mu)/Km)/2
    return kappa, lambd


def get_k_l(df_RECO, eta_bins, phi_bins):
    print("getting initial kappa and lambda")
    kappa_table=np.zeros([len(eta_bins)-1,len(phi_bins)-1])
    lambda_table=np.zeros([len(eta_bins)-1,len(phi_bins)-1])

    for i in range(len(eta_bins)-1):
        #filter for eta bin
        eta_1_filter=(eta_bins[i]<df_RECO["eta_1"]) & (df_RECO["eta_1"]<=eta_bins[i+1])
        eta_2_filter=(eta_bins[i]<df_RECO["eta_2"]) & (df_RECO["eta_2"]<=eta_bins[i+1])
        df_RECO_e1=df_RECO[eta_1_filter]
        df_RECO_e2=df_RECO[eta_2_filter]

        for j in range(len(phi_bins)-1):
            #filter for phi bin
            phi_1_filter=(phi_bins[j]<df_RECO_e1["phi_1"]) & (df_RECO_e1["phi_1"]<=phi_bins[j+1])
            phi_2_filter=(phi_bins[j]<df_RECO_e2["phi_2"]) & (df_RECO_e2["phi_2"]<=phi_bins[j+1])
            df_RECO_ep1=df_RECO_e1[phi_1_filter]
            df_RECO_ep2=df_RECO_e2[phi_2_filter]

            pt1=np.array(df_RECO_ep1.pt_1)[0]
            pt1_cor=np.array(df_RECO_ep1.pt_1_mean_roccor)[0]
            pt2=np.array(df_RECO_ep2.pt_2)[0]
            pt2_cor=np.array(df_RECO_ep2.pt_2_mean_roccor)[0]
            
        
            kappa_table[i][j]+=pt1*pt2*(pt1_cor+pt2_cor)/(pt1_cor*pt2_cor*(pt1+pt2))
            lambda_table[i][j]+=1/pt2_cor-kappa_table[i][j]/pt2

            #print(kappa_table[i][j],lambda_table[i][j], pt1, pt1_cor, 1/(kappa_table[i][j]/pt1-lambda_table[i][j]))

    return kappa_table, lambda_table
