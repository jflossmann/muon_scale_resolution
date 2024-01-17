import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import uproot
from tqdm import tqdm
import json
import os
from scipy.stats import crystalball


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
    return {"val": X2[m], "x": K[m]}

def crystal_ball(x, alpha, n, mu, sigma):  # this is the parametrization from wikipedia
    mu = np.zeros(len(mu))
    x = x - np.ones(len(x)) * np.mean(x)
    return  crystalball.pdf(x=-x, beta=alpha, m=n, loc=-mu, scale=sigma)  # scipy uses different parameter names

def chi2(k, df_gen, reco_hist, bins=50, rang=[81,101]):

    #calc weight
    x1, x2         = df_gen.x1, df_gen.x2
    alpha1, alpha2 = df_gen.CB_alpha_1, df_gen.CB_alpha_2
    n1, n2         = df_gen.CB_n_1, df_gen.CB_n_2
    mu1, mu2       = df_gen.CB_mean_1, df_gen.CB_mean_2
    sig1, sig2     = df_gen.CB_sigma_1*df_gen.std1, df_gen.CB_sigma_2*df_gen.std2
    
    w1=crystal_ball(x1, alpha1, n1, mu1, sig1*k)/crystal_ball(x1, alpha1, n1, mu1, sig1)
    w2=crystal_ball(x2, alpha2, n2, mu2, sig2*k)/crystal_ball(x2, alpha2, n2, mu2, sig2)
    
    weights=w1*w2
    weights=np.nan_to_num(weights)
    weights = weights/np.mean(weights)
    
    #make histogram
    gen_hist=np.histogram(df_gen.mass_Z_smeared, bins=bins, range=rang, weights=weights)

    #calc chi2/me/mse

    #x2
    x2=np.sum((gen_hist[0]-reco_hist[0])**2/reco_hist[0])

    #me
    #me=np.sum(np.abs(gen_hist[0]-reco_hist[0]))/bins

    #mse
    #me=np.sum((gen_hist[0]-reco_hist[0])**2)/bins
    
    return x2

def residual_correction(samples, abseta_bins, hdir, pdir):
    #dict for results
    k_results={}
    bins=60
    rang=[80,102]
    n=8
    #get gen data
    file=uproot.open(samples["GEN"]["GEN"])
    tree=file["Events"]
    variables=tree.keys()
    df_GEN=tree.arrays(variables, library="pd")

    #cut events with 81<M_z_smeared<101
    mask_M_Z=(df_GEN.genmass_Z>rang[0]) & (df_GEN.genmass_Z<rang[1])
    df_GEN_m=df_GEN[mask_M_Z]
    df_GEN_m["genabseta_1"]=np.abs(df_GEN_m.geneta_1)
    df_GEN_m["genabseta_2"]=np.abs(df_GEN_m.geneta_2)

    #calc "x"
    df_GEN_m["x1"]=df_GEN_m.genpt_1_smeared/df_GEN_m["PT_1_COR"] - 1
    df_GEN_m["x2"]=df_GEN_m.genpt_2_smeared/df_GEN_m["PT_2_COR"] - 1

    plt.hist(df_GEN_m["x1"],range=[-0.2,0.2], bins=50,histtype="step", label="x1")
    plt.hist(df_GEN_m["x2"],range=[-0.2,0.2], bins=50,histtype="step", label="x2")
    plt.legend()
    os.makedirs(f"{pdir}residual/fits", exist_ok=True)
    plt.savefig(f"{pdir}residual/fits/x_hist_all.png")
    plt.clf()


    
    #filter events with very poorly reconstructed muons
    bad_1_filter = (np.abs(df_GEN_m["x1"]) < 2*df_GEN_m["CB_sigma_1"]*df_GEN_m["std1"])
    bad_2_filter = (np.abs(df_GEN_m["x2"]) < 2*df_GEN_m["CB_sigma_2"]*df_GEN_m["std2"])
    
    df_GEN_m=df_GEN_m[bad_1_filter]
    df_GEN_m=df_GEN_m[bad_2_filter]
    

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

                    abseta_1_filter_reco=(abseta_bins[i]<df_RECO["abseta_1"]) & (df_RECO["abseta_1"]<=abseta_bins[i+1])
                    abseta_2_filter_reco=(abseta_bins[i]<df_RECO["abseta_2"]) & (df_RECO["abseta_2"]<=abseta_bins[i+1])

                    df_RECO_me=df_RECO_m[abseta_1_filter_reco & abseta_2_filter_reco]

                    abseta_1_filter_gen=(abseta_bins[i]<df_GEN_m["genabseta_1"]) & (df_GEN_m["genabseta_1"]<=abseta_bins[i+1])
                    abseta_2_filter_gen=(abseta_bins[i]<df_GEN_m["genabseta_2"]) & (df_GEN_m["genabseta_2"]<=abseta_bins[i+1])

                    df_GEN_me=df_GEN_m[abseta_1_filter_gen & abseta_2_filter_gen]

                    plt.hist(df_GEN_me["x1"],range=[-0.05,0.05], bins=50,histtype="step", label="x1")
                    plt.hist(df_GEN_me["x2"],range=[-0.05,0.05], bins=50,histtype="step", label="x2")
                    plt.legend()
                    plt.savefig(f"{pdir}residual/fits/x_hist_{i}.png")
                    plt.clf()

                    n_gen=len(df_GEN_me)
                    n_reco=len(df_RECO_me)
                    scale_weights_reco=n_gen/n_reco*np.ones(n_reco)
                    hist_reco=list(plt.hist(df_RECO_me["MASS_Z_COR"],bins=bins, weights=scale_weights_reco,  range=rang, color="k", alpha=0.5, label="reco"))
                    
                    #fit
                    k=minfinder(fun=chi2, args=(df_GEN_me, hist_reco, bins, rang), bounds=[0.8,10.0], n=10, N=n)

                    K.append(k["x"])

                    x1, x2         = df_GEN_me.x1, df_GEN_me.x2
                    alpha1, alpha2 = df_GEN_me.CB_alpha_1, df_GEN_me.CB_alpha_2
                    n1, n2         = df_GEN_me.CB_n_1, df_GEN_me.CB_n_2
                    mu1, mu2       = df_GEN_me.CB_mean_1, df_GEN_me.CB_mean_2
                    sig1, sig2     = df_GEN_me.CB_sigma_1*df_GEN_me.std1, df_GEN_me.CB_sigma_2*df_GEN_me.std2

                    print("Error check: x mean", np.mean(x1), np.mean(x2))
                    print("Error check: x stds", np.std(x1), np.std(x2))
                    print("Error check: x stds", np.std(x1)/np.sqrt(len(x1)), np.std(x2)/np.sqrt(len(x2)))

                    w1=crystal_ball(x1, alpha1, n1, mu1, sig1*K[i])/crystal_ball(x1, alpha1, n1, mu1, sig1)
                    w2=crystal_ball(x2, alpha2, n2, mu2, sig2*K[i])/crystal_ball(x2, alpha2, n2, mu2, sig2)
                    
                    weights=w1*w2
                    weights=np.nan_to_num(weights)
                    weights = weights/np.mean(weights)

                    #errorchecking
                    if i == 0:
                        print(K[i])
                        Ks = [0.8,1,1.2,1.4, 1.6, 1.8, 5.0, 20.]

                        #plot fit for different ks
                        for K_test in Ks:
                            w1=crystal_ball(x1, alpha1, n1, mu1, sig1*K_test)/crystal_ball(x1, alpha1, n1, mu1, sig1)
                            w2=crystal_ball(x2, alpha2, n2, mu2, sig2*K_test)/crystal_ball(x2, alpha2, n2, mu2, sig2)
                        
                            weights=w1*w2
                            weights=np.nan_to_num(weights)
                            weights = weights/np.mean(weights)
                            resid=plt.hist(df_GEN_me["mass_Z_smeared"], bins=bins, range=rang,  weights=weights, histtype="step",label="k="+str(K_test))
                            _chi2 = np.sum( (resid[0] - hist_reco[0])**2 / hist_reco[0] )
                            print("alternative chi squareds", K_test, _chi2)
                        plt.legend()
                        plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_test_fit.png")
                        plt.clf()


                        plt.hist2d(x1,df_GEN_me.mass_Z_smeared,bins=40,range=[[-.05,.05],rang],label="k="+str(K_test))
                        plt.colorbar()
                        plt.ylabel("mass")
                        plt.xlabel("x1")

                        plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_x1mass_2d.png")
                        plt.clf()

                        plt.hist2d(x2,df_GEN_me.mass_Z_smeared,bins=40,range=[[-.05,.05],rang],label="k="+str(K_test))
                        plt.colorbar()
                        plt.ylabel("mass")
                        plt.xlabel("x2")

                        plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_x2mass_2d.png")
                        plt.clf()

                        plt.hist2d(x1,x2,bins=40,range=[[-.05,.05],[-.05,.05]],label="k="+str(K_test))
                        plt.colorbar()
                        plt.xlabel("x1")
                        plt.ylabel("x2")

                        plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_x1x2_2d.png")
                        plt.clf()

                        for K_test in Ks:
                            w1=crystal_ball(x1, alpha1, n1, mu1, sig1*K_test)/crystal_ball(x1, alpha1, n1, mu1, sig1)
                            w2=crystal_ball(x2, alpha2, n2, mu2, sig2*K_test)/crystal_ball(x2, alpha2, n2, mu2, sig2)
                        
                            weights=w1*w2
                            weights=np.nan_to_num(weights)
                            weights = weights/np.mean(weights)
                            m=np.mean(weights)

                            plt.hist2d(w1,x1,bins=50,range=[[0.5,1.5],[-0.05,0.05]],label="k="+str(K_test))
                            plt.colorbar()
                            plt.ylabel("x1")
                            plt.xlabel("weight1")

                            plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_{K_test}_weights1_2d.png")
                            plt.clf()

                            plt.hist2d(w2,x2,bins=50,range=[[0.5,1.5],[-0.05,0.05]],label="k="+str(K_test))
                            plt.colorbar()
                            plt.ylabel("x2")
                            plt.xlabel("weight2")

                            plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_{K_test}_weights2_2d.png")
                            plt.clf()

                            plt.hist2d(weights,df_GEN_me.mass_Z_smeared,bins=40,range=[[0.2,1.5],rang],label="k="+str(K_test))
                            plt.colorbar()
                            plt.ylabel("mass")
                            plt.xlabel("weight")

                            plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_{K_test}_weights_2d.png")
                            plt.clf()
                            print("Weight details:", K_test)
                            print("Mean, std", np.mean(weights), np.std(weights))
                            print("Min, max", min(weights), max(weights))

                        #for K_test in [0.5,0.75,0.95,1.05,1.25,1.5]:
                        #    w1=crystal_ball(x1, alpha1, n1, mu1, sig1*K_test)/crystal_ball(x1, alpha1, n1, mu1, sig1)
                        #    w2=crystal_ball(x2, alpha2, n2, mu2, sig2*K_test)/crystal_ball(x2, alpha2, n2, mu2, sig2)
                       # 
                       #     weights=(w1+w2)/2
                       #     weights=np.nan_to_num(weights)
                       #     plt.hist2d(weights,x1,bins=50,range=[0.5,1.5],histtype="step",label="k="+str(K_test)+" max:"+str(round(max(weights))))
                       #     plt.legend()
                       #     plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_weights.png")
                       #     plt.clf()

                    i_max=np.argmax(weights)
                    print("pt1_gen_s: ", df_GEN_me.iloc[i_max]["genpt_1_smeared"], 
                          " pt1_reco: ", df_GEN_me.iloc[i_max]["PT_1_COR"], 
                          " geneta_1: ", df_GEN_me.iloc[i_max]["geneta_1"],
                          " x1: ", df_GEN_me.iloc[i_max]["x1"],
                          " w1: ", w1[i_max],
                          " sig:", df_GEN_me.iloc[i_max]["CB_sigma_1"]*df_GEN_me.iloc[i_max]["std1"])
                    print("pt2_gen_s: ", df_GEN_me.iloc[i_max]["genpt_2_smeared"], 
                          " pt2_reco: ", df_GEN_me.iloc[i_max]["PT_2_COR"], 
                          " geneta_2: ", df_GEN_me.iloc[i_max]["geneta_2"],
                          " x2: ", df_GEN_me.iloc[i_max]["x2"],
                          " w2: ", w2[i_max],
                          " sig:", df_GEN_me.iloc[i_max]["CB_sigma_2"]*df_GEN_me.iloc[i_max]["std2"])
                    
                
                        #plot CB_distr
                        #x=np.linspace(-0.1,0.1,500)
                        #for K_test in [0.95,1,1.01,1.02]:
                        #    plt.plot(x,crystal_ball(x, np.mean(alpha1), np.mean(n1), np.mean(mu1), np.mean(sig1)*K_test),label="k="+str(K_test))
                        #plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_CB.png")
                        #plt.clf()
                        #x=np.linspace(-0.2,0.1,500)
                        #for K_test in [0.95,1,1.01,1.02]:
                        #    plt.plot(x,crystal_ball(x, np.mean(alpha1), np.mean(n1), np.mean(mu1), np.mean(sig1)*K_test)/crystal_ball(x, np.mean(alpha1), np.mean(n1), np.mean(mu1), np.mean(sig1)),label="k="+str(K_test))
                        #plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_CB_ratio.png")
                        #plt.clf()

                        #plt.plot(x,crystal_ball(x, df_GEN_me.iloc[i_max]["CB_alpha_1"], df_GEN_me.iloc[i_max]["CB_n_1"], df_GEN_me.iloc[i_max]["CB_mean_1"], df_GEN_me.iloc[i_max]["CB_sigma_1"]*df_GEN_me.iloc[i_max]["std1"]*K_test))
                        #plt.plot(x,crystal_ball(x, df_GEN_me.iloc[i_max]["CB_alpha_1"], df_GEN_me.iloc[i_max]["CB_n_1"], df_GEN_me.iloc[i_max]["CB_mean_1"], df_GEN_me.iloc[i_max]["CB_sigma_1"]*df_GEN_me.iloc[i_max]["std1"]))
                        #plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_max_weight_CB.png")
                        #plt.clf()
                    
                    
                    w1=crystal_ball(x1, alpha1, n1, mu1, sig1*K[i])/crystal_ball(x1, alpha1, n1, mu1, sig1)
                    w2=crystal_ball(x2, alpha2, n2, mu2, sig2*K[i])/crystal_ball(x2, alpha2, n2, mu2, sig2)
                    
                    weights=w1*w2
                    weights=np.nan_to_num(weights)
                    weights = weights/np.mean(weights)
                    resid=plt.hist(df_GEN_me["mass_Z_smeared"], bins=bins, range=rang,  weights=weights, histtype="step", color="b",label="residual")
                    
                    
                    smear=plt.hist(df_GEN_me["mass_Z_smeared"],bins=bins, range=rang, histtype="step",color="r", label="smeared")
                    x2_smeared=np.sum((smear[0]-hist_reco[0])**2/hist_reco[0])
                    x2_residual=np.sum((resid[0]-hist_reco[0])**2/hist_reco[0])
                    plt.text(91, 0.001, s=f"s:{round(x2_smeared,5)}, r:{round(x2_residual,5)}")
                    plt.legend()
                    plt.savefig(f"{pdir}residual/fits/{subtyp}_eta_{i}_fit.png")
                    plt.clf()

                k_results[typ]=K
            
    json_object=json.dumps(k_results, indent=4)
    with open(f"{hdir}/residual_results.json", "w") as outfile:
        outfile.write(json_object)
