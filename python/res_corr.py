import uproot
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

import pandas as pd
pd.set_option('mode.chained_assignment', None)
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ROOT
import os
import json
from array import array


def get_res_correction(ntuples_gen, pull_bins, abseta_bins, nl_bins, pt_bins, pdir , hdir, do_plot):
    #read data
    pdir = pdir+'resolution/'
    os.makedirs(pdir, exist_ok=True)
    if do_plot:
        os.makedirs(pdir+'CB_fits/', exist_ok=True)
        os.makedirs(pdir+'pol_fits/', exist_ok=True)
    ROOT.gROOT.SetBatch(1)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    file=uproot.open(ntuples_gen)
    tree=file["Events"]
    variables=tree.keys()
    df=tree.arrays(variables, library="pd")

    #calculate abseta
    df["abseta_1"]=np.abs(df.eta_1)
    df["abseta_2"]=np.abs(df.eta_2)

    #calculate R
    df["R_1"]=df["genpt_1"]/df["pt_1"]
    df["R_2"]=df["genpt_2"]/df["pt_2"]


    #make histogram of distribution in abseta, nl
    if do_plot:
        plt.hist2d(df["abseta_1"],df["nTrkLayers_1"],bins=[abseta_bins, nl_bins])
        plt.colorbar()
        plt.xlabel('|eta|')
        plt.ylabel('Number of Layers')
        plt.savefig(f"./{pdir}/nL_abseta.png")
        plt.savefig(f"./{pdir}/nL_abseta.pdf")
        plt.clf()

        #make x axis for polynom plot later
        x_ax=np.linspace(min(pt_bins),max(pt_bins),200)

    
    fit_results={"pull":{}, "poly":{}}
    hist_std=[]
    hist_std_err=[]

    #iterate over eta bins
    for i in tqdm(range(len(abseta_bins)-1)):
        
        for f in fit_results:
            fit_results[f]["eta_"+str(i)]={}

        hist_std.append([])
        hist_std_err.append([])
        

        eta_1_filter=(abseta_bins[i]<df["eta_1"]) & (df["eta_1"]<=abseta_bins[i+1])
        eta_2_filter=(abseta_bins[i]<df["eta_2"]) & (df["eta_2"]<=abseta_bins[i+1])

        df_e1=df[eta_1_filter]
        df_e2=df[eta_2_filter]

        #iterate over nl bins
        for j in range(len(nl_bins)-1):

            for f in fit_results:
                fit_results[f]["eta_"+str(i)]["nL_"+str(j)]={}

            hist_std[i].append([])
            hist_std_err[i].append([])


            nl_1_filter=(nl_bins[j]<df_e1["nTrkLayers_1"]) & (df_e1["nTrkLayers_1"]<=nl_bins[j+1])
            nl_2_filter=(nl_bins[j]<df_e2["nTrkLayers_2"]) & (df_e2["nTrkLayers_2"]<=nl_bins[j+1])

            df_en1=df_e1[nl_1_filter]
            df_en2=df_e2[nl_2_filter]

            #calculate pull distribution
            R_1=df_en1["R_1"]
            R_2=df_en2["R_2"]
            R=pd.concat([R_1,R_2])

            std=np.std(R)
            pull=(R-np.mean(R))/std
            pull_hist=np.histogram(pull, bins=pull_bins)[0]

            # Create a ROOT TH1 object
            hist = ROOT.TH1D("hist", "PullHistogram", len(pull_bins) - 1, pull_bins)
            # Fill the histogram with data
            for k, value in enumerate(pull_hist):
                hist.SetBinContent(k+1, value)

            if hist.Integral() == 0:
                continue

            hist.Scale(1./hist.Integral())

            x = ROOT.RooRealVar("x", "m_vis (GeV)", -5, 5)
            x.setBins(10000,"cache")
            x.setMin("cache", -10)
            x.setMax("cache", 10)

            mean = ROOT.RooRealVar("mean", "mean", 0, -1, 1)
            sigma = ROOT.RooRealVar("sigma", "sigma", hist.GetStdDev(), 0, 5)
            n = ROOT.RooRealVar("n", "n", 10, 0, 1000)
            alpha = ROOT.RooRealVar("alpha", "alpha", 1, 0, 5)
            cb = ROOT.RooCrystalBall("cb", "CrystalBall", x, mean, sigma,
                                                sigma, alpha, n, alpha, n)
            roohist = ROOT.RooDataHist("roohist", "", ROOT.RooArgSet(x), hist)
            fitResult = cb.fitTo(roohist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintEvalErrors(-1))

            fit_parameters = [mean.getVal(), sigma.getVal()*2, n.getVal(), alpha.getVal()]
            parameter_errors = [mean.getError(), sigma.getError()*2, n.getError(), alpha.getError()]

            # save the fit results

            for k, (param, error) in enumerate(zip(fit_parameters, parameter_errors)):
                #print(f"Parameter {k}: {param} +/- {error}")
                fit_results["pull"]["eta_"+str(i)]["nL_"+str(j)][f"Parameter {k}"]={"value":param, "error":error}
                
            
            if do_plot:
                # Create a canvas for plotting
                c1 = ROOT.TCanvas("c1", "Fitted Histogram", 800, 600)
                frame = x.frame()
                roohist.plotOn(frame, ROOT.RooFit.DrawOption("B"), ROOT.RooFit.FillStyle(0), ROOT.RooFit.FillColor(ROOT.kBlue))
                cb.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue))
                frame.Draw()

                # Display the canvas
                c1.Update()
                c1.Modified()

                # Optionally, save the plot as an image
                c1.SaveAs(f"./{pdir}CB_fits/eta{i}_nL{j}.png")
                    
            
            # bin in pt_gen for polynomial fit
            for k in range(len(pt_bins)-1):
                pt_1_filter=(pt_bins[k]<df_en1["genpt_1"]) & (df_en1["genpt_1"]<=pt_bins[k+1])
                pt_2_filter=(pt_bins[k]<df_en2["genpt_2"]) & (df_en2["genpt_2"]<=pt_bins[k+1])

                df_enp1=df_en1[pt_1_filter]
                df_enp2=df_en2[pt_2_filter]

                R_1=df_enp1["R_1"]
                R_2=df_enp2["R_2"]
                R=pd.concat([R_1,R_2])
                if not len(R) == 0:
                    hist_std[i][j].append(np.std(R))
                    hist_std_err[i][j].append(np.std(R)/np.sqrt(2*len(R)))
           
            hist = ROOT.TH1D("hist", "PolyHistogram", len(pt_bins) - 1, array('d', pt_bins))
            for k, value in enumerate(hist_std[i][j]):
                hist.SetBinContent(k+1, value)
                hist.SetBinError(k+1, hist_std_err[i][j][k])

            # Define the second-order polynomial function
            polynomial = ROOT.TF1("polynomial", "[0] + [1]*x + [2]*x*x")

            # Set initial parameter values
            polynomial.SetParameter(0, 0.01)  # Coefficient for the constant term
            polynomial.SetParameter(1, 5e-5)   # Coefficient for the linear term
            polynomial.SetParameter(2, 5e-5)   # Coefficient for the quadratic term

            # Fit the histogram with the plynomial function
            hist.Fit("polynomial")
            # Get the fit results
            fit_parameters = [polynomial.GetParameter(i) for i in range(3)]
            parameter_errors = [polynomial.GetParError(i) for i in range(3)]

            # Print the fit results
            for k, (param, error) in enumerate(zip(fit_parameters, parameter_errors)):
                #print(f"Parameter {k}: {param} +/- {error}")
                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)][f"Parameter {k}"]={"value":param, "error":error}
            
            if do_plot:
            # Create a canvas for plotting
                c1 = ROOT.TCanvas("c1", "Fitted Histogram", 800, 600)

                # Plot the histogram
                hist.SetMarkerStyle(ROOT.kFullCircle)
                hist.Draw("ep")

                # Plot the fitted function
                polynomial.SetLineColor(ROOT.kRed)  # Set the color to red
                polynomial.Draw("same")

                # Display the canvas
                #c1.Update()
                #c1.Modified()

                # Optionally, save the plot as an image
                c1.SaveAs(f"./{pdir}pol_fits/eta{i}_nL{j}.png")

    #save the fit results as json file
    fit_results["bins"]={"abseta":list(abseta_bins), "nl": list(nl_bins)}


    json_object=json.dumps(fit_results, indent=4)
    with open(f"{hdir}/fit_results_res.json", "w") as outfile:
        outfile.write(json_object)  
    
    



def apply_res_corr(ntuples_gen, hdir, pdir, do_plot, do_binwise_plot=False):
    use_CB_smear=True

    pdir = pdir+'resolution/'
    #open fit_result file
    #open json, read as dict
    with open(f"{hdir}/fit_results_res.json", 'r') as openfile: fit_results = json.load(openfile)
    
    #get bins from
    abseta_bins=fit_results["bins"]["abseta"]
    nl_bins=fit_results["bins"]["nl"]

    #open data
    file=uproot.open(ntuples_gen)
    tree=file["Events"]
    variables=tree.keys()
    df=tree.arrays(variables, library="pd")

    df_new1=pd.DataFrame()

    #loop over bins, fill df_new1
    
    print("getting corrections for pt 1")



    for i in tqdm(range(len(abseta_bins)-1)):
        eta_1_filter=(abseta_bins[i]<df["eta_1"]) & (df["eta_1"]<=abseta_bins[i+1])
        df_e1=df[eta_1_filter]

        for j in range(len(nl_bins)-1):
            nl_1_filter=(nl_bins[j]<df_e1["nTrkLayers_1"]) & (df_e1["nTrkLayers_1"]<=nl_bins[j+1])
            df_en1=df_e1[nl_1_filter]

            if len(df_en1)>0:
                
                #get polynomial fit parameters
                poly_par=np.array([fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 0"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 2"]["value"]])
                pt1=df_en1["genpt_1"]
                #calc correction

                #calc std from parabula
                std1=pt1*pt1*poly_par[2]+pt1*poly_par[1]+poly_par[0]

                #smear pt
                if use_CB_smear:
                    #get CB sigma for current bin:
                    CB_sigma=fit_results["pull"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"]

                    df_en1["genpt_1_smeared"]=np.array(df_en1["pt_1"]) + np.random.normal(0, abs(std1), size=len(df_en1))*CB_sigma
                else:
                    df_en1["genpt_1_smeared"]=np.array(df_en1["pt_1"]) + np.random.normal(0, abs(std1), size=len(df_en1))#*CB_sigma

                df_new1=pd.concat([df_new1, df_en1])


                if do_plot and do_binwise_plot:
                    bins=50
                    rang=[0,200]
                    fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
                    plt.subplots_adjust(wspace=0, hspace=0)
                    ax0.set_xlim(left=rang[0],right=rang[1])
                    ax0.set_xticks([])
                    h_s=ax0.hist(df_en1["genpt_1_smeared"],bins=bins, range=rang, histtype="step",density=True, label="smeared")
                    h_r=ax0.hist(df_en1["pt_1"], bins=bins, histtype="step", range=rang, density=True, label="MC_reco")
                    ax0.legend()

                    if use_CB_smear:
                        ax0.set_title(f"with CB-smearing, CB_σ={round(CB_sigma,5)}")
                    else:
                        ax0.set_title("without CB-smearing")

                    ax1.set_ylim(bottom=0.85,top=1.15)
                    ax1.set_xlim(left=rang[0],right=rang[1])
                    ax1.plot(np.linspace(rang[0],rang[1],bins),h_s[0]/h_r[0],'.')
                    ax1.set_xlabel("pt_µ1 (GeV)")
                    ax1.set_ylabel("ratio")
                    ax1.grid(True)

                    if use_CB_smear:
                        plt.savefig(f"{pdir}binwise_pt_Distribution/abseta{i}_nL{j}_1_wCB.png")
                    else:
                        plt.savefig(f"{pdir}binwise_pt_Distribution/abseta{i}_nL{j}_1_woCB.png")
                    
                    plt.clf()
        

    df_new2=pd.DataFrame()
    df=df_new1
    print("getting corrections for pt 2")
    for i in tqdm(range(len(abseta_bins)-1)):
        eta_2_filter=(abseta_bins[i]<df["eta_2"]) & (df["eta_2"]<=abseta_bins[i+1])
        df_e2=df[eta_2_filter]

        for j in range(len(nl_bins)-1):
            nl_2_filter=(nl_bins[j]<df_e2["nTrkLayers_2"]) & (df_e2["nTrkLayers_2"]<=nl_bins[j+1])
            df_en2=df_e2[nl_2_filter]

            if len(df_en2)>0:
                
                #get polynomial fit parameters
                poly_par=np.array([fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 0"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 2"]["value"]] )
                
                pt2=df_en2["genpt_2"]
                #calc correction
                #calc std from parabula
                std2=pt2*pt2*poly_par[2]+pt2*poly_par[1]+poly_par[0]

                #smear pt 
                if use_CB_smear:
                    #get CB sigma for current bin:
                    CB_sigma=fit_results["pull"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"]
                    df_en2["genpt_2_smeared"]=np.array(df_en2["pt_2"]) + np.random.normal(0, abs(std2), size=len(df_en2))*CB_sigma 
                else:
                    df_en2["genpt_2_smeared"]=np.array(df_en2["pt_2"]) + np.random.normal(0, abs(std2), size=len(df_en2))
                df_new2=pd.concat([df_new2, df_en2])

                if do_plot and do_binwise_plot:
                    bins=50
                    rang=[0,200]
                    fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
                    ax0.set_xlim(left=rang[0],right=rang[1])
                    ax0.set_xticks([])
                    plt.subplots_adjust(wspace=0, hspace=0)
                    h_s=ax0.hist(df_en2["genpt_1_smeared"],bins=bins, range=rang, histtype="step",density=True, label="smeared")
                    h_r=ax0.hist(df_en2["pt_1"], bins=bins, histtype="step", range=rang, density=True, label="MC_reco")
                    ax0.legend()
                    if use_CB_smear:
                        ax0.set_title(f"with CB-smearing, CB_σ={round(CB_sigma,5)}")
                    else:
                        ax0.set_title("without CB-smearing")
                    
                    ax1.set_ylim(bottom=0.85,top=1.15)
                    ax1.set_xlim(left=rang[0],right=rang[1])
                    ax1.plot(np.linspace(rang[0],rang[1],bins),h_s[0]/h_r[0],'.')
                    ax1.set_xlabel("pt_µ2 (GeV/c)")
                    ax1.set_ylabel("ratio")
                    ax1.grid(True)

                    if use_CB_smear:
                        plt.savefig(f"{pdir}binwise_pt_Distribution/abseta{i}_nL{j}_2_wCB.png")
                    else:
                        plt.savefig(f"{pdir}binwise_pt_Distribution/abseta{i}_nL{j}_2_woCB.png")
                    plt.clf()
    

    #calculate invarriant mass
    
    df_new2["mass_Z_smeared_Dori"]=np.sqrt( 2*df_new2.genpt_1_smeared*df_new2.genpt_2_smeared*(np.cosh(df_new2.eta_1-df_new2.eta_2)-np.cos(df_new2.phi_1-df_new2.phi_2)) )

    if do_plot:
        bins=100
        rang=[50,130]
        fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
        plt.subplots_adjust(wspace=0, hspace=0)
        h_s=ax0.hist(df_new2["mass_Z_smeared_Dori"],bins=bins, range=rang, histtype="step",density=True, label="smeared")
        h_r=ax0.hist(df_new2["mass_Z"],bins=bins, histtype="step", range=rang, density=True, label="MC_reco")
        #ax0.hist(df_new2["mass_Z_smeared_Jost"],bins=bins, range=rang, histtype="step",density=True, label="MC_gen") 
        chi2=np.sum((h_r[0]-h_s[0])**2/h_s[0])
        ndf=bins
        ax0.annotate(f"χ²/NDF = {round(chi2,6)}/{ndf}",xy=(100,400), xycoords="figure pixels")
        ax0.set_xlim(left=rang[0],right=rang[1])
        ax0.set_xticks([])
        ax0.legend()
        ax0.set_xlabel("M_µµ (GeV)")
        ax0.set_title("without CB-smearing")
        ax1.set_ylim(bottom=0.85,top=1.15)
        ax1.set_xlim(left=rang[0],right=rang[1])
        ax1.plot(np.linspace(rang[0],rang[1],bins),h_s[0]/h_r[0])
        ax1.set_xlabel("M_µµ (GeV)")
        ax1.set_ylabel("ratio")
        ax1.grid(True)
        #plt.savefig(f"{pdir}Z_mass_comparison.pdf")
        if use_CB_smear:
            plt.savefig(f"{pdir}Z_mass_comparison_wCB.png")
        else:
            plt.savefig(f"{pdir}Z_mass_comparison_woCB.png")
        plt.clf()

    #save data
    print(f"saving corrected gen to {ntuples_gen.replace('.root', '_corr.root')}")
    data={key: df_new2[key].values for key in df_new2.columns}
    rdf = ROOT.RDF.MakeNumpyDataFrame(data)
    rdf.Snapshot("Events", ntuples_gen.replace('.root', '_corr.root'))
    print("done")

    



def make_CB_plot(hdir,pdir):
    with open(f"{hdir}/fit_results_res.json", 'r') as openfile: fit_results = json.load(openfile)
    abseta_bins=fit_results["bins"]["abseta"]
    nl_bins=fit_results["bins"]["nl"]

    CB_sigma=np.empty([len(abseta_bins)-1, len(nl_bins)-1])
    CB_sigma[:]=np.nan
    for i,I in enumerate(fit_results["pull"]):
        for j,J in enumerate(fit_results["pull"][I]):
            if not fit_results["pull"][I][J]=={}:
                CB_sigma[i][j]=fit_results["pull"][I][J]["Parameter 1"]["value"]
    
    im=plt.imshow(np.flip(CB_sigma.T,axis=0))
    plt.colorbar(im)
    plt.yticks([1,3,5,7,9],["16","14","12","10","8"])
    plt.ylabel("Number of layers")
    plt.xticks([-0.5,2,4.5,7,9.5],["0.0","0.5","1.0","1.5","2"])
    plt.xlabel("|eta|")
    plt.savefig(f"{pdir}sigma_CB.png")
    plt.clf()
