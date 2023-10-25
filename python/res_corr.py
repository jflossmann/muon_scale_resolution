import uproot
#import warnings
#warnings.simplefilter(action='ignore', category=UserWarning)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ROOT



def get_res_correction(ntuples_gen, pull_bins, abseta_bins, nl_bins, pt_bins, pdir ,do_plot):
    #read data
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

    #iterate over bins to fill histogram
    fit_par_pull=[]
    hist_std=[]

    #iterate over eta bins
    for i in tqdm(range(len(abseta_bins)-1)):
        
        hist_std.append([])
        fit_par_pull.append([])

        eta_1_filter=(abseta_bins[i]<df["eta_1"]) & (df["eta_1"]<=abseta_bins[i+1])
        eta_2_filter=(abseta_bins[i]<df["eta_2"]) & (df["eta_2"]<=abseta_bins[i+1])

        df_e1=df[eta_1_filter]
        df_e2=df[eta_2_filter]

        #iterate over pt bins
        for j in range(len(nl_bins)-1):
            hist_std[i].append([])


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

            # Define the Crystal Ball function
            crystal_ball = ROOT.TF1("crystal_ball", "[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4])", pull_bins[0], pull_bins[-1])
            # Set initial parameter values
            crystal_ball.SetParameter(0, 1000)  # Amplitude
            crystal_ball.SetParameter(1, 0.0)  # Mean
            crystal_ball.SetParameter(2, std)  # Sigma
            crystal_ball.SetParameter(3, 1.0)  # Alpha
            crystal_ball.SetParameter(4, 1.0)  # N

            # Fit the histogram with the Crystal Ball function
            hist.Fit("crystal_ball")
            fit_parameters = [crystal_ball.GetParameter(i) for i in range(5)]
            parameter_errors = [crystal_ball.GetParError(i) for i in range(5)]

            # Print the fit results
            for k, (param, error) in enumerate(zip(fit_parameters, parameter_errors)):
                print(f"Parameter {i}: {param} +/- {error}")
            
            if do_plot:
                # Create a canvas for plotting
                c1 = ROOT.TCanvas("c1", "Fitted Histogram", 800, 600)

                # Plot the histogram
                hist.Draw()

                # Plot the fitted function
                crystal_ball.SetLineColor(ROOT.kRed)  # Set the color to red
                crystal_ball.Draw("same")

                # Display the canvas
                c1.Update()
                c1.Modified()

                # Optionally, save the plot as an image
                c1.SaveAs(f"./{pdir}/CB/eta{i}_nL{j}_pull_fit.png")
                    
            
            #bin in pt_reco for polynomial fit
            #for k in range(len(pt_bins)-1):
#                pt_1_filter=(pt_bins[k]<df_en1["pt_1"]) & (df_en1["pt_1"]<=pt_bins[k+1])
 #               pt_2_filter=(pt_bins[k]<df_en2["pt_2"]) & (df_en2["pt_2"]<=pt_bins[k+1])
#
 #               df_enp1=df_en1[pt_1_filter]
  #              df_enp2=df_en2[pt_2_filter]
#
 #               R_1=df_enp1["R_1"]
  #              R_2=df_enp2["R_2"]
   #             R=pd.concat([R_1,R_2])
    #            hist_std[i][j].append(np.std(R))
     #       
      #      poly_hist=np.histogram(pull, bins=pull_bins)[0]
       #     hist = ROOT.TH1D("hist", "PolyHistogram", len(pt_bins) - 1, pt_bins)
        #    for k, value in enumerate(pull_hist):
         #       hist.SetBinContent(k+1, value)
#
 #           # Define the second-order polynomial function
  #          polynomial = ROOT.TF1("polynomial", "[0] + [1]*x + [2]*x*x")
#
 #           # Set initial parameter values
  #          polynomial.SetParameter(0, 10.0)  # Coefficient for the constant term
   #         polynomial.SetParameter(1, 2.0)   # Coefficient for the linear term
    #        polynomial.SetParameter(2, 0.5)   # Coefficient for the quadratic term
#
 #           # Fit the histogram with the polynomial function
  #          hist.Fit("polynomial")
   #         # Get the fit results
    #        fit_parameters = [polynomial.GetParameter(i) for i in range(3)]
     #       parameter_errors = [polynomial.GetParError(i) for i in range(3)]
#
 #           # Print the fit results
  #          for k, (param, error) in enumerate(zip(fit_parameters, parameter_errors)):
   #             print(f"Parameter {i}: {param} +/- {error}")
    #        
     #       if do_plot:
      #      # Create a canvas for plotting
       #         c1 = ROOT.TCanvas("c1", "Fitted Histogram", 800, 600)
#
 #               # Plot the histogram
  #              hist.Draw()
#
 #               # Plot the fitted function
  #              polynomial.SetLineColor(ROOT.kRed)  # Set the color to red
   #             polynomial.Draw("same")
#
 #               # Display the canvas
  #              c1.Update()
   #             c1.Modified()
#
 #               # Optionally, save the plot as an image
  #              c1.SaveAs(f"./{pdir}/Poly/histogram_fit.png")
