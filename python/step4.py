import numpy as np
from tqdm import tqdm
import json
import os
import ROOT
from array import array
from python.plot import roofit_mass, roofit_masscb, plot_ratio, plot_ratio2
from python.apply_corrections import step1, step2, step3, step4



def residual_correction(samples, abseta_bins, hdir, pdir, weight, doplot):

    m_bins = np.linspace(86, 96, 41)

    #read data into dataframe
    df_data = ROOT.RDataFrame("Events", samples['DATA']['DATA'])
    df_data = df_data.Filter("pt_1>26 && pt_2>26 && pt_1 < 500 && pt_2<500")
    df_data = df_data.Define("bs_weight", "(int)(poisson())")
    df_data = df_data.Define("weight", weight)
    df_data = step1(df_data, hdir, 'DATA')
    df_data = step2(df_data, hdir, 'DATA')
    df_data = step3(df_data, hdir, 'DATA')
    df_data = df_data.Define("abseta_1", "abs(eta_1)")
    df_data = df_data.Define("abseta_2", "abs(eta_2)")
    df_data = df_data.Filter("abseta_1 < 2.4 && abseta_2 < 2.4")

    h_mass_data_1 = df_data.Histo2D(
        (
            'h_mass_data_1', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_1",
        "mass_Z_step3", 
        "weight"
    )
    h_mass_data_2 = df_data.Histo2D(
        (
            'h_mass_data_2', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_2",
        "mass_Z_step3", 
        "weight"
    )

    df_gen = ROOT.RDataFrame("Events", samples["GEN"]["GEN"])
    df_gen = df_gen.Filter("pt_1>26 && pt_2>26 && pt_1 < 500 && pt_2<500")
    df_gen = df_gen.Define("bs_weight", "(int)(poisson())")
    df_gen = df_gen.Define("weight", weight)
    df_gen = step1(df_gen, hdir, 'GEN')
    df_gen = step2(df_gen, hdir, 'GEN')
    df_gen = df_gen.Define("abseta_1", "abs(eta_1)").Define("abseta_2", "abs(eta_2)")

    k_hist = ROOT.TH2D(
        'k_hist', '',
        len(abseta_bins)-1, array('d', abseta_bins),
        5, array('d', [-2.5, -1.5, -.5, .5, 1.5, 2.5]) # 0 is nom, +/- 1,2 are up/dn 
    )

    k_init = [0.9, 1.0, 1.1, 1.2, 1.3]

    # initialize values
    for eta in range(len(abseta_bins)-1):
        for var in range(5):
            k_hist.SetBinContent(eta+1, var+1, k_init[var])

    tf = ROOT.TFile(f'{hdir}step4_k_it.root', 'recreate')
    k_hist.Write()
    tf.Close()

    niterations = 5

    chi2s = {}

    for i in tqdm(range(niterations)):
        os.makedirs(f'{pdir}testit_{i}/', exist_ok=True)

        ROOT.gROOT.ProcessLine(f'TFile* tf{i} = TFile::Open("{hdir}step4_k_it.root", "READ");')
        ROOT.gROOT.ProcessLine(f'TH2D *k_hist{i} = (TH2D*) tf{i}->Get("k_hist");')

        if i == 0:
            vars = range(5)
            for var in vars:
                chi2s[var] = [1000 for j in range(len(abseta_bins)-1)]
        else:
            vars = [1,3]


        for var in vars:
            
            df_gen = df_gen.Define(
                f"genpt_1_smeared_it{i}_var{var}",
                f'double pt1;\
                Int_t etabin1 = h_results_cb->GetXaxis()->FindBin(abs(eta_1));\
                Int_t nlbin1 = h_results_cb->GetYaxis()->FindBin(nTrkLayers_1);\
                double mean_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 1);\
                double sig_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 2);\
                double n_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 3);\
                double alpha_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 4);\
                double k1 = k_hist{i}->GetBinContent(k_hist{i}->GetXaxis()->FindBin(abs(eta_1)), {var+1});\
                double rndm_cb1 = (double) cb_rndm(mean_cb1, k1 * sig_cb1, alpha_cb1, n_cb1);\
                double sig_poly1_a = h_results_poly->GetBinContent(etabin1, nlbin1, 1);\
                double sig_poly1_b = h_results_poly->GetBinContent(etabin1, nlbin1, 2);\
                double sig_poly1_c = h_results_poly->GetBinContent(etabin1, nlbin1, 3);\
                double sig_poly1 = sig_poly1_a + sig_poly1_b * genpt_1 + sig_poly1_c * genpt_1*genpt_1;\
                if (sig_poly1 < 0) sig_poly1 = 0;\
                pt1 = 1. / (1./genpt_1 * ( 1 + sig_poly1 * rndm_cb1));\
                return pt1;'
            )

            df_gen = df_gen.Define(
                f"genpt_2_smeared_it{i}_var{var}",
                f"double pt2;\
                Int_t etabin2 = h_results_cb->GetXaxis()->FindBin(abs(eta_2));\
                Int_t nlbin2 = h_results_cb->GetYaxis()->FindBin(nTrkLayers_2);\
                double mean_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 1);\
                double sig_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 2);\
                double n_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 3);\
                double alpha_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 4);\
                double k2 = k_hist{i}->GetBinContent(k_hist{i}->GetXaxis()->FindBin(abs(eta_2)), {var+1});\
                double rndm_cb2 = (double) cb_rndm(mean_cb2, k2*sig_cb2, alpha_cb2, n_cb2);\
                double sig_poly2_a = h_results_poly->GetBinContent(etabin2, nlbin2, 1);\
                double sig_poly2_b = h_results_poly->GetBinContent(etabin2, nlbin2, 2);\
                double sig_poly2_c = h_results_poly->GetBinContent(etabin2, nlbin2, 3);\
                double sig_poly2 = sig_poly2_a + sig_poly2_b * genpt_2 + sig_poly2_c * genpt_2*genpt_2;\
                if (sig_poly2 < 0) sig_poly2 = 0;\
                pt2 = 1. / (1./genpt_2 * ( 1 + sig_poly2 * rndm_cb2));\
                return pt2;"
            )

            df_gen = df_gen.Define(
                f"genmass_Z_smeared_{i}_var{var}",
                f"sqrt(2 * genpt_1_smeared_it{i}_var{var} * genpt_2_smeared_it{i}_var{var} * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)))"
            )
        
            # Create histograms
            h_gen_mass_smeared_1 = df_gen.Histo2D(
                (
                    f'h_gen_mass_smeared_1_k_{i}_{var}', '',
                    len(abseta_bins) - 1, array('d', abseta_bins),
                    len(m_bins) - 1, array('d', m_bins)
                ),
                "abseta_1",
                f"genmass_Z_smeared_{i}_var{var}",
                "weight"
            )

            h_gen_mass_smeared_2 = df_gen.Histo2D(
                (
                    f'h_gen_mass_smeared_2_k_{i}_var{var}', '',
                    len(abseta_bins) - 1, array('d', abseta_bins),
                    len(m_bins) - 1, array('d', m_bins)
                ),
                "abseta_2",
                f"genmass_Z_smeared_{i}_var{var}",
                "weight"
            )

            for eta in range(len(abseta_bins)-1):

                h_data_1 = h_mass_data_1.ProjectionY('h_reco_1', eta+1, eta+1)
                h_data_2 = h_mass_data_2.ProjectionY('h_reco_2', eta+1, eta+1)
                h_data = h_data_1.Clone(f'h_data_eta_{eta}')
                h_data.Add(h_data_2)
                h_data.Scale(1./h_data.Integral())

                h_gen_1 = h_gen_mass_smeared_1.ProjectionY(f'h_gen_1_{i}', eta+1, eta+1)
                h_gen_2 = h_gen_mass_smeared_2.ProjectionY(f'h_gen_2_{i}', eta+1, eta+1)
                h_gen = h_gen_1.Clone(f'h_gen_{i}_eta_{eta}')
                h_gen.Add(h_gen_2)
                h_gen.Scale(1./h_gen.Integral())

                chi2s[var][eta] = h_data.Chi2Test(h_gen, "WW CHI2/NDF")

                if doplot:
                    plot_ratio(
                        hists={
                            'gen': h_gen,
                            'mc': h_data,
                            'dt': h_data
                        }, 
                        title='',
                        outfile=f'{pdir}testit_{i}/etabin_{eta}_var_{var}',
                        text=[f'k: {round(k_hist.GetBinContent(eta+1, var+1),5)}','',''],
                        #xrange=[80, 102],
                        ratio_range = [0.95, 1.05],
                        labels={
                            'gen': 'smeared genmass',
                            'mc': 'data',
                            'dt': 'data',
                        }
                    )

        print(chi2s)

        for eta in range(len(abseta_bins)-1):
            if (chi2s[3][eta] > chi2s[0][eta]): # cut away two rightmost values 
                k_hist.SetBinContent(eta+1, 5, k_hist.GetBinContent(eta+1, 3))
                k_hist.SetBinContent(eta+1, 3, k_hist.GetBinContent(eta+1, 2))
                chi2s[4][eta] = chi2s[2][eta]
                chi2s[2][eta] = chi2s[1][eta]


            elif chi2s[1][eta] > chi2s[4][eta]: # cut away two leftmost values 
                k_hist.SetBinContent(eta+1, 1, k_hist.GetBinContent(eta+1, 3))
                k_hist.SetBinContent(eta+1, 3, k_hist.GetBinContent(eta+1, 4))
                chi2s[0][eta] = chi2s[2][eta]
                chi2s[2][eta] = chi2s[3][eta]

            
            else: # cut away right- and leftmost values 
                k_hist.SetBinContent(eta+1, 1, k_hist.GetBinContent(eta+1, 2))
                k_hist.SetBinContent(eta+1, 5, k_hist.GetBinContent(eta+1, 4))   
                chi2s[4][eta] = chi2s[3][eta]
                chi2s[0][eta] = chi2s[1][eta]
 

            k_hist.SetBinContent(eta+1, 2, 0.5 * (k_hist.GetBinContent(eta+1, 1) + k_hist.GetBinContent(eta+1, 3)))
            k_hist.SetBinContent(eta+1, 4, 0.5 * (k_hist.GetBinContent(eta+1, 3) + k_hist.GetBinContent(eta+1, 5)))

        tf = ROOT.TFile(f'{hdir}step4_k_it.root', 'recreate')
        k_hist.Write()
        tf.Close()




def plot_closure(samples, hdir, pdir, weight, m_bins):
    pdir += 'residual/'
    hists = []
    for typ in tqdm(samples):
        first = 1
        for subtyp in samples[typ]:
            df = ROOT.RDataFrame("Events", samples[typ][subtyp])
            df = df.Filter("pt_1>26 && pt_2>26 && pt_1 < 500 && pt_2<500")
            df = df.Define("bs_weight", "(int)(poisson())")
            df = df.Define("weight", weight)

            df = step1(df, hdir, typ)
            df = step2(df, hdir, typ)
            df = step3(df, hdir, typ)
            df = step4(df, hdir, typ)
            
            if typ == 'GEN':
                h_names = ['genmass_Z_step4', "genmass_Z_step3", "mass_Z_step4"]
            else: 
                h_names = ['mass_Z_step4']

            for h in h_names:
                hists.append(
                    df.Histo1D(
                        (
                            f'{h}_{subtyp}', '',
                            len(m_bins)-1, array('d', m_bins),
                        ),
                        h,
                        "weight"
                    )
                )

    tf = ROOT.TFile(f'{hdir}step4_closure.root', 'recreate')
    for h in hists:
        h.Write()
    tf.Close()

    tf = ROOT.TFile(f'{hdir}step4_closure.root', 'read')
    print(tf.ls())
    h_gen =  tf.Get('genmass_Z_step4_GEN')
    h_gen.Sumw2()
    h_gen.Scale(1./h_gen.Integral())


    h_tmp_data = tf.Get(f'mass_Z_step4_DATA')
    print(h_tmp_data.Integral())
    h_tmp_bkg = tf.Get('mass_Z_step4_WW')
    print(h_tmp_bkg.Integral())
    h_tmp_bkg.Add(tf.Get('mass_Z_step4_WZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_step4_ZZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_step4_TT'))
    h_tmp_mc = h_tmp_bkg.Clone('mass_Z_step4_mc')
    h_tmp_mc.Add(tf.Get('mass_Z_step4_SIG'))
    h_tmp_data.Scale(h_tmp_mc.Integral()/h_tmp_data.Integral())

    h_data = h_tmp_data.Clone('mass_Z_step4_data')
    h_data.Add(h_tmp_bkg, -1)
    h_data.Sumw2()
    h_data.Scale(1./h_data.Integral())

    # MC
    h_mc = tf.Get(f'mass_Z_step4_SIG')
    h_mc.Sumw2()
    h_mc.Scale(1./h_mc.Integral())
    # first show impact of gensmearing to MC
    plot_ratio(
        hists={
            'gen': h_mc,
            'mc': h_gen,
            'dt': h_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_SIG',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.95, 1.05],
        labels={
            'gen': 'k-smeared genmass',
            'mc': 'k-smeared reco-mass',
            'dt': 'data'
        }
    )

    tf.Close()