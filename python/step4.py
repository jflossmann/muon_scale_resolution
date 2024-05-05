import numpy as np
from tqdm import tqdm
import json
import os
import ROOT
from array import array
import glob
from python.plot import roofit_mass, roofit_masscb, plot_ratio, plot_ratio2
from python.apply_corrections import step1, step2, step3, step4


def load_df(samples, typ, weight, hdir, ):
    allfiles = []
    for s in samples[typ]:
        allfiles += glob.glob(samples[typ][s])

    df = ROOT.RDataFrame("Events", allfiles)
    df = df.Define("bs_weight", "(int)(poisson())")
    df = df.Define("weight", weight)

    df = step1(df, hdir, typ)
    df = step2(df, hdir, typ)
    if typ != 'GEN':
        df = step3(df, hdir, typ)
    
    df = df.Define("abseta_1", "abs(eta_1)")
    df = df.Define("abseta_2", "abs(eta_2)")
    df = df.Filter("abseta_1 < 2.4 && abseta_2 < 2.4")

    return df



def makeHists(df, abseta_bins, m_bins, rfile, typ, var, hname):
    tf = ROOT.TFile(rfile, 'update')
    h_mass_1 = df.Histo2D(
        (
            'h_mass_1', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_1",
        var, 
        "weight"
    )
    h_mass_2 = df.Histo2D(
        (
            'h_mass_2', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_2",
        var,
        "weight"
    )

    for eta in range(len(abseta_bins)-1):

        h_1 = h_mass_1.ProjectionY('h_1', eta+1, eta+1)
        h_2 = h_mass_2.ProjectionY('h_2', eta+1, eta+1)
        h_sum = h_1.Clone(f'{hname}_{typ}_eta_{eta}')
        h_sum.Add(h_2)

        if typ == 'DATA':
            h_bkg = tf.Get(f"{hname}_BKG_eta_{eta}")
            h_sig = tf.Get(f"{hname}_SIG_eta_{eta}")

            h_sig.Add(h_bkg)
            scale = h_sum.Integral() / h_sig.Integral()

            h_sum.Add(h_bkg, -scale)

        h_sum.Write()

    tf.Close()


def initialize_ks(k_init, rfile, abseta_bins, typ):
    tf = ROOT.TFile(rfile, 'update')

    k_hist = ROOT.TH2D(
        f'k_hist_{typ}', '',
        len(abseta_bins)-1, array('d', abseta_bins),
        5, array('d', [-2.5, -1.5, -.5, .5, 1.5, 2.5]), # 0 is nom, +/- 1,2 are up/dn 
    )
    # initialize values
    for eta in range(len(abseta_bins)-1):
        for var in range(5):
                k_hist.SetBinContent(eta+1, var+1, k_init[var])
    k_hist.Write()
    tf.Close()


def iterate(i, df_gen, rfile_ks, rfile_hists, abseta_bins, m_bins, chi2s, pdir, doplot, typ):
    if doplot:
        os.makedirs(f'{pdir}testit_{i}/', exist_ok=True)

    ROOT.gROOT.ProcessLine(f'TFile* tf{typ}{i} = TFile::Open("{rfile_ks}", "READ");')
    ROOT.gROOT.ProcessLine(f'TH2D *k_hist{typ}{i} = (TH2D*) tf{typ}{i}->Get("k_hist_{typ}");')

    tf_k = ROOT.TFile(rfile_ks, 'read')
    k_hist = tf_k.Get(f"k_hist_{typ}")
    k_hist.SetDirectory(ROOT.nullptr)
    tf_k.Close()

    if i == 0:
        vars = range(5)
        for var in vars:
            chi2s[var] = [1000 for j in range(len(abseta_bins)-1)]
    else:
        vars = [1,3]

    for var in vars:
        
        df_gen = df_gen.Define(
            f"genpt_1_smeared_{typ}_it{i}_var{var}",
            f'double pt1;\
            double k1 = k_hist{typ}{i}->GetBinContent(k_hist{typ}{i}->GetXaxis()->FindBin(abs(eta_1)), {var+1});\
            pt1 = genpt_1 / (1 + k1 * (genpt_1/genpt_1_step2 - 1));\
            return pt1;'
        )

        df_gen = df_gen.Define(
            f"genpt_2_smeared_{typ}_it{i}_var{var}",
            f"double pt2;\
            double k2 = k_hist{typ}{i}->GetBinContent(k_hist{typ}{i}->GetXaxis()->FindBin(abs(eta_2)), {var+1});\
            pt2 = genpt_2 / (1 + k2 * (genpt_2/genpt_2_step2 - 1));\
            return pt2;"
        )

        df_gen = df_gen.Define(
            f"genmass_Z_smeared_{typ}_{i}_var{var}",
            f"sqrt(2 * genpt_1_smeared_{typ}_it{i}_var{var} * genpt_2_smeared_{typ}_it{i}_var{var} * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)))"
        )

        makeHists(df_gen, abseta_bins, m_bins, rfile_hists, 'GEN', f"genmass_Z_smeared_{typ}_{i}_var{var}", f"h_gen_{typ}_it{i}_var{var}")

        tf = ROOT.TFile(rfile_hists, 'read')
        for eta in range(len(abseta_bins)-1):
            # TODO scale histograms
            h_gen = tf.Get(f"h_gen_{typ}_it{i}_var{var}_GEN_eta_{eta}")
            h_gen.Scale(1./h_gen.Integral())

            h_reco = tf.Get(f"mass_Z_{typ}_eta_{eta}")
            h_reco.Scale(1./h_reco.Integral())

            chi2s[var][eta] = h_reco.Chi2Test(h_gen, "WW CHI2/NDF")

            if doplot:
                plot_ratio(
                    hists={
                        'gen': h_gen,
                        'mc': h_reco,
                        'dt': h_reco
                    }, 
                    title='',
                    outfile=f'{pdir}testit_{i}/{typ}_etabin_{eta}_var_{var}',
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

    tf = ROOT.TFile(rfile_ks, 'update')
    tf.Delete(k_hist.GetName()+';1')
    k_hist.Write()
    tf.Close()

    return chi2s



def residual_correction(samples, abseta_bins, hdir, pdir, weight, doplot):

    m_bins = np.linspace(86, 96, 41)

    # definition and recreation of root files
    rfile_hists = f'{hdir}step4_hists.root'
    rfile_ks = f'{hdir}step4_k.root'
    
    tf = ROOT.TFile(rfile_hists, 'recreate')
    tf = ROOT.TFile(rfile_ks, 'recreate')
    tf.Close()

    #read data into dataframe
    df_data = load_df(samples, 'DATA', weight, hdir)
    df_sig = load_df(samples, 'SIG', weight, hdir)
    df_bkg = load_df(samples, 'BKG', weight, hdir)
    df_gen = load_df(samples, 'GEN', weight, hdir)

    # make reco hists
    makeHists(df_bkg, abseta_bins, m_bins, rfile_hists, 'BKG', 'mass_Z_step3', 'mass_Z')
    makeHists(df_sig, abseta_bins, m_bins, rfile_hists, 'SIG', 'mass_Z_step3', 'mass_Z')
    makeHists(df_data, abseta_bins, m_bins, rfile_hists, 'DATA', 'mass_Z_step3', 'mass_Z')

    niterations = 5

    k_init = [0.9, 1.0, 1.1, 1.2, 1.3]

    for typ in ["DATA", "SIG"]:
        chi2s = {}
        initialize_ks(k_init, rfile_ks, abseta_bins, typ)

        for i in tqdm(range(niterations)):
            chi2s = iterate(i, df_gen, rfile_ks, rfile_hists, abseta_bins, m_bins, chi2s, pdir, doplot, typ)


# TODO: add residual k-factor for scale or replace by fit


def plot_closure(samples, hdir, pdir, weight, m_bins):
    pdir += 'residual/'
    os.makedirs(pdir, exist_ok=True)
    hists = []
    # """
    for typ in tqdm(samples):
        for subtyp in samples[typ]:
            df = ROOT.RDataFrame("Events", samples[typ][subtyp])
            # df = df.Filter("pt_1>26 && pt_2>26 && pt_1 < 500 && pt_2<500")
            df = df.Define("bs_weight", "(int)(poisson())")
            df = df.Define("weight", weight)

            df = step1(df, hdir, typ)
            df = step2(df, hdir, typ)
            df = step3(df, hdir, typ)
            df = step4(df, hdir, typ)
            
            if typ == 'GEN':
                h_names = ['genmass_Z_step4', "genmass_Z_step3", "genmass_Z", "mass_Z_step4"]
            else: 
                h_names = ['mass_Z_step4', 'mass_Z_step3', "mass_Z_step2"]

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
    # """
    tf = ROOT.TFile(f'{hdir}step4_closure.root', 'read')
    print(tf.ls())
    h_gen_step4 =  tf.Get('genmass_Z_step4_GEN')
    h_gen_step4.Sumw2()
    h_gen_step4.Scale(1./h_gen_step4.Integral())

    h_gen_step3 =  tf.Get('genmass_Z_step3_GEN')
    h_gen_step3.Sumw2()
    h_gen_step3.Scale(1./h_gen_step3.Integral())


    h_tmp_data = tf.Get(f'mass_Z_step4_DATA')
    print(h_tmp_data.Integral())
    h_tmp_bkg = tf.Get('mass_Z_step4_WW')
    print(h_tmp_bkg.Integral())
    h_tmp_bkg.Add(tf.Get('mass_Z_step4_WZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_step4_ZZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_step4_TT'))
    h_tmp_mc = h_tmp_bkg.Clone('mass_Z_step4_mc')
    h_tmp_sig = tf.Get('mass_Z_step4_SIG').Clone('mass_Z_tmp_sig')
    h_tmp_mc.Add(h_tmp_sig)
    h_tmp_data.Scale(h_tmp_mc.Integral()/h_tmp_data.Integral())

    h_data = h_tmp_data.Clone('mass_Z_step4_data')
    h_data.Add(h_tmp_bkg, -1)
    h_data.Sumw2()
    h_data.Scale(1./h_data.Integral())

    # MC

    h_mc_step4 = tf.Get(f'mass_Z_step4_SIG')
    h_mc_step4.Sumw2()
    h_mc_step4.Scale(1./h_mc_step4.Integral())

    h_mc_step3 = tf.Get(f'mass_Z_step3_SIG')
    h_mc_step3.Sumw2()
    h_mc_step3.Scale(1./h_mc_step3.Integral())

    # first show impact of gensmearing to MC
    plot_ratio(
        hists={
            'gen': h_gen_step3,
            'mc': h_gen_step4,
            'dt': h_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_GEN',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.9, 1.1],
        labels={
            'gen': 'presemeared genmass',
            'mc': 'k-smeared genmass',
            'dt': 'data'
        }
    )

    plot_ratio(
        hists={
            'gen': h_mc_step3,
            'mc': h_mc_step4,
            'dt': h_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_RECO',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.9, 1.1],
        labels={
            'gen': 'reco-mass',
            'mc': 'k-smeared reco-mass',
            'dt': 'data'
        }
    )

    tf.Close()