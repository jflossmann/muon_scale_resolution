import numpy as np
from tqdm import tqdm
import json
import os
import ROOT
from array import array
import glob
from python.plot import roofit_mass, roofit_masscb, plot_ratio, plot_ratio2
from python.apply_corrections import step1, step2, step3, step4, step5


def load_df(samples, typ, weight, hdir, ):
    allfiles = []
    for s in samples[typ]:
        allfiles += glob.glob(samples[typ][s])

    df = ROOT.RDataFrame("Events", allfiles)
    df = df.Define("bs_weight", "(int)(poisson())")
    df = df.Define("weight", weight)

    df = step1(df, hdir, typ)
    df = step2(df, hdir, typ)
    df = step3(df, hdir, typ)
    df = step4(df, hdir, typ)
    
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


def residual_correction(samples, abseta_bins, hdir, pdir, weight, doplot):

    pdir = f'{pdir}step5/'
    os.makedirs(pdir, exist_ok=True)

    m_bins = np.linspace(86, 96, 41)

    # definition and recreation of root files
    rfile_hists = f'{hdir}step5_hists.root'
    rfile_corr = f'{hdir}step5_corrections.root'
    
    tf = ROOT.TFile(rfile_hists, 'recreate')
    tf = ROOT.TFile(rfile_corr, 'recreate')
    tf.Close()

    #read data into dataframe
    df_data = load_df(samples, 'DATA', weight, hdir)
    df_sig = load_df(samples, 'SIG', weight, hdir)
    df_bkg = load_df(samples, 'BKG', weight, hdir)
    df_gen = load_df(samples, 'GEN', weight, hdir)

    # make reco hists
    makeHists(df_bkg, abseta_bins, m_bins, rfile_hists, 'BKG', 'mass_Z_step4', 'mass_Z')
    makeHists(df_sig, abseta_bins, m_bins, rfile_hists, 'SIG', 'mass_Z_step4', 'mass_Z')
    makeHists(df_data, abseta_bins, m_bins, rfile_hists, 'DATA', 'mass_Z_step4', 'mass_Z')
    makeHists(df_gen, abseta_bins, m_bins, rfile_hists, 'GEN', 'genmass_Z_step4', 'mass_Z')
    makeHists(df_gen, abseta_bins, m_bins, rfile_hists, 'GEN', 'genmass_Z', 'mass_Z_original')

    tf = ROOT.TFile(rfile_hists, 'read')

    plot_dict = {
        "DATA": ["mass_Z"],
        "SIG": ["mass_Z"],
        "GEN": ["mass_Z", "mass_Z_original"],
    }

    for typ in plot_dict.keys():

        for hname in plot_dict[typ]:

            plotdir = f'{pdir}{typ}_{hname}/'
            os.makedirs(plotdir, exist_ok=True)

            hist_corr = ROOT.TH1D(
                f'means_{hname}_{typ}', '',
                len(abseta_bins)-1, array('d', abseta_bins)
            )

            for eta in range(len(abseta_bins)-1):
                # fit all histograms
                hist = tf.Get(f"{hname}_{typ}_eta_{eta}")
                hist.Scale(10000/hist.Integral()) # improves fit quality
                print(hist, f"mass_Z_{typ}_eta_{eta}")
                results, errors = roofit_masscb(hist, hist, plot=f'{plotdir}{eta}', fitf='bwxcb', tag='', massrange = [86, 96])

                hist_corr.SetBinContent(eta+1, results[0])

                

            tf_c = ROOT.TFile(rfile_corr, 'update')
            hist_corr.Write()
            tf_c.Close()





def plot_closure(samples, hdir, pdir, weight, m_bins):
    pdir += 'step5/'
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
            df = step5(df, hdir, typ)
            
            if typ == 'GEN':
                h_names = ['genmass_Z_step5', "genmass_Z_step4", "genmass_Z"]
            else: 
                h_names = ['mass_Z_step5', 'mass_Z_step4', 'mass_Z_step3', "mass_Z_step2"]

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

    tf = ROOT.TFile(f'{hdir}step5_closure.root', 'recreate')
    for h in hists:
        h.Write()
    tf.Close()
    # """
    tf = ROOT.TFile(f'{hdir}step5_closure.root', 'read')
    print(tf.ls())
    h_gen_step4 =  tf.Get('genmass_Z_step5_GEN')
    h_gen_step4.Sumw2()
    h_gen_step4.Scale(1./h_gen_step4.Integral())

    h_gen_step3 =  tf.Get('genmass_Z_step4_GEN')
    h_gen_step3.Sumw2()
    h_gen_step3.Scale(1./h_gen_step3.Integral())


    h_tmp_data = tf.Get(f'mass_Z_step5_DATA')
    print(h_tmp_data.Integral())
    h_tmp_bkg = tf.Get('mass_Z_step5_WW')
    print(h_tmp_bkg.Integral())
    h_tmp_bkg.Add(tf.Get('mass_Z_step5_WZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_step5_ZZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_step5_TT'))
    h_tmp_mc = h_tmp_bkg.Clone('mass_Z_step5_mc')
    h_tmp_sig = tf.Get('mass_Z_step5_SIG').Clone('mass_Z_tmp_sig')
    h_tmp_mc.Add(h_tmp_sig)
    h_tmp_data.Scale(h_tmp_mc.Integral()/h_tmp_data.Integral())

    h_data = h_tmp_data.Clone('mass_Z_step5_data')
    h_data.Add(h_tmp_bkg, -1)
    h_data.Sumw2()
    h_data.Scale(1./h_data.Integral())

    # MC

    h_mc_step4 = tf.Get(f'mass_Z_step5_SIG')
    h_mc_step4.Sumw2()
    h_mc_step4.Scale(1./h_mc_step4.Integral())

    h_mc_step3 = tf.Get(f'mass_Z_step4_SIG')
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
            'gen': 'gen step4',
            'mc': 'gen step5',
            'dt': 'data step5'
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
            'gen': 'reco step4',
            'mc': 'reco step5',
            'dt': 'data'
        }
    )

    tf.Close()