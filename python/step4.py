import numpy as np
from tqdm import tqdm
import json
import os
import ROOT
from array import array
from python.plot import roofit_mass, roofit_masscb, plot_ratio, plot_ratio2
from python.apply_corrections import step1, step2, step3, step4


def residual_correction(samples, abseta_bins, hdir, pdir, weight):

    # here the binning does not need to be changed for syst eval,
    # since only part of range is used in fit anyway
    m_bins = np.linspace(60, 120, 120)
    
    #get gen data
    df_gen = ROOT.RDataFrame("Events", samples["GEN"]["GEN"])
    df_gen = df_gen.Define("bs_weight", "(int)(poisson())")
    df_gen = df_gen.Define("weight", weight)
    df_gen = step1(df_gen, hdir, 'GEN')
    df_gen = step2(df_gen, hdir, 'GEN')
    df_gen = df_gen.Define("abseta_1", "abs(eta_1)").Define("abseta_2", "abs(eta_2)")

    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step2_fitresults.root", "read");')
    ROOT.gROOT.ProcessLine(f'TH3D* h_results_cb = (TH3D*)tf->Get("h_results_cb");')
    ROOT.gROOT.ProcessLine(f'TH3D* h_results_poly = (TH3D*)tf->Get("h_results_poly");')

    h_gen_mass_1 = df_gen.Histo2D(
        (
            'h_gen_mass_1', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_1",
        "genmass_Z", 
        "weight"
    )
    h_gen_mass_2 = df_gen.Histo2D(
        (
            'h_gen_mass_2', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_2",
        "genmass_Z", 
        "weight"
    )
    h_gen_mass_smeared_1 = df_gen.Histo2D(
        (
            'h_gen_mass_smeared_1', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_1",
        "genmass_Z_smeared", 
        "weight"
    )
    h_gen_mass_smeared_2 = df_gen.Histo2D(
        (
            'h_gen_mass_smeared_2', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_2",
        "genmass_Z_smeared", 
        "weight"
    )

    for typ in samples:
        if typ == "SIG" or typ=="DATA":
            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")
                
                #read data into dataframe
                df_reco = ROOT.RDataFrame("Events", samples[typ][subtyp])
                df_reco = df_reco.Define("bs_weight", "(int)(poisson())")
                df_reco = df_reco.Define("weight", weight)
                df_reco = step1(df_reco, hdir, typ)
                df_reco = step3(df_reco, hdir, typ)
                df_reco = df_reco.Define("abseta_1", "abs(eta_1)")
                df_reco = df_reco.Define("abseta_2", "abs(eta_2)")
                df_reco = df_reco.Filter("abseta_1 < 2.4 && abseta_2 < 2.4")

                h_mass_reco_1 = df_reco.Histo2D(
                    (
                        'h_mass_reco_1', '',
                        len(abseta_bins)-1, array('d', abseta_bins),
                        len(m_bins)-1, array('d', m_bins)
                    ),
                    "abseta_1",
                    "mass_Z_roccor_it", 
                    "weight"
                )
                h_mass_reco_2 = df_reco.Histo2D(
                    (
                        'h_mass_reco_2', '',
                        len(abseta_bins)-1, array('d', abseta_bins),
                        len(m_bins)-1, array('d', m_bins)
                    ),
                    "abseta_2",
                    "mass_Z_roccor_it", 
                    "weight"
                )

                tf = ROOT.TFile(f'{hdir}step4_{typ}.root', 'recreate')

                for eta in range(len(abseta_bins)-1):

                    h_gen_1 = h_gen_mass_1.ProjectionY('h_gen_1', eta+1, eta+1)
                    h_gen_2 = h_gen_mass_2.ProjectionY('h_gen_2', eta+1, eta+1)
                    h_gen = h_gen_1.Clone(f'h_gen_eta_{eta}')
                    h_gen.Add(h_gen_2)
                    h_gen.Scale(1./h_gen.Integral())

                    h_gen_smeared_1 = h_gen_mass_smeared_1.ProjectionY('h_gen_smeared_1', eta+1, eta+1)
                    h_gen_smeared_2 = h_gen_mass_smeared_2.ProjectionY('h_gen_smeared_2', eta+1, eta+1)
                    h_gen_smeared = h_gen_smeared_1.Clone(f'h_gen_smeared_eta_{eta}')
                    h_gen_smeared.Add(h_gen_smeared_2)
                    h_gen_smeared.Scale(1./h_gen_smeared.Integral())

                    h_reco_1 = h_mass_reco_1.ProjectionY('h_reco_1', eta+1, eta+1)
                    h_reco_2 = h_mass_reco_2.ProjectionY('h_reco_2', eta+1, eta+1)
                    h_reco = h_reco_1.Clone(f'h_reco_eta_{eta}')
                    h_reco.Add(h_reco_2)
                    h_reco.Scale(1./h_reco.Integral())

                    h_reco.Write()
                    h_gen.Write()
                    h_gen_smeared.Write()
                tf.Close()


def perform_fits(samples, abseta_bins, hdir, pdir, m_bins, doplot):

    massrange = [m_bins[0], m_bins[-1]]

    hists = []
    pdir = pdir+'residual/residual/'
    os.makedirs(pdir, exist_ok=True)
    for typ in samples:
        if typ == "DATA" or typ=='SIG':

            hists.append(
                ROOT.TH1D(
                    f'h_k_{typ}', '',
                    len(abseta_bins)-1, array('d', abseta_bins)
                )
            )
            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")  
                tf = ROOT.TFile(f'{hdir}step4_{typ}.root', 'read')
                os.makedirs(f'{pdir}template_fits/', exist_ok=True)

                for eta in range(len(abseta_bins)-1):
                    h_reco = tf.Get(f'h_reco_eta_{eta}')
                    h_gen_smeared = tf.Get(f'h_gen_smeared_eta_{eta}')
                    h_gen = tf.Get(f'h_gen_eta_{eta}')

                    print(h_gen, h_gen_smeared, h_reco)

                    if doplot:
                        genplot = f'{pdir}template_fits/gensmeared_{eta}'
                        recoplot = f'{pdir}template_fits/{typ}_{eta}'
                    else:
                        genplot=False
                        recoplot=False

                    results_gen, errors_gen = roofit_mass(h_gen, h_gen_smeared, plot=genplot, fitf='template', tag='', massrange=massrange)
                    # results_reco, errors_reco = roofit_mass(h_gen, h_reco, plot=recoplot, fitf='template', tag='', massrange=massrange)

                    # results_gen, errors_gen = roofit_masscb(h_gen, h_gen_smeared, plot=genplot, fitf='template', tag='', massrange=massrange)
                    results_reco, errors_reco = roofit_masscb(h_gen, h_reco, plot=recoplot, fitf='template', tag='', massrange=massrange)


                    hists[-1].SetBinContent(eta+1, results_reco[1]/results_gen[1])

                    print(results_gen[1], results_reco[1])
                tf.Close()
    
    tf = ROOT.TFile(f'{hdir}step4_k.root', 'recreate')
    for h in hists:
        h.Write()
    tf.Close()



def plot_closure(samples, hdir, pdir, weight, m_bins):
    pdir += 'residual/'
    hists = []

    for typ in tqdm(samples):
        first = 1
        for subtyp in samples[typ]:
            df = ROOT.RDataFrame("Events", samples[typ][subtyp])
            df = df.Define("bs_weight", "(int)(poisson())")
            df = df.Define("weight", weight)

            df = step1(df, hdir, typ)
            df = step2(df, hdir, typ)
            df = step3(df, hdir, typ)
            df = step4(df, hdir, typ)
            
            if typ == 'GEN':
                h_names = ['genmass_Z_smeared_SIG', 'genmass_Z_smeared_DATA', 'genmass_Z_smeared', "mass_Z_smeared"]
            else: 
                h_names = ['mass_Z_roccor_it']

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
    h_gen =  tf.Get('genmass_Z_smeared_GEN')
    h_gen.Sumw2()
    h_gen.Scale(1./h_gen.Integral())

    # MC
    h_mc = tf.Get(f'mass_Z_roccor_it_SIG')
    h_mc.Sumw2()
    h_mc.Scale(1./h_mc.Integral())

    h_gen_smeared = tf.Get(f'genmass_Z_smeared_SIG_GEN')
    h_gen_smeared.Sumw2()
    h_gen_smeared.Scale(1./h_gen_smeared.Integral())

    # first show impact of gensmearing to MC
    plot_ratio2(
        hists={
            'gen': h_gen_smeared,
            'mc': h_mc,
            'gen_smeared': h_gen
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_SIG',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.95, 1.05],
        labels={
            'gen': 'smeared genmass',
            'mc': 'reconstructed Mass',
            'gen_smeared': 'presmeared genmass'
        }
    )

    tf.Close()

    tf = ROOT.TFile(f'{hdir}step4_closure.root', 'read')

    # Then show impact of gensmearing to Data
    h_tmp_data = tf.Get(f'mass_Z_roccor_it_DATA')
    h_tmp_bkg = tf.Get('mass_Z_roccor_it_WW')
    h_tmp_bkg.Add(tf.Get('mass_Z_roccor_it_WZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_roccor_it_ZZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_roccor_it_TT'))
    h_tmp_mc = h_tmp_bkg.Clone('mass_Z_roccor_it_mc')
    h_tmp_mc.Add(tf.Get('mass_Z_roccor_it_SIG'))
    h_tmp_data.Scale(h_tmp_mc.Integral()/h_tmp_data.Integral())

    h_data = h_tmp_data.Clone('mass_Z_roccor_it_data')
    h_data.Add(h_tmp_bkg, -1)
    h_data.Sumw2()
    h_data.Scale(1./h_data.Integral())

    h_gen = tf.Get('genmass_Z_smeared_GEN')
    h_gen.Sumw2()
    h_gen.Scale(1./h_gen.Integral())

    h_gen_smeared = tf.Get(f'genmass_Z_smeared_DATA_GEN')
    h_gen_smeared.Sumw2()
    h_gen_smeared.Scale(1./h_gen_smeared.Integral())
   
    plot_ratio(
        hists={
            'gen': h_gen,
            'mc': h_gen_smeared,
            'dt': h_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_DATA',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.95, 1.05],
        labels={
            'gen': 'presmeared genmass',
            'mc': 'smeared genmass',
            'dt': 'data'
        }
    ) 

    # finally, show smeared reco mc vs data
    h_gen =  tf.Get('genmass_Z_smeared_GEN')
    h_gen.Sumw2()
    h_gen.Scale(1./h_gen.Integral())

    h_mc = tf.Get(f'mass_Z_smeared_GEN')
    h_mc.Sumw2()
    h_mc.Scale(1./h_mc.Integral())


    plot_ratio(
        hists={
            'gen': h_gen,
            'mc': h_mc,
            'dt': h_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_SIG_DATA',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.95, 1.05],
        labels={
            'gen': 'presmeared genmass',
            'mc': 'smeared reco Mass',
            'dt': 'data'
        }
    )