import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import uproot
from tqdm import tqdm
import json
import os
import ROOT
from array import array
from python.plot import roofit_mass, plot_ratio

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

def chi2_2(k2, k1, z1, z2, df_gen, reco_hist, bins=50, rang=[81,101]):
    #reconstruct Z-Mass with smearing
    MZ=np.sqrt( 2*(1+k1*z1)*(1+k2*z2)*df_gen.genpt_1_smeared*df_gen.genpt_2_smeared*(np.cosh(df_gen.geneta_1-df_gen.geneta_2)-np.cos(df_gen.genphi_1-df_gen.genphi_2)) )
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

    #dict for results
    k_results={}
    m_bins = np.linspace(60, 120, 120)
    step2=False
    
    n=8
    #get gen data
    df_gen = ROOT.RDataFrame("Events", samples["GEN"]["GEN"])
    # df_gen = df_gen.Define("abseta_1", "abs(eta_1)"). Define("abseta_2", "abs(eta_2)")
    # df_gen = df_gen.Filter("abseta_1 < 2.4 && abseta_2 < 2.4")

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
        "genmass_Z"
    )
    h_gen_mass_2 = df_gen.Histo2D(
        (
            'h_gen_mass_2', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_2",
        "genmass_Z"
    )
    h_gen_mass_smeared_1 = df_gen.Histo2D(
        (
            'h_gen_mass_smeared_1', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_1",
        "genmass_Z_smeared"
    )
    h_gen_mass_smeared_2 = df_gen.Histo2D(
        (
            'h_gen_mass_smeared_2', '',
            len(abseta_bins)-1, array('d', abseta_bins),
            len(m_bins)-1, array('d', m_bins)
        ),
        "abseta_2",
        "genmass_Z_smeared"
    )

    for typ in samples:
        if typ == "SIG" or typ=="DATA":
            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")
                
                #read data into dataframe
                df_reco = ROOT.RDataFrame("Events", samples[typ][subtyp])
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
                    "mass_Z_roccor_it"
                )
                h_mass_reco_2 = df_reco.Histo2D(
                    (
                        'h_mass_reco_1', '',
                        len(abseta_bins)-1, array('d', abseta_bins),
                        len(m_bins)-1, array('d', m_bins)
                    ),
                    "abseta_2",
                    "mass_Z_roccor_it"
                )

                tf = ROOT.TFile(f'{hdir}step4_{typ}.root', 'recreate')

                for eta in range(len(abseta_bins)-1):

                    h_gen_1 = h_gen_mass_1.ProjectionY('h_gen_1', eta+1, eta+1)
                    h_gen_2 = h_gen_mass_1.ProjectionY('h_gen_2', eta+1, eta+1)
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

def perform_fits(samples, abseta_bins, hdir, pdir):

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

                    results_gen, errors_gen = roofit_mass(h_gen, h_gen_smeared, plot=f'{pdir}template_fits/gensmeared_{eta}', fitf='template', tag='')
                    results_reco, errors_reco = roofit_mass(h_gen, h_reco, plot=f'{pdir}template_fits/{typ}_{eta}', fitf='template', tag='')

                    hists[-1].SetBinContent(eta+1, results_reco[1]/results_gen[1])

                    print(results_gen[1], results_reco[1])
                tf.Close()
    
    tf = ROOT.TFile(f'{hdir}step4_k.root', 'recreate')
    for h in hists:
        h.Write()
    tf.Close()



def apply_res_corr(samples, hdir, pdir):
    pdir += 'residual/'
    m_bins = np.linspace(80, 102, 44)
    hists = []

    for typ in tqdm(samples):
        first = 1
        for subtyp in samples[typ]:
            df = ROOT.RDataFrame("Events", samples[typ][subtyp])
            
            if typ == 'GEN':
                ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step4_k.root", "read");')

                for dtsg in ['DATA', 'SIG']:

                    ROOT.gROOT.ProcessLine(f'TH1D* h_k_{dtsg} = (TH1D*)tf->Get("h_k_{dtsg}");')

                    df = df.Define(
                        f"genpt_1_smeared_{dtsg}", 
                        f"genpt_1 * (1 + h_k_{dtsg}->GetBinContent(h_k_{dtsg}->FindBin(abs(eta_1))) * (genpt_1_smeared/genpt_1 - 1))"
                    )
                    df = df.Define(
                        f"genpt_2_smeared_{dtsg}", 
                        f"genpt_2 * (1 + h_k_{dtsg}->GetBinContent(h_k_{dtsg}->FindBin(abs(eta_2))) * (genpt_2_smeared/genpt_2 - 1))"
                    )
                    df = df.Define(
                        f"genmass_Z_smeared_{dtsg}",
                        f"sqrt(2 * genpt_1_smeared_{dtsg} * genpt_2_smeared_{dtsg} * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)));"
                    )

                h_names = ['genmass_Z_smeared_SIG', 'genmass_Z_smeared_DATA', 'genmass_Z_smeared']

            else: 
                h_names = ['mass_Z_roccor_it']

            df = df.Define("weight", "zPtWeight*genWeight/sumwWeight*xsec*sf_id*sf_iso")

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

    h_dt = tf.Get(f'genmass_Z_smeared_SIG_GEN')
    h_dt.Sumw2()
    h_dt.Scale(1./h_dt.Integral())


    plot_ratio(
        hists={
            'gen': h_gen,
            'mc': h_mc,
            'dt': h_dt
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_SIG',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.85, 1.15],
        labels={
            'gen': 'presmeared genmass',
            'mc': 'reconstructed Mass',
            'dt': 'smeared genmass'
        }
    )

    tf.Close()

    tf = ROOT.TFile(f'{hdir}step4_closure.root', 'read')

    # DATA
    h_tmp_data = tf.Get(f'mass_Z_roccor_it_DATA')
    h_tmp_bkg = tf.Get('mass_Z_roccor_it_WW')
    h_tmp_bkg.Add(tf.Get('mass_Z_roccor_it_WZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_roccor_it_ZZ'))
    h_tmp_bkg.Add(tf.Get('mass_Z_roccor_it_TT'))
    h_tmp_mc = h_tmp_bkg.Clone('mass_Z_roccor_it_mc')
    h_tmp_mc.Add(tf.Get('mass_Z_roccor_it_SIG'))
    h_tmp_data.Scale(h_tmp_mc.Integral()/h_tmp_data.Integral())

   
    plot_ratio(
        hists={
            'gen': h_tmp_bkg,
            'mc': h_tmp_mc,
            'dt': h_tmp_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_test',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.85, 1.15],
        labels={
            'gen': 'background',
            'mc': 'MC total',
            'dt': 'data'
        }
    ) 

    h_gen =  tf.Get('genmass_Z_smeared_GEN')
    h_gen.Sumw2()
    h_gen.Scale(1./h_gen.Integral())

    h_data = h_tmp_data.Clone('h_DATA-BKG')
    h_data.Add(h_tmp_bkg, -1)
    h_data.Sumw2()
    h_data.Scale(1./h_data.Integral())

    h_dt = tf.Get(f'genmass_Z_smeared_DATA_GEN')
    h_dt.Sumw2()
    h_dt.Scale(1./h_dt.Integral())


    plot_ratio(
        hists={
            'gen': h_gen,
            'mc': h_data,
            'dt': h_dt
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_DATA',
        text=['','',''],
        #xrange=[80, 102],
        ratio_range = [0.85, 1.15],
        labels={
            'gen': 'presmeared genmass',
            'mc': 'reconstructed Mass',
            'dt': 'smeared genmass'
        }
    )