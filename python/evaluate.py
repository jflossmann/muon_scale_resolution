import ROOT
from python.apply_corrections import step1, step2, step3, step4, uncertainties
from python.plot import plot_ratio, plot_ratio_updn
import os
import numpy as np
import glob
from tqdm import tqdm
from array import array

def plot_closure(samples, hdir, pdir, weight, m_bins):
    print(samples)
    pdir += 'uncertainties/'
    m_bins = np.linspace(80, 102, 45)
    os.makedirs(pdir, exist_ok=True)
    hists = []
    # """
    for typ in tqdm(samples):
        allfiles = []
        for s in samples[typ]:
            allfiles += glob.glob(samples[typ][s])
        print(allfiles)
        df = ROOT.RDataFrame("Events", allfiles)

        df = df.Define("bs_weight", "(int)(poisson())")
        df = df.Define("weight", weight)

        df = step1(df, hdir, typ)
        df = step2(df, hdir, typ)
        df = step3(df, hdir, typ)
        df = step4(df, hdir, typ)
        df = uncertainties(df, hdir, typ)
            
        if typ == 'GEN':
            h_names = ['genmass_Z_step4']
        elif typ == 'SIG':
            h_names = ['mass_Z_step4', 'mass_Z_step4_scaleup', 'mass_Z_step4_scaledn', 'mass_Z_step4_resolup', 'mass_Z_step4_resoldn']
        else: 
            h_names = ['mass_Z_step4']

        for h in h_names:
            hists.append(
                df.Histo1D(
                    (
                        f'{h}_{typ}', '',
                        len(m_bins)-1, array('d', m_bins),
                    ),
                    h,
                    "weight"
                )
            )

    tf = ROOT.TFile(f'{hdir}uncertainties.root', 'recreate')
    for h in hists:
        h.Write()
    tf.Close()
    # """
    tf = ROOT.TFile(f'{hdir}uncertainties.root', 'read')

    # data distributions - bkg contributions
    h_data = tf.Get(f'mass_Z_step4_DATA')
    h_tmp_bkg = tf.Get('mass_Z_step4_BKG')
    h_tmp_sig = tf.Get('mass_Z_step4_SIG')
    scale = h_data.Integral() / (h_tmp_sig.Integral() + h_tmp_bkg.Integral())
    h_tmp_bkg.Scale(scale)

    h_data.Add(h_tmp_bkg, -1)
    h_data.Scale(1./h_data.Integral())

    # MC
    h_mc_step4 = tf.Get(f'mass_Z_step4_SIG')
    h_mc_step4.Scale(1./h_mc_step4.Integral())

    h_mc_step4_scaleup = tf.Get(f'mass_Z_step4_scaleup_SIG')
    h_mc_step4_scaleup.Scale(1./h_mc_step4_scaleup.Integral())

    h_mc_step4_scaledn = tf.Get("mass_Z_step4_scaledn_SIG")
    h_mc_step4_scaledn.Scale(1./h_mc_step4_scaledn.Integral())

    h_mc_step4_resolup = tf.Get(f'mass_Z_step4_resolup_SIG')
    h_mc_step4_resolup.Scale(1./h_mc_step4_resolup.Integral())

    h_mc_step4_resoldn = tf.Get("mass_Z_step4_resoldn_SIG")
    h_mc_step4_resoldn.Scale(1./h_mc_step4_resoldn.Integral())

    # first show impact of gensmearing to MC
    plot_ratio_updn(
        hists=[h_data.Clone(), h_mc_step4.Clone(), h_mc_step4_scaleup.Clone(), h_mc_step4_scaledn.Clone()],
        title='',
        outfile=f'{pdir}scale_updn.pdf',
        text=['','',''],
        xrange=[80, 102],
        ratiorange = [0.9, 1.1],
        labels=['data', 'mc', 'mc scale up', 'mc scale down']
    )

    plot_ratio_updn(
        hists=[h_data, h_mc_step4, h_mc_step4_resolup, h_mc_step4_resoldn],
        title='',
        outfile=f'{pdir}resol_updn.pdf',
        text=['','',''],
        xrange=[80, 102],
        ratiorange = [0.9, 1.1],
        labels=['data', 'mc', 'mc resolution up', 'mc resolution down']
    )

    tf.Close()