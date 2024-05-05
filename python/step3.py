import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ROOT
from array import array
from python.plot import plot_ratio
import os
from time import time
from python.apply_corrections import step1, step2
import glob


def load_df(samples, typ, hdir, weight):
    allfiles = []
    for s in samples[typ]:
        allfiles += glob.glob(samples[typ][s])

    df = ROOT.RDataFrame("Events", allfiles)
    df = df.Define("bs_weight", "(int)(poisson())")
    df = df.Define("weight", weight)

    df = step1(df, hdir, typ)
    df = step2(df, hdir, typ)

    if typ == 'GEN':
        df = df.Filter(f"genmass_Z_step2 > 86 && genmass_Z_step2 < 96")

    df = df.Filter("abs(eta_1) < 2.4 && abs(eta_2) < 2.4")
    df = df.Define('masspt_1', 'mass_Z_step2 * pt_1_step2')
    df = df.Define('masspt_2', 'mass_Z_step2 * pt_2_step2')

    return df


def make_hists(df, typ, rfile, eta_bins, phi_bins, mass_bins, variables):

    print(f"Making {variables} histograms for {typ} sample...")

    # print(df.GetColumnNames())

    h_n = df.Histo3D(
        (
            f'h_{variables[0].split("_step")[0]}_n_{typ}', '',
            len(eta_bins)-1, array('d', eta_bins),
            len(phi_bins)-1, array('d', phi_bins),
            len(mass_bins)-1, array('d', mass_bins)
        ),
        'eta_1', 
        'phi_1',
        variables[0],
        "weight"
    )
    h_p = df.Histo3D(
        (
            f'h_{variables[1].split("_step")[0]}_p_{typ}', '',
            len(eta_bins)-1, array('d', eta_bins),
            len(phi_bins)-1, array('d', phi_bins),
            len(mass_bins)-1, array('d', mass_bins)
        ),
        'eta_2', 
        'phi_2',
        variables[1],
        "weight"
    )

    means_n = h_n.Project3DProfile(option='yx')
    means_p = h_p.Project3DProfile(option='yx')

    means = ROOT.TH2D(
        f'h_{variables[2].split("_step")[0]}_{typ}', '',
        len(eta_bins)-1, array('d', eta_bins),
        len(phi_bins)-1, array('d', phi_bins)
    )
    for eta in range(len(eta_bins)-1):
        for phi in range(len(phi_bins)-1):
            means.SetBinContent(
                eta+1, phi+1, 
                .5*(means_n.GetBinContent(eta+1, phi+1) + means_p.GetBinContent(eta+1, phi+1))
                )

    tf = ROOT.TFile(rfile, 'update')
    h_n.Write()
    h_p.Write()
    means_n.Write()
    means_p.Write()
    means.Write()
    tf.Close()

    return


def normalize_backgrounds(rf_step3, eta_bins, phi_bins):
    tf = ROOT.TFile(rf_step3, 'update')

    for np in ['n', 'p']:
        h_bkg = tf.Get(f"h_mass_Z_{np}_BKG")
        h_sig = tf.Get(f"h_mass_Z_{np}_SIG")
        h_bkg_norm = tf.Get(f"h_means_mass_Z_BKG").Clone(f"bkg_normalization_{np}")

        for eta in range(len(eta_bins)-1):
            for phi in range(len(phi_bins)-1):
                n_bkg = h_bkg.ProjectionZ(
                    "bkg_z", eta+1, eta+1, phi+1, phi+1
                ).Integral()
                n_sig = h_sig.ProjectionZ(
                    "sig_z", eta+1, eta+1, phi+1, phi+1
                ).Integral()

                h_bkg_norm.SetBinContent(
                    eta+1, phi+1, n_bkg / (n_sig + n_bkg)
                )
                if eta + phi == 0:
                    print(n_bkg, n_sig)

        h_bkg_norm.Write()
    tf.Close()


    
def update_corrections(rf_step3, rf_step3_k, eta_bins, phi_bins, i, typ):
    tf_hists = ROOT.TFile(rf_step3, 'read')
    tf_k = ROOT.TFile(rf_step3_k, 'read')

    reco_means_n = tf_hists.Get(f"h_mass_Z_n_{typ}_pyx")
    reco_means_p = tf_hists.Get(f"h_mass_Z_p_{typ}_pyx")
    reco_means_mpt_n = tf_hists.Get(f"h_masspt_1_n_{typ}_pyx")
    reco_means_mpt_p = tf_hists.Get(f"h_masspt_2_p_{typ}_pyx")
    if typ == 'DATA':
        bkg_means_n = tf_hists.Get(f"h_mass_Z_n_BKG_pyx")
        bkg_means_p = tf_hists.Get(f"h_mass_Z_p_BKG_pyx")
        bkg_means_mpt_n = tf_hists.Get(f"h_masspt_1_n_BKG_pyx")
        bkg_means_mpt_p = tf_hists.Get(f"h_masspt_2_p_BKG_pyx")
        bkg_normalizations_n = tf_hists.Get(f"bkg_normalization_n")
        bkg_normalizations_p = tf_hists.Get(f"bkg_normalization_p")

    gen_means_n = tf_hists.Get(f"h_genmass_Z_n_GEN_pyx")
    gen_means_p = tf_hists.Get(f"h_genmass_Z_p_GEN_pyx")

    h_kappa = tf_k.Get(f"M_{typ}")
    h_kappa.SetDirectory(ROOT.nullptr)
    h_lambd = tf_k.Get(f"A_{typ}")
    h_lambd.SetDirectory(ROOT.nullptr)

    reco_means = tf_hists.Get(f"h_means_mass_{typ}")
    reco_means.SetDirectory(ROOT.nullptr)


    for eta in range(len(eta_bins)-1):
        for phi in range(len(phi_bins)-1):
            reco_mean_n = reco_means_n.GetBinContent(eta+1, phi+1)
            reco_mean_p = reco_means_p.GetBinContent(eta+1, phi+1)

            reco_mean_mpt_n = reco_means_mpt_n.GetBinContent(eta+1, phi+1)
            reco_mean_mpt_p = reco_means_mpt_p.GetBinContent(eta+1, phi+1)

            if typ=='DATA':
                bkg_mean_n = bkg_means_n.GetBinContent(eta+1, phi+1)
                bkg_mean_p = bkg_means_p.GetBinContent(eta+1, phi+1)
                bkg_mean_mpt_n = bkg_means_mpt_n.GetBinContent(eta+1, phi+1)
                bkg_mean_mpt_p = bkg_means_mpt_p.GetBinContent(eta+1, phi+1)

                bkg_norm_n = bkg_normalizations_n.GetBinContent(eta+1, phi+1)
                bkg_norm_p = bkg_normalizations_p.GetBinContent(eta+1, phi+1)
                reco_mean_n = (reco_mean_n - bkg_norm_n * bkg_mean_n) / (1 - bkg_norm_n)
                reco_mean_p = (reco_mean_p - bkg_norm_p * bkg_mean_p) / (1 - bkg_norm_p)
                reco_mean_mpt_n = (reco_mean_mpt_n - bkg_norm_n * bkg_mean_mpt_n) / (1 - bkg_norm_n)
                reco_mean_mpt_p = (reco_mean_mpt_p - bkg_norm_p * bkg_mean_mpt_p) / (1 - bkg_norm_p)

            gen_mean_n = gen_means_n.GetBinContent(eta+1, phi+1)
            gen_mean_p = gen_means_p.GetBinContent(eta+1, phi+1)

            reco_means.SetBinContent(
                i+1, eta+1, phi+1,
                .5*(reco_mean_n + reco_mean_p)
                )

            V_n = -2 * (gen_mean_n - reco_mean_n)
            V_p = -2 * (gen_mean_p - reco_mean_p)

            if eta==3 and phi == 3:
                print(V_n, V_p)

            M_n = 2 * reco_mean_n
            M_p = 2 * reco_mean_p

            K_n = - reco_mean_mpt_n
            K_p = reco_mean_mpt_p

            kappa = (V_p / K_p - V_n / K_n) / (M_p / K_p - M_n / K_n) + 1
            lambd = ((V_p - M_p * (kappa-1)) / K_p + (V_n - M_n * (kappa-1)) / K_n) / 2.

            kappa_it = h_kappa.GetBinContent(eta+1, phi+1) * kappa
            lambd_it = h_lambd.GetBinContent(eta+1, phi+1) * kappa + lambd

            h_kappa.SetBinContent(eta+1, phi+1, kappa_it)
            h_lambd.SetBinContent(eta+1, phi+1, lambd_it)

    tf_hists.Close()
    tf_k.Close()

    tf_hists = ROOT.TFile(rf_step3, 'update')
    reco_means.Write()
    tf_hists.Close()

    tf_k = ROOT.TFile(rf_step3_k, 'update')
    tf_k.Delete(h_kappa.GetName()+';1')
    tf_k.Delete(h_lambd.GetName()+';1')
    h_kappa.Write()
    h_lambd.Write()
    tf_k.Close()


def iterative_correction(samples, eta_bins, phi_bins, mass_bins, hdir, pdir, iterationsteps, weight):
    masspt_bins = np.linspace(0, 1e4, 1000)

    pdir += 'iterative/'

    rf_step3_hists = f'{hdir}step3_masshists.root'
    rf_step3_corr = f'{hdir}step3_correction.root'

    tf = ROOT.TFile(rf_step3_hists, 'recreate')
    tf = ROOT.TFile(rf_step3_corr, 'recreate')
    tf = ROOT.TFile(f'{hdir}step3_closure.root', 'recreate')
    tf.Close()

    # save gen hists
    gen_df = load_df(samples, 'GEN', hdir, weight)
    make_hists(gen_df, 'GEN', rf_step3_hists, eta_bins, phi_bins, mass_bins, ["genmass_Z_step2", "genmass_Z_step2", "means_genmass_Z_step2"])

    # get initial histograms for sig, bkg and data
    bkg_df = load_df(samples, 'BKG', hdir, weight)
    make_hists(bkg_df, 'BKG', rf_step3_hists, eta_bins, phi_bins, mass_bins, ["mass_Z_step2", "mass_Z_step2", "means_mass_Z_step2"])
    make_hists(bkg_df, 'BKG', rf_step3_hists, eta_bins, phi_bins, masspt_bins, ["masspt_1", "masspt_2", "means_masspt"])

    sig_df = load_df(samples, 'SIG', hdir, weight)
    make_hists(sig_df, 'SIG', rf_step3_hists, eta_bins, phi_bins, mass_bins, ["mass_Z_step2", "mass_Z_step2", "means_mass_Z_step2"])
    make_hists(sig_df, 'SIG', rf_step3_hists, eta_bins, phi_bins, masspt_bins, ["masspt_1", "masspt_2", "means_masspt"])

    data_df = load_df(samples, 'DATA', hdir, weight)
    make_hists(data_df, 'DATA', rf_step3_hists, eta_bins, phi_bins, mass_bins, ["mass_Z_step2", "mass_Z_step2", "means_mass_Z_step2"])
    make_hists(data_df, 'DATA', rf_step3_hists, eta_bins, phi_bins, masspt_bins, ["masspt_1", "masspt_2", "means_masspt"])

    dfs = {
        "SIG": sig_df,
        "DATA": data_df
    }

    normalize_backgrounds(rf_step3_hists, eta_bins, phi_bins)


    for typ in samples:
        if typ == "SIG" or typ=="DATA":
            mass = f'mass_Z_step3'

            df_reco = dfs[typ]
                
            reco_means = ROOT.TH3D(
                f'h_means_mass_{typ}', '',
                iterationsteps, array('d', np.linspace(-0.5, iterationsteps-0.5, iterationsteps+1)),
                len(eta_bins)-1, array('d', eta_bins),
                len(phi_bins)-1, array('d', phi_bins)
            )
            tf = ROOT.TFile(rf_step3_hists, 'update')
            reco_means.Write()

            tf.Close()

            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")

                tf = ROOT.TFile(f'{hdir}step1_C.root', 'read')
                h_kappa = tf.Get(f'M_{typ}')
                h_kappa.SetDirectory(ROOT.nullptr)
                h_lambd = tf.Get(f'A_{typ}')
                h_lambd.SetDirectory(ROOT.nullptr)
                tf.Close()
                tf_corr = ROOT.TFile(rf_step3_corr, 'update')
                h_kappa.Write()
                h_lambd.Write()
                tf_corr.Close()

                #read data into dataframe
                df_reco = df_reco.Define(
                    'pt_1_step3', 'pt_1_step2'
                ).Define(
                    'pt_2_step3', 'pt_2_step2'
                ).Define(
                    'mass_Z_step3', 'mass_Z_step2'
                )


                for i in tqdm(range(iterationsteps)):
                    df_reco = df_reco.Redefine(
                        f"mass_Z_step3",
                        f"sqrt( 2 * pt_1_step3 * pt_2_step3 * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)) )"
                    )
                    df_reco_f = df_reco.Filter(f'{mass} > 86 && {mass} < 96')
                    df_reco_f = df_reco_f.Redefine(f"masspt_1", f"{mass} * pt_1_step3")
                    df_reco_f = df_reco_f.Redefine(f"masspt_2", f"{mass} * pt_2_step3")

                    make_hists(df_reco_f, typ, rf_step3_hists, eta_bins, phi_bins, mass_bins, [mass, mass, f"mean_{mass}"])
                    make_hists(df_reco_f, typ, rf_step3_hists, eta_bins, phi_bins, masspt_bins, ["masspt_1", "masspt_2", "mean_masspt"])

                    update_corrections(rf_step3_hists, rf_step3_corr, eta_bins, phi_bins, i, typ)

                    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{rf_step3_corr}", "READ");')
                    ROOT.gROOT.ProcessLine(f'TH2D* h_kappa_{i}_{typ} = (TH2D*)tf->Get("M_{typ}");')
                    ROOT.gROOT.ProcessLine(f'h_kappa_{i}_{typ}->SetDirectory(nullptr);')
                    ROOT.gROOT.ProcessLine(f'TH2D* h_lambd_{i}_{typ} = (TH2D*)tf->Get("A_{typ}");')
                    ROOT.gROOT.ProcessLine(f'h_lambd_{i}_{typ}->SetDirectory(nullptr);')
                    ROOT.gROOT.ProcessLine(f'tf->Close();')

                    # application of corrections
                    df_reco = df_reco.Redefine(
                        f"pt_1_step3",
                        f"double pt;\
                        pt = 1./ (h_kappa_{i}_{typ}->GetBinContent( h_kappa_{i}_{typ}->GetXaxis()->FindBin(eta_1) , h_kappa_{i}_{typ}->GetYaxis()->FindBin(phi_1) ) / pt_1 - \
                        h_lambd_{i}_{typ}->GetBinContent( h_lambd_{i}_{typ}->GetXaxis()->FindBin(eta_1) , h_lambd_{i}_{typ}->GetYaxis()->FindBin(phi_1) ));\
                        return pt;"
                    )
                    df_reco = df_reco.Redefine(
                        f"pt_2_step3",
                        f"1./ (h_kappa_{i}_{typ}->GetBinContent( h_kappa_{i}_{typ}->GetXaxis()->FindBin(eta_2) , h_kappa_{i}_{typ}->GetYaxis()->FindBin(phi_2) ) / pt_2 + \
                        h_lambd_{i}_{typ}->GetBinContent( h_lambd_{i}_{typ}->GetXaxis()->FindBin(eta_2) , h_lambd_{i}_{typ}->GetYaxis()->FindBin(phi_2) ))"
                    )

                    df_reco_f = df_reco.Filter(f'{mass} > 86 && {mass} < 96')
                    df_reco_f = df_reco.Filter(f'abs(eta_1) < 2.4 && abs(eta_2) < 2.4')
            
            rang = np.linspace(86, 96, 60)
            h_mass_it = df_reco.Histo1D((f'h_mass_it_{typ}', '', len(rang)-1, array('d', rang)), 'mass_Z_step3', "weight")
            h_mass = df_reco.Histo1D((f'h_mass_{typ}', '', len(rang)-1, array('d', rang)), 'mass_Z_step2', "weight")
            h_mass_gen = gen_df.Histo1D(('h_mass_gen', '', len(rang)-1, array('d', rang)), 'genmass_Z_step2', "weight")
            h_mass_bkg = bkg_df.Histo1D(('h_mass_bkg', '', len(rang)-1, array('d', rang)), 'mass_Z_step2', "weight")

            tf = ROOT.TFile(f'{hdir}step3_closure.root', 'update')
            h_mass_it.Write()
            h_mass.Write()
            h_mass_gen.Write()
            h_mass_bkg.Write()
            tf.Close()


def plot_closure(hdir, pdir, samples, eta_bins, phi_bins, iterationsteps):
    pdir = pdir+'iterative/'
    os.makedirs(pdir, exist_ok=True)

    tf_closure = ROOT.TFile(f'{hdir}step3_closure.root', 'read')
    tf_means = ROOT.TFile(f'{hdir}step3_masshists.root', 'read')

    h_mass_sig_it = tf_closure.Get('h_mass_it_SIG')
    h_mass_sig_it.Sumw2()
    h_mass_sig_it_c = h_mass_sig_it.Clone()
    h_mass_sig_it.Scale(1./h_mass_sig_it.Integral())

    h_mass_data_it = tf_closure.Get('h_mass_it_DATA')
    h_mass_data_it.Sumw2()
    h_mass_bkg = tf_closure.Get('h_mass_bkg')
    h_mass_bkg.Sumw2()
    scale_todata = h_mass_data_it.Integral()/(h_mass_bkg.Integral()+h_mass_sig_it_c.Integral())
    h_mass_bkg.Scale(scale_todata)
    h_mass_data_it.Add(h_mass_bkg, -1)    
    h_mass_data_it.Scale(1./h_mass_data_it.Integral())

    h_mass_data = tf_closure.Get('h_mass_DATA')
    h_mass_data.Sumw2()
    h_mass_bkg = tf_closure.Get('h_mass_bkg')
    h_mass_bkg.Sumw2()
    scale_todata = h_mass_data.Integral()/(h_mass_bkg.Integral()+h_mass_sig_it_c.Integral())
    h_mass_bkg.Scale(scale_todata)
    h_mass_data.Add(h_mass_bkg, -1)    
    h_mass_data.Scale(1./h_mass_data.Integral())

    h_mass_sig = tf_closure.Get('h_mass_SIG')
    h_mass_sig.Sumw2()
    h_mass_sig.Scale(1./h_mass_sig.Integral())

    h_mass_gen = tf_closure.Get('h_mass_gen')
    h_mass_gen.Sumw2()
    h_mass_gen.Scale(1./h_mass_gen.Integral())


    data_means = tf_means.Get(f'h_means_mass_DATA')
    sig_means = tf_means.Get(f'h_means_mass_SIG')
    gen_means = tf_means.Get(f'h_means_genmass_Z_GEN')

    # before correction
    plot_ratio(
        hists={
            'gen': h_mass_gen,
            'mc': h_mass_sig,
            'dt': h_mass_data
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_before',
        text=['','',''],
        #xrange=[60, 120],
        labels={
            'gen': 'gen mc',
            'mc': 'reco mc',
            'dt': 'data'
        },
        ratio_range=[0.9, 1.1]
    )

    # after correction
    plot_ratio(
        hists={
            'gen': h_mass_gen,
            'mc': h_mass_sig_it,
            'dt': h_mass_data_it
        }, 
        title='',
        outfile=f'{pdir}Z_mass_comparison_after',
        text=['','',''],
        #xrange=[60, 120],
        labels={
            'gen': 'gen mc',
            'mc': 'it. corr. reco mc',
            'dt': 'it. corr. data'
        },
        ratio_range=[0.9, 1.1]
    )
                    

    for typ in samples:
        if typ =='DATA' or typ=='SIG':
            os.makedirs(f"{pdir}iterative/binwise/{typ}", exist_ok=True)
            #make binwise plot
            for i in range(len(eta_bins)-1):
                for j in range(len(phi_bins)-1):
                    reco=[]
                    gen=[]
                    iterations = []
                    for n in range(iterationsteps):
                        if typ == 'DATA':
                            reco.append(data_means.GetBinContent(n+1, i+1, j+1))
                        else:
                            reco.append(sig_means.GetBinContent(n+1, i+1, j+1))

                        gen.append(gen_means.GetBinContent(i+1, j+1))
                        iterations.append(n)

                    plt.plot(iterations, reco,"x",c="blue",label="RECO")
                    plt.plot(iterations, gen,"--",c="darkblue",label="GEN")
                    plt.xlabel("iteration")
                    plt.ylabel("mean(M_µµ)")
                    plt.title(f"eta[{eta_bins[i]}, {eta_bins[i+1]}], phi[{round(phi_bins[i],1)}, {round(phi_bins[i+1],1)}]")
                    plt.legend()
                    plt.savefig(f"{pdir}iterative/binwise/{typ}/{typ}iteration_mean_eta{i}_phi{j}.png")
                    plt.clf()