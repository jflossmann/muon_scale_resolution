import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ROOT
from array import array
from python.plot import plot_ratio
import os
from time import time
from python.apply_corrections import step1, step2

def iterative_correction(samples, eta_bins, phi_bins, mass_bins, hdir, pdir, iterationsteps, weight):
    masspt_bins = np.linspace(0, 1e4, 1000)
    q_bins = [-1, 0, 1]
    pdir += 'iterative/'

    print(samples['GEN']['GEN'])

    df_gen = ROOT.RDataFrame("Events", samples['GEN']['GEN'])
    df_gen = df_gen.Define("bs_weight", "(int)(poisson())")
    df_gen = df_gen.Define("weight", weight)
    df_gen = step1(df_gen, hdir, 'GEN')
    df_gen = step2(df_gen, hdir, 'GEN')
    df_gen = df_gen.Filter('genmass_Z_smeared > 86 && genmass_Z_smeared < 96')

    h_gen_n = df_gen.Histo3D(
        (
            'h_gen_n', '',
            len(eta_bins)-1, array('d', eta_bins),
            len(phi_bins)-1, array('d', phi_bins),
            len(mass_bins)-1, array('d', mass_bins)
        ),
        'eta_1',
        'phi_1',
        'genmass_Z_smeared',
        "weight"
    )
    h_gen_p = df_gen.Histo3D(
        (
            'h_gen_p', '',
            len(eta_bins)-1, array('d', eta_bins),
            len(phi_bins)-1, array('d', phi_bins),
            len(mass_bins)-1, array('d', mass_bins)
        ),
        'eta_2',
        'phi_2',
        'genmass_Z_smeared',
        "weight"
    )

    gen_means_n = h_gen_n.Project3DProfile(option='yx')
    gen_means_p = h_gen_p.Project3DProfile(option='yx')

    gen_means = ROOT.TH2D(
        f'gen_means', '',
        len(eta_bins)-1, array('d', eta_bins),
        len(phi_bins)-1, array('d', phi_bins)
    )
    for eta in range(len(eta_bins)-1):
        for phi in range(len(phi_bins)-1):
            gen_means.SetBinContent(
                eta+1, phi+1, 
                .5*(gen_means_n.GetBinContent(eta+1, phi+1) + gen_means_p.GetBinContent(eta+1, phi+1))
                )

    for typ in samples:
        if typ == "SIG" or typ=="DATA":
            # if typ == "SIG":

            reco_means = ROOT.TH3D(
                f'reco_means_{typ}', '',
                iterationsteps, array('d', np.linspace(-0.5, iterationsteps-0.5, iterationsteps+1)),
                len(eta_bins)-1, array('d', eta_bins),
                len(phi_bins)-1, array('d', phi_bins)
            )

            for subtyp in samples[typ]:
                print(f"now processing {subtyp}")

                tf = ROOT.TFile(f'{hdir}step1_C.root', 'read')
                h_kappa = tf.Get(f'M_{typ}')
                h_kappa.SetDirectory(ROOT.nullptr)
                h_lambd = tf.Get(f'A_{typ}')
                h_lambd.SetDirectory(ROOT.nullptr)
                tf.Close()

                #read data into dataframe
                df_reco = ROOT.RDataFrame("Events", samples[typ][subtyp])
                df_reco = df_reco.Define("bs_weight", "(int)(poisson())")
                df_reco = df_reco.Define("weight", weight)
                df_reco = step1(df_reco, hdir, typ)
                df_reco = df_reco.Define(f'mass_Z_roccor_it', 'mass_Z_roccor')
                df_reco = df_reco.Define('pt_1_roccor_it', 'pt_1_roccor')
                df_reco = df_reco.Define('pt_2_roccor_it', 'pt_2_roccor')
                df_reco = df_reco.Define('masspt_1', 'mass_Z_roccor * pt_1_roccor')
                df_reco = df_reco.Define('masspt_2', 'mass_Z_roccor * pt_2_roccor')


                mass = f'mass_Z_roccor_it'

                df_reco_f = df_reco.Filter(f'{mass} > 86 && {mass} < 96')
                df_reco_f = df_reco.Filter(f'abs(eta_1) < 2.4 && abs(eta_2) < 2.4')

                for i in tqdm(range(iterationsteps)):

                    h_reco_n = df_reco_f.Histo3D(
                        (
                            f'h_reco_n_it{i}', '',
                            len(eta_bins)-1, array('d', eta_bins),
                            len(phi_bins)-1, array('d', phi_bins),
                            len(mass_bins)-1, array('d', mass_bins)
                        ),
                        'eta_1',
                        'phi_1',
                        mass,
                        "weight"
                    )
                    h_reco_p = df_reco_f.Histo3D(
                        (
                            f'h_reco_p_it{i}', '',
                            len(eta_bins)-1, array('d', eta_bins),
                            len(phi_bins)-1, array('d', phi_bins),
                            len(mass_bins)-1, array('d', mass_bins)
                        ),
                        'eta_2',
                        'phi_2',
                        mass,
                        "weight"
                    )

                    df_reco_f = df_reco_f.Redefine(f"masspt_1", f"{mass} * pt_1_roccor_it")
                    df_reco_f = df_reco_f.Redefine(f"masspt_2", f"{mass} * pt_2_roccor_it")

                    h_reco_mpt_n = df_reco_f.Histo3D(
                        (
                            f'h_reco_mpt_n_it{i}', '',
                            len(eta_bins)-1, array('d', eta_bins),
                            len(phi_bins)-1, array('d', phi_bins),
                            len(masspt_bins)-1, array('d', masspt_bins)
                        ),
                        'eta_1',
                        'phi_1',
                        f'masspt_1',
                        "weight"
                    )
                    h_reco_mpt_p = df_reco_f.Histo3D(
                        (
                            f'h_reco_mpt_p_it{i}', '',
                            len(eta_bins)-1, array('d', eta_bins),
                            len(phi_bins)-1, array('d', phi_bins),
                            len(masspt_bins)-1, array('d', masspt_bins)
                        ),
                        'eta_2',
                        'phi_2',
                        f'masspt_2',
                        "weight"
                    )

                    reco_means_n = h_reco_n.Project3DProfile(option='yx')
                    reco_means_p = h_reco_p.Project3DProfile(option='yx')
                    reco_means_mpt_n = h_reco_mpt_n.Project3DProfile(option='yx')
                    reco_means_mpt_p = h_reco_mpt_p.Project3DProfile(option='yx')

                    for eta in range(len(eta_bins)-1):
                        for phi in range(len(phi_bins)-1):
                            reco_mean_n = reco_means_n.GetBinContent(eta+1, phi+1)
                            reco_mean_p = reco_means_p.GetBinContent(eta+1, phi+1)
                            gen_mean_n = gen_means_n.GetBinContent(eta+1, phi+1)
                            gen_mean_p = gen_means_p.GetBinContent(eta+1, phi+1)

                            reco_means.SetBinContent(
                                i+1, eta+1, phi+1,
                                .5*(reco_mean_n + reco_mean_p)
                                )

                            V_n = -2 * (gen_mean_n - reco_mean_n)
                            V_p = -2 * (gen_mean_p - reco_mean_p)

                            if eta==7 and phi==2: print(V_n, V_p)

                            M_n = 2 * reco_mean_n
                            M_p = 2 * reco_mean_p

                            K_n = - reco_means_mpt_n.GetBinContent(eta+1, phi+1)
                            K_p = reco_means_mpt_p.GetBinContent(eta+1, phi+1)

                            kappa = (V_p / K_p - V_n / K_n) / (M_p / K_p - M_n / K_n) + 1
                            lambd = ((V_p - M_p * (kappa-1)) / K_p + (V_n - M_n * (kappa-1)) / K_n) / 2.

                            kappa_it = h_kappa.GetBinContent(eta+1, phi+1) * kappa
                            lambd_it = h_lambd.GetBinContent(eta+1, phi+1) * kappa + lambd

                            h_kappa.SetBinContent(eta+1, phi+1, kappa_it)
                            h_lambd.SetBinContent(eta+1, phi+1, lambd_it)

                    tf = ROOT.TFile(f'{hdir}step3_it.root', 'recreate')
                    h_kappa.Write()
                    h_lambd.Write()
                    tf.Close()

                    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step3_it.root", "READ");')
                    ROOT.gROOT.ProcessLine(f'TH2D* h_kappa_{i}_{typ} = (TH2D*)tf->Get("M_{typ}");')
                    ROOT.gROOT.ProcessLine(f'h_kappa_{i}_{typ}->SetDirectory(nullptr);')
                    ROOT.gROOT.ProcessLine(f'TH2D* h_lambd_{i}_{typ} = (TH2D*)tf->Get("A_{typ}");')
                    ROOT.gROOT.ProcessLine(f'h_lambd_{i}_{typ}->SetDirectory(nullptr);')
                    ROOT.gROOT.ProcessLine(f'tf->Close();')

                    # application of corrections
                    df_reco = df_reco.Redefine(
                        f"pt_1_roccor_it",
                        f"double pt;\
                        pt = 1./ (h_kappa_{i}_{typ}->GetBinContent( h_kappa_{i}_{typ}->GetXaxis()->FindBin(eta_1) , h_kappa_{i}_{typ}->GetYaxis()->FindBin(phi_1) ) / pt_1 - \
                        h_lambd_{i}_{typ}->GetBinContent( h_lambd_{i}_{typ}->GetXaxis()->FindBin(eta_1) , h_lambd_{i}_{typ}->GetYaxis()->FindBin(phi_1) ));\
                        return pt;"
                    )
                    df_reco = df_reco.Redefine(
                        f"pt_2_roccor_it",
                        f"1./ (h_kappa_{i}_{typ}->GetBinContent( h_kappa_{i}_{typ}->GetXaxis()->FindBin(eta_2) , h_kappa_{i}_{typ}->GetYaxis()->FindBin(phi_2) ) / pt_2 + \
                        h_lambd_{i}_{typ}->GetBinContent( h_lambd_{i}_{typ}->GetXaxis()->FindBin(eta_2) , h_lambd_{i}_{typ}->GetYaxis()->FindBin(phi_2) ))"
                    )
                    df_reco = df_reco.Redefine(
                        f"mass_Z_roccor_it",
                        f"sqrt( 2 * pt_1_roccor_it * pt_2_roccor_it * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)) )"
                    )

                    df_reco_f = df_reco.Filter(f'{mass} > 86 && {mass} < 96')
                    df_reco_f = df_reco.Filter(f'abs(eta_1) < 2.4 && abs(eta_2) < 2.4')

                tf = ROOT.TFile(f'{hdir}step3_it_{typ}.root', 'recreate')
                h_kappa.Write()
                h_lambd.Write()
                tf.Close()
            
            rang = np.linspace(86, 96, 60)
            h_mass_it = df_reco.Histo1D(('h_mass_it', '', len(rang)-1, array('d', rang)), 'mass_Z_roccor_it', "weight")
            h_mass = df_reco.Histo1D(('h_mass', '', len(rang)-1, array('d', rang)), 'mass_Z_roccor', "weight")
            h_mass_gen = df_gen.Histo1D(('h_mass_gen', '', len(rang)-1, array('d', rang)), 'genmass_Z_smeared', "weight")

            tf = ROOT.TFile(f'{hdir}step3_closure_{typ}.root', 'recreate')
            h_mass_it.Write()
            h_mass.Write()
            h_mass_gen.Write()
            reco_means.Write()
            gen_means.Write()
            tf.Close()


def plot_closure(hdir, pdir, samples, eta_bins, phi_bins, iterationsteps):
    pdir = pdir+'iterative/'

    tf_sig = ROOT.TFile(f'{hdir}step3_closure_SIG.root', 'read')
    tf_data = ROOT.TFile(f'{hdir}step3_closure_DATA.root', 'read')

    h_mass_data_it = tf_data.Get('h_mass_it').Clone('h_mass_data_it')
    h_mass_data_it.Sumw2()
    h_mass_data_it.Scale(1./h_mass_data_it.Integral())

    h_mass_sig_it = tf_sig.Get('h_mass_it').Clone('h_mass_sig_it')
    h_mass_sig_it.Sumw2()
    h_mass_sig_it.Scale(1./h_mass_sig_it.Integral())

    h_mass_data = tf_data.Get('h_mass').Clone('h_mass_data')
    h_mass_data.Sumw2()
    h_mass_data.Scale(1./h_mass_data.Integral())

    h_mass_sig = tf_sig.Get('h_mass').Clone('h_mass_sig')
    h_mass_sig.Sumw2()
    h_mass_sig.Scale(1./h_mass_sig.Integral())

    h_mass_gen = tf_data.Get('h_mass_gen')
    h_mass_gen.Sumw2()
    h_mass_gen.Scale(1./h_mass_gen.Integral())

    data_means = tf_data.Get(f'reco_means_DATA')
    sig_means = tf_sig.Get(f'reco_means_SIG')
    gen_means = tf_data.Get(f'gen_means')

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