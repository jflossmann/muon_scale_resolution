import ROOT
from array import array
from tqdm import tqdm
import numpy as np
import os
import glob

# TODO: check correlation matrices somehow


# years = ['2022C', '2022D', '2022E', '2022F', '2022G']
years = ['2022_nlo', '2022EE_nlo']

systs = ['stat', 'bin_step3', 'bin_step4']

def get_stuff(hname, eta_bins, phi_bins, var_bins, files):
    hists = {}

    if "step3" in files[0]:
        hname_new = hname+"3"
    else:
        hname_new = hname

    for am in ["M_", "A_"]:


        hist_all = ROOT.TH3D(
            am + hname_new+'_all', '',
            len(eta_bins)-1, array('d', eta_bins),
            len(phi_bins)-1, array('d', phi_bins),
            len(var_bins[am])-1, array('d', var_bins[am])
        )

        for f in files:

            if not os.path.exists(f):
                print(f) 
                continue

            tf = ROOT.TFile(f, 'read')
            hist = tf.Get(am + hname)
            hist.SetDirectory(ROOT.nullptr)
            tf.Close()

            for eta in range(len(eta_bins)-1):
                for phi in range(len(phi_bins)-1):

                    e = (eta_bins[eta] + eta_bins[eta+1])/2.
                    p = (phi_bins[phi] + phi_bins[phi+1])/2.

                    hist_all.Fill(
                        e, p,
                        hist.GetBinContent(eta+1, phi+1)
                    )

        h_tmp = hist_all.Project3DProfile("yx")
        h_tmp.Sumw2()
        h_tmp.SetErrorOption("s")

        hists[am] = h_tmp

    return hists


# get correlations
def get_correlations_scale(hname, eta_bins, phi_bins, var_bins, var1, var2, files1, files2, pdir, hdir):
    if "step3" in files1[0]:
        hname_new = hname+"3"
    else:
        hname_new = hname

    hist_correlations = ROOT.TH2D(
        'correlation'+hname_new, '',
        len(eta_bins)-1, array('d', eta_bins),
        len(phi_bins)-1, array('d', phi_bins)
    )
    for eta in range(len(eta_bins)-1):
        for phi in range(len(phi_bins)-1):

            hist_corr_tmp = ROOT.TH2D(
                f'hist_corr_tmp{eta}{phi}', '',
                len(var_bins[var1])-1, array('d', var_bins[var1]),
                len(var_bins[var2])-1, array('d', var_bins[var2])
            )

            if not len(files1)==len(files2):
                print("Number of inputs files must be of same length")
                continue

            for fno in range(len(files1)):
                tf_corr1 = ROOT.TFile(files1[fno], 'read')
                tf_corr2 = ROOT.TFile(files2[fno], 'read')

                hist1 = tf_corr1.Get(var1+hname)
                hist2 = tf_corr2.Get(var2+hname)

                if not hist1 or not hist2:
                    print(f"Histgrams {var1+hname} or {var2+hname} not found in file: {f}")
                    continue

                hist_corr_tmp.Fill(
                    hist1.GetBinContent(eta+1, phi+1),
                    hist2.GetBinContent(eta+1, phi+1)
                )

                tf_corr1.Close()
                tf_corr2.Close()

            hist_correlations.SetBinContent(
                eta+1, phi+1,
                hist_corr_tmp.GetCorrelationFactor()
            )

            del hist_corr_tmp
    
    tf = ROOT.TFile(f'{hdir}{var1}{var2}{hname}.root', 'recreate')
    hist_correlations.Write()
    tf.Close()

    c = ROOT.TCanvas()
    hist_correlations.Draw("colz")
    hist_correlations.GetXaxis().SetTitle("#eta")
    hist_correlations.GetYaxis().SetTitle("#phi")
    hist_correlations.SetStats(0)

    cmsTex=ROOT.TLatex()
    cmsTex.SetTextFont(42)
    cmsTex.SetTextSize(0.025)
    cmsTex.SetNDC()
    #if CMSONLY:
    cmsTex.SetTextSize(0.035)
    cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

    c.SaveAs(f'{pdir}{var1}{var2}{hname}.pdf')

    return


# get correlations
def get_correlations_resol(abseta_bins, var_bins, var, files, pdir, hdir):
    hist_correlations = ROOT.TH1D(
        'correlation_k_DATA_SIG', '',
        len(abseta_bins)-1, array('d', abseta_bins),
    )
    for abseta in range(len(abseta_bins)-1):

        hist_corr_tmp = ROOT.TH2D(
            f'hist_corr_tmp{abseta}', '',
            len(var_bins[var])-1, array('d', var_bins[var]),
            len(var_bins[var])-1, array('d', var_bins[var])
        )

        for fno in range(len(files)):
            tf_corr1 = ROOT.TFile(files[fno], 'read')
            tf_corr2 = ROOT.TFile(files[fno], 'read')

            hist1 = tf_corr1.Get(var+'DATA')
            hist2 = tf_corr2.Get(var+'SIG')

            if not hist1 or not hist2:
                print(f"Histgrams {var}DATA or {var}SIG not found in file: {files[fno]}")
                continue

            hist_corr_tmp.Fill(
                hist1.GetBinContent(abseta+1, 3),
                hist2.GetBinContent(abseta+1, 3)
            )

            tf_corr1.Close()
            tf_corr2.Close()

        hist_correlations.SetBinContent(
            abseta+1,
            hist_corr_tmp.GetCorrelationFactor()
        )

        del hist_corr_tmp
    
    tf = ROOT.TFile(f'{hdir}{var}DATA_SIG.root', 'recreate')
    hist_correlations.Write()
    tf.Close()

    c = ROOT.TCanvas()
    hist_correlations.Draw("colz")
    hist_correlations.GetXaxis().SetTitle("#eta")
    hist_correlations.GetYaxis().SetTitle("#phi")
    hist_correlations.SetStats(0)

    cmsTex=ROOT.TLatex()
    cmsTex.SetTextFont(42)
    cmsTex.SetTextSize(0.025)
    cmsTex.SetNDC()
    #if CMSONLY:
    cmsTex.SetTextSize(0.035)
    cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

    c.SaveAs(f'{pdir}{var}DATA_SIG.pdf')

    return


for year in years:
    for syst in systs:
        files_1 = glob.glob(f'../hists/{year}/condor/{syst}/*/step1_C.root')
        files_3 = glob.glob(f'../hists/{year}/condor/{syst}/*/step3_correction.root')
        files_4 = glob.glob(f'../hists/{year}/condor/{syst}/*/step4_k.root')

        pdir = f'../plots/{year}/systs/{syst}/'
        os.makedirs(pdir, exist_ok=True)
        hdir = f'../hists/{year}/systs/'
        os.makedirs(hdir, exist_ok=True)

        ROOT.gROOT.SetBatch(1)

        eta_bins = [-2.4, -2.1, -1.85, -0.4, 0, 0.4, 1.85, 2.1, 2.4]
        phi_bins = np.linspace(-3.2, 3.2, 17)
        var_bins = {
            "M_": np.linspace(0.995, 1.005, 1000),
            "A_": np.linspace(-5e-4, 5e-4, 1000),
            "k_hist_": np.linspace(.9, 1.4, 1000)
        }
        n_bs_bins = 18 # 400
        abseta_bins = np.linspace(0, 2.4, 13)


        # make correlation plots
        # get_correlations_scale('DATA', eta_bins, phi_bins, var_bins, 'M_', 'A_', files_3, files_3, pdir, hdir)
        # get_correlations_scale('SIG', eta_bins, phi_bins, var_bins, 'M_', 'A_', files_3, files_3, pdir, hdir)
        get_correlations_resol(abseta_bins, var_bins, 'k_hist_', files_4, pdir, hdir)
        # get_correlations('DATA', eta_bins, phi_bins, var_bins, 'M_', 'k_', files_3, files_4, pdir)
        # get_correlations('SIG', eta_bins, phi_bins, var_bins, 'M_', 'k_', files_3, files_4, pdir)
        # get_correlations('DATA', eta_bins, phi_bins, var_bins, 'A_', 'k_', files_3, files_4, pdir)
        # get_correlations('SIG', eta_bins, phi_bins, var_bins, 'A_', 'k_', files_3, files_4, pdir)

        c = ROOT.TCanvas()
        leg = ROOT.TLegend()
        colors = {
            "DATA": ROOT.kRed,
            "SIG": ROOT.kBlue
        }
        labels = {
            "DATA": "Data",
            "SIG": "MC"
        }
        # get resolution correction factors k
        tfname = f'{hdir}{syst}.root'
        tf_systs = ROOT.TFile(tfname, 'recreate')
        tf_systs.Close()
        for dtsg in ["SIG", "DATA"]:
            hist_all = ROOT.TH2D(
                f'k_all_{dtsg}', '',
                len(abseta_bins)-1, array('d', abseta_bins),
                len(var_bins["k_hist_"])-1, array('d', var_bins["k_hist_"])
            )
            for f in files_4:
                tf = ROOT.TFile(f, 'read')
                hist = tf.Get("k_hist_"+dtsg)
                hist.SetDirectory(ROOT.nullptr)
                tf.Close()
                for eta in range(len(abseta_bins)-1):
                    e = (abseta_bins[eta] + abseta_bins[eta+1])/2.

                    hist_all.Fill(e, hist.GetBinContent(eta+1, 3))

            print(hist_all)  

            h_tmp = hist_all.ProfileX()
            h_tmp.Sumw2()
            h_tmp.SetErrorOption("s")

            h_tmp.SetStats(0)

            h_tmp.GetYaxis().SetRangeUser(.95, 1.15)

            h_tmp.SetLineColor(colors[dtsg])
            h_tmp.SetMarkerStyle(8)
            h_tmp.SetMarkerColor(colors[dtsg])
            h_tmp.SetTitle('')
            h_tmp.GetXaxis().SetTitle('#eta')
            h_tmp.GetYaxis().SetTitle("k")
            h_tmp.Draw("e same")

            leg.AddEntry(h_tmp, dtsg)

            tf_systs = ROOT.TFile(tfname, 'update')
            h_tmp.Write()
            tf_systs.Close()

        leg.Draw("same")
        cmsTex=ROOT.TLatex()
        cmsTex.SetTextFont(42)
        cmsTex.SetTextSize(0.025)
        cmsTex.SetNDC()
        #if CMSONLY:
        cmsTex.SetTextSize(0.035)
        cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

        c.SaveAs(f'{pdir}k.pdf')

        hist_dict = {}
        for dtsig in ["DATA", "SIG"]:
            hist_dict[f"{dtsig}_step1"] = get_stuff(
                hname=dtsig,
                eta_bins=eta_bins,
                phi_bins=phi_bins,
                var_bins=var_bins,
                files=files_1
            )
            hist_dict[f"{dtsig}_step3"] = get_stuff(
                hname=dtsig,
                eta_bins=eta_bins,
                phi_bins=phi_bins,
                var_bins=var_bins,
                files=files_3
            )



        # colors = [ROOT.kMagenta, ROOT.kRed, ROOT.kGreen, ROOT.kBlue]
        colors = [ROOT.kRed, ROOT.kBlue]
        # labels = ['data 1/pt corr', 'data it. corr', 'mc 1/pt corr', 'mc it. corr']
        labels = ['data it. corr', 'mc it. corr']
        markers = [8, 8]
        title = {
            "M_": "#kappa",
            "A_": "#lambda"
        }

        for eta in range(1, len(eta_bins)):

            for am in ["M_", "A_"]:
                c1 = ROOT.TCanvas()
                leg = ROOT.TLegend()
                i = 0

                for dtsig in ["DATA", "SIG"]:
                    # for step in ["step1", f"step3"]:
                    for step in ["step3"]:

                        if eta == 1:
                            tf_systs = ROOT.TFile(tfname, 'update')
                            hist_dict[f"{dtsig}_{step}"][am].Write()
                            tf_systs.Close()

                        h_tmp = hist_dict[f"{dtsig}_{step}"][am].ProjectionY(f'h_{dtsig}_{step}', eta, eta)

                        h_tmp.SetStats(0)

                        if "M" in am:
                            h_tmp.GetYaxis().SetRangeUser(.995, 1.005)
                        else:
                            h_tmp.GetYaxis().SetRangeUser(-2.5e-4, 2.5e-4)

                        h_tmp.SetLineColor(colors[i])
                        h_tmp.SetMarkerStyle(markers[i])
                        h_tmp.SetMarkerColor(colors[i])
                        h_tmp.SetTitle('')
                        h_tmp.GetXaxis().SetTitle('#phi')
                        h_tmp.GetYaxis().SetTitle(title[am])
                        h_tmp.Draw("e same")

                        leg.AddEntry(h_tmp, labels[i])

                        i += 1

                leg.Draw("same")
                c1.SaveAs(f'{pdir}{am}eta_{eta}.pdf')
