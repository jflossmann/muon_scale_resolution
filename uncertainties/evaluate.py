import ROOT
from array import array
from tqdm import tqdm
import numpy as np
import os

pdir = 'plots/scale/'
os.makedirs(pdir, exist_ok=True)
hdir = 'hists/scale/'
os.makedirs(hdir, exist_ok=True)

ROOT.gROOT.SetBatch(1)

eta_bins = [-2.4, -2.1, -1.85, -0.4, 0, 0.4, 1.85, 2.1, 2.4]
phi_bins = np.linspace(-3.2, 3.2, 17)
var_bins = {
    "M_": np.linspace(0.995, 1.005, 1000),
    "A_": np.linspace(-5e-4, 5e-4, 10000),
    "k_": np.linspace(.9, 1.4, 1000)
}
n_bs_bins = 200
abseta_bins = np.linspace(0, 2.4, 13)

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

for dtsg in ["SIG", "DATA"]:
    hist_all = ROOT.TH2D(
        f'k_all_{dtsg}', '',
        len(abseta_bins)-1, array('d', abseta_bins),
        len(var_bins["k_"])-1, array('d', var_bins["k_"])
    )
    for i in range(n_bs_bins):
        tf = ROOT.TFile(f'../hists/condor/{i}/step4_k.root', 'read')
        hist = tf.Get("h_k_"+dtsg)
        hist.SetDirectory(ROOT.nullptr)
        tf.Close()
        for eta in range(len(abseta_bins)-1):
            e = (abseta_bins[eta] + abseta_bins[eta+1])/2.

            hist_all.Fill(e, hist.GetBinContent(eta+1))        

    h_tmp = hist_all.ProfileX()
    h_tmp.Sumw2()
    h_tmp.SetErrorOption("s")

    h_tmp.SetStats(0)

    h_tmp.GetYaxis().SetRangeUser(.95, 1.3)

    h_tmp.SetLineColor(colors[dtsg])
    h_tmp.SetMarkerStyle(8)
    h_tmp.SetMarkerColor(colors[dtsg])
    h_tmp.SetTitle('')
    h_tmp.GetXaxis().SetTitle('#eta')
    h_tmp.GetYaxis().SetTitle("k")
    h_tmp.Draw("e same")

    leg.AddEntry(h_tmp, dtsg)

    i += 1

leg.Draw("same")
cmsTex=ROOT.TLatex()
cmsTex.SetTextFont(42)
cmsTex.SetTextSize(0.025)
cmsTex.SetNDC()
#if CMSONLY:
cmsTex.SetTextSize(0.035)
cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

c.SaveAs(f'{pdir}k.pdf')

c = ROOT.TCanvas()
hist_tmp = hist_all.ProjectionY("hist_tmp", 5, 5)
hist_tmp.Draw()
c.SaveAs('test.pdf')




def get_stuff(nfiles, fname, hname, eta_bins, phi_bins, var_bins):
    hists = {}

    for am in ["M_", "A_"]:
        if "step3" in fname:
            hname_new = am + hname+"3"
        else:
            hname_new = am + hname

        hist_all = ROOT.TH3D(
            hname_new+'_all', '',
            len(eta_bins)-1, array('d', eta_bins),
            len(phi_bins)-1, array('d', phi_bins),
            len(var_bins[am])-1, array('d', var_bins[am])
        )

        for i in range(n_bs_bins):
            path = f'../hists/condor/{i}/{fname}'
            if not os.path.exists(path): continue

            tf = ROOT.TFile(path, 'read')
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


hist_dict = {}
for dtsig in ["DATA", "SIG"]:
    for step in ["step1_C", f"step3_it_{dtsig}"]:
        hist_dict[f"{dtsig}_{step}"] = get_stuff(
            nfiles=200,
            fname=f'{step}.root',
            hname=dtsig,
            eta_bins=eta_bins,
            phi_bins=phi_bins,
            var_bins=var_bins
        )



colors = [ROOT.kMagenta, ROOT.kRed, ROOT.kGreen, ROOT.kBlue]
labels = ['data 1/pt corr', 'data it. corr', 'mc 1/pt corr', 'mc it. corr']
markers = [4, 8, 4, 8]
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
            for step in ["step1_C", f"step3_it_{dtsig}"]:

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
