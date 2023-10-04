import ROOT
from array import array



def plot_ratio(hists, title, outfile, text=['','',''], xrange=[80,102]):
    hists['mc'].Scale(hists['dt'].Integral()/hists['mc'].Integral())
    c = ROOT.TCanvas("c", title, 800, 700)
    ROOT.gROOT.SetBatch(1)
    pad1 = ROOT.TPad("pad1", "pad1", 0,0.3,1,1)
    pad1.SetBottomMargin(0)
    pad1.Draw()
    pad1.cd()
    hists['dt'].SetMarkerStyle(21)
    hists['dt'].SetMarkerSize(.8)
    hists['mc'].SetStats(0)
    hists['mc'].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    hists['mc'].GetYaxis().SetRangeUser(0.0001, 1.2* hists['dt'].GetMaximum())
    hists['dt'].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    hists['mc'].Draw("same hist")
    hists['dt'].DrawCopy("same")
    hists['mc'].SetLineWidth(2)
    hists['mc'].SetTitle("dilepton mass")

    legend = ROOT.TLegend(0.12, 0.7, 0.3, 0.88)
    legend.AddEntry('mc', 'DY #rightarrow #mu#mu')
    legend.AddEntry('dt', 'Muon Data', "lep")
    legend.SetBorderSize(0)
    legend.Draw('same')

    #print(hists['dt'], hists['mc'])
    #test_ad = hists['dt'].AndersonDarlingTest(hists['mc'], "D")
    test_chi2 = hists['dt'].Chi2Test(hists['mc'], "WW CHI2/NDF")
    #test_ks = hists['dt'].KolmogorovTest(hists['mc'], "WW")
    cmsTex=ROOT.TLatex()
    cmsTex.SetTextFont(42)
    cmsTex.SetTextSize(0.025)
    cmsTex.SetNDC()
    cmsTex.SetTextSize(0.035)
    cmsTex.DrawLatex(0.11,0.915,'#bf{CMS} #it{Preliminary}')
    #cmsTex.DrawLatex(0.745, 0.92, '{} data events'.format(evts))
    cmsTex.DrawLatex(0.67, 0.915, '656.9 pb^{-1} (2022, 13.6 TeV)')
    
    #cmsTex.DrawLatex(0.7, 0.85, 'A-D = {}'.format(test_ad))
    cmsTex.DrawLatex(0.64, 0.86, text[0])
    cmsTex.DrawLatex(0.68, 0.82, text[1])
    cmsTex.DrawLatex(0.68, 0.78, text[2])
    cmsTex.DrawLatex(0.69, 0.65, 'chi2/NDF = {}'.format(round(test_chi2,3)))
    #cmsTex.DrawLatex(0.7, 0.75, 'K-S = {}'.format(test_ks))
    #pad1.SetLogy(1)
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0,0,1,0.3)
    pad2.SetTopMargin(.05)
    pad2.SetBottomMargin(.15)
    pad2.Draw()
    pad2.cd()
    hists['dt'].SetStats(0)
    hists['dt'].Divide(hists['mc'])
    #hists['dt'].SetMarkerStyle(20)
    hists['dt'].Draw("ep")
    hists['dt'].GetYaxis().SetRangeUser(0.8,1.2)
    hists['dt'].GetXaxis().SetLabelSize(0.07)
    hists['dt'].GetYaxis().SetLabelSize(0.07)
    hists['dt'].GetYaxis().SetTitle("dt/mc")
    hists['dt'].GetYaxis().SetTitleOffset(0.4)
    hists['dt'].GetXaxis().SetTitleOffset(1)
    hists['dt'].GetYaxis().SetTitleSize(0.07)
    hists['dt'].GetXaxis().SetTitleSize(0.07)
    hists['dt'].SetTitle("")
    line = ROOT.TLine(xrange[0], 1, xrange[1], 1)
    line.SetLineWidth(2)
    line.Draw("same")
    c.cd()
    c.SaveAs(outfile+'.pdf')    
    c.SaveAs(outfile+".png")
    return test_chi2


def plot_2d_ratio(hists, outfile, binsx, binsy):
    ROOT.gStyle.SetOptStat(0)
    h_ratio = ROOT.TH2D(
        'h_ratio', '', 
        len(binsx)-1, array('f', binsx),
        len(binsy)-1, array('f', binsy)
    )

    for i in range(len(binsx)-1):
        for j in range(len(binsy)-1):
            # print(hists['dt'].GetBinContent(i+1, j+1) / hists['mc'].GetBinContent(i+1, j+1))
            h_ratio.SetBinContent(
                i+1, j+1,
                hists['dt'].GetBinContent(i+1, j+1) / hists['mc'].GetBinContent(i+1, j+1)
            )
    c = ROOT.TCanvas()
    h_ratio.Draw('COLZ')
    c.SaveAs(outfile + '.pdf')
    c.SaveAs(outfile + '.png')


def hist_ntuples(
    ntuples,
    var, nbins, low, up, 
    hdir, fname, option="update"
    ):  
    #make histogram of any variable in ntuple (e.g. "mass_Z" )
    hists=[]
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s])
        h_info=(var+"_"+s, var+" "+s, nbins, low, up)
        h = rdf.Histo1D(h_info, var, 'zPtWeight')
        h.Scale(1./h.Integral())
        hists.append(h)
    
    tf = ROOT.TFile(f"{hdir}{fname}.root", option)
    for h in hists:
        h.Write()
    tf.Close()



def plot_hists(hfile, hists, outfile, binsx=[], binsy=[], dim=1):
    tf = ROOT.TFile(hfile, "READ")
    h_mc = tf.Get(hists['MC'])
    h_dt = tf.Get(hists['DATA'])

    if dim==1:
        plot_ratio(
            hists = {
                'mc': h_mc,
                'dt': h_dt,
            },
            title='',
            outfile=outfile
        )
    elif dim==2:
        plot_2d_ratio(
            hists={
                'mc': h_mc,
                'dt': h_dt
            },
            outfile=outfile,
            binsx=binsx,
            binsy=binsy
        )


def plot_stuff(pdir, eta_bins, phi_bins):
    # Z mass
    plot_hists(
        hfile='hists/mass_z.root',
        hists={
            'MC': "mass_Z_MC",
            'DATA': "mass_Z_DATA"
        },
        outfile=f"{pdir}mass_z"
    )
    plot_hists(
        hfile='hists/mass_z.root',
        hists={
            'MC': "mass_Z_roccor_MC",
            'DATA': "mass_Z_roccor_DATA"
        },
        outfile=f"{pdir}mass_z_roccor"
    )

    # one over pT
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_neg_pxy",
            'DATA': "h_oneOverPt_DATA_neg_pxy"
        },
        outfile=f"{pdir}oneOverPt_neg",
        dim=2,
        binsx = phi_bins,
        binsy = eta_bins 
    )
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_pos_pxy",
            'DATA': "h_oneOverPt_DATA_pos_pxy"
        },
        outfile=f"{pdir}oneOverPt_pos",
        dim=2,
        binsx = phi_bins,
        binsy = eta_bins 
    )

    plot_hists(
        hfile='hists/oneOverPt_roccor.root',
        hists={
            'MC': "h_oneOverPt_MC_neg_roccor_pxy",
            'DATA': "h_oneOverPt_DATA_neg_roccor_pxy"
        },
        outfile=f"{pdir}oneOverPt_neg_roccor",
        dim=2,
        binsx = phi_bins,
        binsy = eta_bins 
    )
    plot_hists(
        hfile='hists/oneOverPt_roccor.root',
        hists={
            'MC': "h_oneOverPt_MC_pos_roccor_pxy",
            'DATA': "h_oneOverPt_DATA_pos_roccor_pxy"
        },
        outfile=f"{pdir}oneOverPt_pos_roccor",
        dim=2,
        binsx = phi_bins,
        binsy = eta_bins 
    )