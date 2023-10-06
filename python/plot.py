import ROOT
from array import array



def plot_ratio(hists, title, outfile, text=['','',''], xrange=[50,130]):
    ROOT.gROOT.SetBatch(1)
    hists['mc'].Scale(hists['dt'].Integral()/hists['mc'].Integral())

    c = ROOT.TCanvas("c", title, 900, 800)
    c.Divide(1,2)
    c.cd(1)
    plotpad = c.GetPad(1)
    plotpad.SetPad(0, 0.21, 1, 1)
    
    hists['mc'].SetStats(0)
    hists['mc'].SetTitle(title)
    hists['mc'].SetLineWidth(2)

    hists['mc'].GetXaxis().SetLabelSize(0)
    hists['mc'].GetXaxis().SetTitleSize(0)
    hists['mc'].GetXaxis().SetRangeUser(xrange[0], xrange[1])

    hists['mc'].GetYaxis().SetRangeUser(0.0001, 1.2* max(hists['dt'].GetMaximum(), hists['mc'].GetMaximum()))
    hists['mc'].GetYaxis().SetTitle('a.u.')
    hists['mc'].GetYaxis().SetTitleSize(0.07)
    hists['mc'].GetYaxis().SetTitleOffset(0.7)
    hists['mc'].Draw("hist")


    hists['dt'].SetMarkerStyle(ROOT.kFullCircle)
    hists['dt'].SetMarkerSize(.8)
    #hists['dt'].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    hists['dt'].Draw("AP0 same")
    ratio = hists['dt'].Clone()


    legend = ROOT.TLegend(0.12, 0.7, 0.3, 0.88)
    legend.AddEntry(hists['mc'], 'DY #rightarrow #mu#mu')
    legend.AddEntry(hists['dt'], 'Muon Data', "lep")
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
    cmsTex.DrawLatex(0.75, 0.915, '(2018, 13 TeV)')
    
    #cmsTex.DrawLatex(0.7, 0.85, 'A-D = {}'.format(test_ad))
    cmsTex.DrawLatex(0.64, 0.86, text[0])
    cmsTex.DrawLatex(0.68, 0.82, text[1])
    cmsTex.DrawLatex(0.68, 0.78, text[2])
    cmsTex.DrawLatex(0.69, 0.65, 'chi2/NDF = {}'.format(round(test_chi2,2)))
    #cmsTex.DrawLatex(0.7, 0.75, 'K-S = {}'.format(test_ks))
    #pad1.SetLogy(1)
    c.cd(2)
    ratiopad = c.GetPad(2)
    ratiopad.SetPad(0, 0, 1, 0.31)
    ratiopad.SetFillStyle(4000)
    ratiopad.SetBottomMargin(.25)

    ratio.SetStats(0)
    ratio.Divide(hists['mc'])
    #hists['dt'].SetMarkerStyle(20)
    ratio.Draw("ep")

    ratio.GetXaxis().SetLabelSize(0.09)
    ratio.GetXaxis().SetTitleOffset(0.7)
    ratio.GetXaxis().SetTitleSize(0.15)
    ratio.GetXaxis().SetTickSize(0.07)
    ratio.GetXaxis().SetTitle('m_#mu#mu (GeV)')

    ratio.GetYaxis().SetRangeUser(0.7,1.3)
    ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetYaxis().SetTitle("Data/MC")
    ratio.GetYaxis().SetTickSize(0.03)
    ratio.GetYaxis().SetTitleOffset(0.3)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.SetTitle("")

    xlim_hi = ratio.GetXaxis().GetXmax()
    xlim_lo = ratio.GetXaxis().GetXmin()
    line = ROOT.TLine(xlim_lo, 1, xlim_hi, 1)
    line.SetLineWidth(2)
    line.Draw("same")

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
            'MC': "h_oneOverPt_MC_neg_pyx",
            'DATA': "h_oneOverPt_DATA_neg_pyx"
        },
        outfile=f"{pdir}oneOverPt_neg",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_pos_pyx",
            'DATA': "h_oneOverPt_DATA_pos_pyx"
        },
        outfile=f"{pdir}oneOverPt_pos",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )

    plot_hists(
        hfile='hists/oneOverPt_roccor.root',
        hists={
            'MC': "h_oneOverPt_MC_neg_roccor_pyx",
            'DATA': "h_oneOverPt_DATA_neg_roccor_pyx"
        },
        outfile=f"{pdir}oneOverPt_neg_roccor",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )
    plot_hists(
        hfile='hists/oneOverPt_roccor.root',
        hists={
            'MC': "h_oneOverPt_MC_pos_roccor_pyx",
            'DATA': "h_oneOverPt_DATA_pos_roccor_pyx"
        },
        outfile=f"{pdir}oneOverPt_pos_roccor",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )
    # median
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_median_neg",
            'DATA': "h_oneOverPt_DATA_median_neg"
        },
        outfile=f"{pdir}oneOverPt_median_neg",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_median_pos",
            'DATA': "h_oneOverPt_DATA_median_pos"
        },
        outfile=f"{pdir}oneOverPt_median_pos",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins
    )
    # peak
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_peak_neg",
            'DATA': "h_oneOverPt_DATA_peak_neg"
        },
        outfile=f"{pdir}oneOverPt_peak_neg",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )
    plot_hists(
        hfile='hists/oneOverPt.root',
        hists={
            'MC': "h_oneOverPt_MC_peak_pos",
            'DATA': "h_oneOverPt_DATA_peak_pos"
        },
        outfile=f"{pdir}oneOverPt_peak_pos",
        dim=2,
        binsx = eta_bins,
        binsy = phi_bins 
    )