import ROOT
from array import array
import os


def plot_ratio(hists, title, outfile, text=['','',''], xrange=None):
    ROOT.gROOT.SetBatch(1)

    c = ROOT.TCanvas("c", title, 900, 800)
    c.Divide(1,2)
    c.cd(1)
    plotpad = c.GetPad(1)
    plotpad.SetFillStyle(4000)
    plotpad.SetPad(0, 0.21, 1, 1)
    
    hists['mc'].SetStats(0)
    hists['mc'].SetTitle(title)
    hists['mc'].SetMarkerStyle(ROOT.kFullCircle)
    hists['mc'].SetMarkerSize(.8)
    hists['mc'].SetMarkerColor(ROOT.kBlack)

    if xrange:
        hists['mc'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['mc'].GetXaxis().SetLabelSize(0)
    hists['mc'].GetXaxis().SetTitleSize(0)
    
    hists['mc'].GetYaxis().SetRangeUser(0., 1.2* max(hists['dt'].GetMaximum(), hists['mc'].GetMaximum(), hists['gen'].GetMaximum()))
    hists['mc'].GetYaxis().SetTitle('a.u.')
    hists['mc'].GetYaxis().SetLabelSize(0.04)
    hists['mc'].GetYaxis().SetTitleSize(0.07)
    hists['mc'].GetYaxis().SetTitleOffset(0.7)
    hists['mc'].Draw("ep")

    if xrange:
        hists['dt'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['dt'].SetMarkerStyle(ROOT.kFullCircle)
    hists['dt'].SetMarkerColor(ROOT.kRed)
    hists['dt'].SetMarkerSize(.8)
    hists['dt'].Draw("ep same")
    ratio_data = hists['dt'].Clone()

    if xrange:
        hists['gen'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['gen'].SetMarkerStyle(ROOT.kFullCircle)
    hists['gen'].SetMarkerColor(ROOT.kGreen)
    hists['gen'].SetMarkerSize(.8)
    hists['gen'].Draw("ep same")
    ratio_gen = hists['gen'].Clone()


    legend = ROOT.TLegend(0.12, 0.7, 0.3, 0.88)
    legend.AddEntry(hists['mc'], 'MC')
    legend.AddEntry(hists['gen'], 'GEN')
    legend.AddEntry(hists['dt'], 'Data')
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
    cmsTex.DrawLatex(0.68, 0.86, text[0])
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

    if xrange:
        ratio_data.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    ratio_data.SetStats(0)
    ratio_data.Divide(hists['mc'])
    #hists['dt'].SetMarkerStyle(20)
    ratio_data.Draw("ep")

    if xrange:
        ratio_gen.GetXaxis().SetRangeUser(xrange[0], xrange[1])
    ratio_gen.Divide(hists['mc'])
    ratio_gen.Draw("ep same")


    ratio_data.GetXaxis().SetLabelSize(0.09)
    ratio_data.GetXaxis().SetTitleOffset(0.7)
    ratio_data.GetXaxis().SetTitleSize(0.15)
    ratio_data.GetXaxis().SetTickSize(0.07)
    ratio_data.GetXaxis().SetTitle('m_#mu#mu (GeV)')

    ratio_data.GetYaxis().SetRangeUser(0.7,1.3)
    ratio_data.GetYaxis().SetLabelSize(0.09)
    ratio_data.GetYaxis().SetTitle("Data/MC")
    ratio_data.GetYaxis().SetTickSize(0.03)
    ratio_data.GetYaxis().SetTitleOffset(0.3)
    ratio_data.GetYaxis().SetTitleSize(0.12)
    ratio_data.SetTitle("")

    xlim_hi = ratio_data.GetXaxis().GetXmax()
    xlim_lo = ratio_data.GetXaxis().GetXmin()
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
    hdir, fname, option="update", corr=''
    ):  
    #make histogram of any variable in ntuple (e.g. "mass_Z" )
    hists={}
    hists["tosave"] = []
    for typ in ntuples:
        hists[typ] = []
        pvar = var+corr
        if typ=='GEN':
            pvar = 'gen'+var
        for sample in ntuples[typ]:
            rdf = ROOT.RDataFrame("Events", ntuples[typ][sample])
            rdf = rdf.Define("weight", "zPtWeight*genWeight*sumwWeight*xsec*sf_id*sf_iso")
            h_info=(pvar+"_"+sample, var+" "+sample, nbins, low, up)
            h = rdf.Histo1D(h_info, pvar, 'weight')
            hists[typ] += [h]

    pvar = var+corr
    # get histograms of Signal MC and Data and GEN
    h_mc = hists['SIG'][0].Clone(pvar+"_SIG")
    h_dt = hists['DATA'][0].Clone(pvar+"_DATA")

    # get histogram of combined backgrounds (not normalized to data)
    h_sum_bkg = hists['BKG'][0].Clone("h_oneOverPt_BKG")
    for i in range(len(hists['BKG'])-1):
        h_tmp = hists['BKG'][i+1].Clone("h_tmp")
        h_sum_bkg.Add(h_tmp)
    
    # calculate combined MC contribution to extract MC->Data normalization SF
    h_mc.Add(h_sum_bkg)
    sf = h_dt.Integral() / h_mc.Integral()

    # scale backgrounds to data
    h_sum_bkg.Scale(sf)
    h_dt.Add(h_sum_bkg, -1)

    hists["tosave"] += [h_dt, hists['SIG'][0], hists['GEN'][0]]        
    
    tf = ROOT.TFile(f"{hdir}{fname}.root", option)
    for h in hists["tosave"]:
        h.Scale(1./h.Integral())
        h.Write()
    tf.Close()



def plot_hists(hfile, hists, outfile, binsx=[], binsy=[], dim=1):
    tf = ROOT.TFile(hfile, "READ")
    print(hists)
    h_mc = tf.Get(hists['SIG'])
    h_dt = tf.Get(hists['DATA'])

    if dim==1:
        h_gen = tf.Get(hists['GEN'])
        plot_ratio(
            hists = {
                'mc': h_mc,
                'dt': h_dt,
                'gen': h_gen
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
    os.makedirs(f'{pdir}mass_z/', exist_ok=True)
    for corr in ['', '_mean_roccor', '_median_roccor']:
        plot_hists(
            hfile='hists/mass_z.root',
            hists={
                'SIG': f"mass_Z{corr}_SIG",
                'DATA': f"mass_Z{corr}_DATA",
                'GEN': "genmass_Z_GEN"
            },
            outfile=f"{pdir}mass_z/mass_z{corr}"
        )

    # one over pT
    hdict = {
        'mean_neg': {
            'SIG': "h_oneOverPt_SIG_neg_pyx",
            'DATA': "h_oneOverPt_DATA_neg_pyx"
        },
        'mean_pos': {
            'SIG': "h_oneOverPt_SIG_pos_pyx",
            'DATA': "h_oneOverPt_DATA_pos_pyx"
        },
        'median_neg': {
            'SIG': "h_oneOverPt_SIG_neg_median",
            'DATA': "h_oneOverPt_DATA_neg_median"
        },
        'median_pos': {
            'SIG': "h_oneOverPt_SIG_pos_median",
            'DATA': "h_oneOverPt_DATA_pos_median"
        }
    }
    os.makedirs(f'{pdir}oneOverPt/', exist_ok=True)
    for h in hdict:
        plot_hists(
            hfile='hists/oneOverPt.root',
            hists=hdict[h],
            outfile=f"{pdir}oneOverPt/oneOverPt_{h}",
            dim=2,
            binsx = eta_bins,
            binsy = phi_bins 
        )
        plot_hists(
            hfile='hists/oneOverPt_mean_roccor.root',
            hists=hdict[h],
            outfile=f"{pdir}oneOverPt/oneOverPt_{h}_mean_roccor",
            dim=2,
            binsx = eta_bins,
            binsy = phi_bins 
        )

    # specific bin
    for corr in ['', '_mean_roccor', '_median_roccor']:
        tf = ROOT.TFile(f"hists/oneOverPt{corr}.root")
        h_mc = tf.Get("h_oneOverPt_SIG_neg").ProjectionZ("h_mc", 3,3,4,4)
        h_mc.Scale(1./h_mc.Integral())
        h_dt = tf.Get("h_oneOverPt_DATA_neg").ProjectionZ("h_dt", 3,3,4,4)
        h_dt.Scale(1./h_dt.Integral())
        h_gen = tf.Get("h_oneOverPt_GEN_neg").ProjectionZ("h_gen", 3,3,4,4)
        h_gen.Scale(1./h_gen.Integral())
        print(h_mc.GetBinContent(2), h_dt.GetBinContent(2), h_gen.GetBinContent(2))
        plot_ratio(
            hists = {
                'mc': h_mc,
                'dt': h_dt,
                'gen': h_gen
            },
            title='',
            outfile=f"{pdir}oneOverPt/bin_3_4{corr}",
            text=[
                f'MC: {round(h_mc.GetMean(),7)}',
                f'DATA: {round(h_dt.GetMean(),7)}',
                f'GEN: {round(h_gen.GetMean(),7)}'
            ],
            xrange=[1, h_mc.GetXaxis().GetNbins()]
        )
        h_yx = tf.Get("h_oneOverPt_SIG_neg_pyx")
        h_yx.SetBinContent(200, 0)
        print(h_yx.GetBinContent(3,4))
