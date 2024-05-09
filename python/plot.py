import ROOT
from array import array
import os


def plot_ratio(hists, title, outfile, text=['','',''], xrange=None, ratio_range = [0.7, 1.3], labels={'mc':'MC', 'dt': 'Data', 'gen': 'Gen'}, run='2022'):
    ROOT.gROOT.SetBatch(1)

    c = ROOT.TCanvas("c", title, 900, 800)
    c.Divide(1,2)
    c.cd(1)
    c.SetLogy(ROOT.kTRUE)
    plotpad = c.GetPad(1)
    plotpad.SetFillStyle(4000)
    plotpad.SetPad(0, 0.21, 1, 1)
    
    hists['mc'].SetStats(0)
    hists['mc'].SetTitle(title)
    hists['mc'].SetFillColor(ROOT.kGray)

    if xrange:
        hists['mc'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['mc'].GetXaxis().SetLabelSize(0)
    hists['mc'].GetXaxis().SetTitleSize(0)
    hists['mc'].SetLineWidth(0)
    
    hists['mc'].GetYaxis().SetRangeUser(0., 1.2* max(hists['dt'].GetMaximum(), hists['mc'].GetMaximum(), hists['gen'].GetMaximum()))
    hists['mc'].GetYaxis().SetTitle('a.u.')
    hists['mc'].GetYaxis().SetLabelSize(0.04)
    hists['mc'].GetYaxis().SetTitleSize(0.07)
    hists['mc'].GetYaxis().SetTitleOffset(0.7)
    hists['mc'].Draw("hist")

    if xrange:
        hists['dt'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['dt'].SetMarkerStyle(ROOT.kFullCircle)
    hists['dt'].SetMarkerColor(ROOT.kBlack)
    hists['dt'].SetMarkerSize(.8)
    hists['dt'].Draw("ep same")
    ratio_data = hists['dt'].Clone()

    if xrange:
        hists['gen'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['gen'].SetLineWidth(2)
    hists['gen'].SetLineColor(ROOT.kBlue)
    hists['gen'].Draw("el same")
    ratio_gen = hists['gen'].Clone()


    legend = ROOT.TLegend()#0.12, 0.7, 0.3, 0.88)
    legend.AddEntry(hists['mc'], labels['mc'])
    legend.AddEntry(hists['gen'], labels['gen'])
    legend.AddEntry(hists['dt'], labels['dt'])
    legend.SetBorderSize(0)
    legend.Draw('same')

    #print(hists['dt'], hists['mc'])
    #test_ad = hists['dt'].AndersonDarlingTest(hists['mc'], "D")
    test_chi2_dt = hists['dt'].Chi2Test(hists['mc'], "WW CHI2/NDF")
    test_chi2_gen = hists['gen'].Chi2Test(hists['mc'], "WW CHI2/NDF")
    #test_ks = hists['dt'].KolmogorovTest(hists['mc'], "WW")
    cmsTex=ROOT.TLatex()
    cmsTex.SetTextFont(42)
    cmsTex.SetTextSize(0.025)
    cmsTex.SetNDC()
    cmsTex.SetTextSize(0.035)
    cmsTex.DrawLatex(0.11,0.915,'#bf{CMS} #it{Preliminary}')
    #cmsTex.DrawLatex(0.745, 0.92, '{} data events'.format(evts))
    cmsTex.DrawLatex(0.71, 0.915, f'({run}, 13.6 TeV)')
    
    #cmsTex.DrawLatex(0.7, 0.85, 'A-D = {}'.format(test_ad))
    cmsTex.DrawLatex(0.68, 0.86, text[0])
    cmsTex.DrawLatex(0.68, 0.82, text[1])
    cmsTex.DrawLatex(0.68, 0.78, text[2])
    cmsTex.DrawLatex(0.69, 0.65, 'chi2/NDF = {}'.format(round(test_chi2_dt,2)))
    cmsTex.DrawLatex(0.69, 0.62, '#color[{}]{{chi2/NDF = {}}}'.format(ROOT.kBlue, round(test_chi2_gen,2)))
    #cmsTex.DrawLatex(0.7, 0.75, 'K-S = {}'.format(test_ks))
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
    ratio_data.GetXaxis().SetTitle('m_{#mu#mu} (GeV)')

    ratio_data.GetYaxis().SetRangeUser(ratio_range[0], ratio_range[1])
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
    return test_chi2_dt



def plot_ratio2(hists, title, outfile, text=['','',''], xrange=None, ratio_range = [0.7, 1.3], labels={'gen':'GEN', 'mc': 'reco MC', 'gen_smeared': 'smeared Gen'}):
    ROOT.gROOT.SetBatch(1)

    c = ROOT.TCanvas("c", title, 900, 800)
    c.Divide(1,2)
    c.cd(1)
    plotpad = c.GetPad(1)
    plotpad.SetFillStyle(4000)
    plotpad.SetPad(0, 0.21, 1, 1)
    
    hists['mc'].SetStats(0)
    hists['mc'].SetTitle(title)
    hists['mc'].SetFillColor(ROOT.kGray)

    if xrange:
        hists['mc'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['mc'].GetXaxis().SetLabelSize(0)
    hists['mc'].GetXaxis().SetTitleSize(0)
    hists['mc'].SetLineWidth(0)
    
    hists['mc'].GetYaxis().SetRangeUser(0., 1.2* max(hists['gen_smeared'].GetMaximum(), hists['mc'].GetMaximum(), hists['gen'].GetMaximum()))
    hists['mc'].GetYaxis().SetTitle('a.u.')
    hists['mc'].GetYaxis().SetLabelSize(0.04)
    hists['mc'].GetYaxis().SetTitleSize(0.07)
    hists['mc'].GetYaxis().SetTitleOffset(0.7)
    hists['mc'].Draw("hist")

    if xrange:
        hists['gen_smeared'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['gen_smeared'].SetLineWidth(2)
    hists['gen_smeared'].SetLineColor(ROOT.kBlue)
    hists['gen_smeared'].Draw("le same")
    ratio_data = hists['gen_smeared'].Clone()

    if xrange:
        hists['gen'].GetXaxis().SetRange(xrange[0], xrange[1])
    hists['gen'].SetLineWidth(2)
    hists['gen'].SetLineColor(ROOT.kRed)
    hists['gen'].Draw("le same")
    ratio_gen = hists['gen'].Clone()


    legend = ROOT.TLegend()
    legend.AddEntry(hists['mc'], labels['mc'])
    legend.AddEntry(hists['gen'], labels['gen'])
    legend.AddEntry(hists['gen_smeared'], labels['gen_smeared'])
    legend.SetBorderSize(0)
    legend.Draw('same')

    #print(hists['gen_smeared'], hists['mc'])
    #test_ad = hists['gen_smeared'].AndersonDarlingTest(hists['mc'], "D")
    test_chi2_gen = hists['gen'].Chi2Test(hists['mc'], "WW CHI2/NDF")
    test_chi2_gensm = hists['gen_smeared'].Chi2Test(hists['mc'], "WW CHI2/NDF")
    #test_ks = hists['gen_smeared'].KolmogorovTest(hists['mc'], "WW")
    cmsTex=ROOT.TLatex()
    cmsTex.SetTextFont(42)
    cmsTex.SetTextSize(0.025)
    cmsTex.SetNDC()
    cmsTex.SetTextSize(0.035)
    cmsTex.DrawLatex(0.15,0.915,'#bf{CMS} #it{Preliminary}')
    #cmsTex.DrawLatex(0.745, 0.92, '{} data events'.format(evts))
    cmsTex.DrawLatex(0.73, 0.915, '(2022, 13.6 TeV)')
    
    #cmsTex.DrawLatex(0.7, 0.85, 'A-D = {}'.format(test_ad))
    cmsTex.DrawLatex(0.68, 0.86, text[0])
    cmsTex.DrawLatex(0.68, 0.82, text[1])
    cmsTex.DrawLatex(0.68, 0.78, text[2])
    cmsTex.DrawLatex(0.69, 0.65, '#color[{}]{{chi2/NDF = {}}}'.format(ROOT.kRed, round(test_chi2_gen)))
    cmsTex.DrawLatex(0.69, 0.62, '#color[{}]{{chi2/NDF = {}}}'.format(ROOT.kBlue, round(test_chi2_gensm,2)))
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
    #hists['gen_smeared'].SetMarkerStyle(20)
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

    ratio_data.GetYaxis().SetRangeUser(ratio_range[0], ratio_range[1])
    ratio_data.GetYaxis().SetLabelSize(0.09)
    ratio_data.GetYaxis().SetTitle("Gen/Reco")
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
    return test_chi2_gen





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
            rdf = rdf.Define("weight", "zPtWeight*genWeight/sumwWeight*xsec*sf_id*sf_iso")
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


def plot_stuff(pdir, eta_bins, phi_bins, corr):
    # Z mass
    os.makedirs(f'{pdir}mass_z/', exist_ok=True)
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



def roofit_mass(hist_mc, hist_dt, plot=False, fitf='bwxcb', tag='', massrange = [80, 102]):
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

    hist_mc.GetXaxis().SetRangeUser(massrange[0], massrange[1])
    hist_dt.GetXaxis().SetRangeUser(massrange[0], massrange[1])

    results, errors, plots = {}, {}, {}
    if not hist_dt:
        if 'mc' in tag:
            hist_mc.Scale(1000) # scaled to fb
    else:
        hist_mc.Scale(hist_dt.Integral()/hist_mc.Integral())
    #hist_mc.Sumw2(True)
    #hist_dt.Sumw2(True)
    
    x = ROOT.RooRealVar("x", "m_vis (GeV)", massrange[0], massrange[1])
    x.setBins(10000,"cache")
    x.setMin("cache",0)
    x.setMax("cache",500)

    if plot:
        c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
        c1.SetGridx()
        c1.SetGridy()
        c1.GetFrame().SetFillColor( 21 )
        c1.GetFrame().SetBorderMode(-1 )
        c1.GetFrame().SetBorderSize( 5 )
        frame = x.frame()
        frame.SetTitle('')

    #ws = ROOT.RooWorkspace('ws')

    roohist_mc = ROOT.RooDataHist('roomc', 'mc hist', ROOT.RooArgSet(x), hist_mc)
    histpdf = ROOT.RooHistPdf("histpdf", "histpdf", x, roohist_mc, 1) # 1 is order of interpolation
    # definition of gaussian distribution
    mean = ROOT.RooRealVar("mean", "mean", 0, 0, 0) # not really Z mass but for the sake of not having too many exceptions in the code
    sigma = ROOT.RooRealVar("sigma", "sigma", 1.5, 0, 5)
    n_sigma=1
    gaus = ROOT.RooGaussian("gaus", "gaus", x, mean, sigma)

    func = ROOT.RooFFTConvPdf("func", "func", x, histpdf, gaus)

    roohist_dt = ROOT.RooDataHist('dt', 'dt', ROOT.RooArgSet(x), hist_dt)

    #print("="*50)
    fitResult = func.fitTo(roohist_dt)#, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintEvalErrors(-1))
    #ws.Import(func)

    if plot:
        roohist_mc.plotOn(frame, ROOT.RooFit.DrawOption("B"), ROOT.RooFit.FillStyle(0), ROOT.RooFit.FillColor(ROOT.kBlue))
        func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
        #chi2 = frame.chiSquare(6)

        roohist_dt.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack))
        #func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack))
        chi2 = frame.chiSquare(6) 

        # print(d)
        results = [mean.getVal(), n_sigma*sigma.getVal(), chi2]
        errors = [mean.getAsymErrorHi(), n_sigma*sigma.getAsymErrorHi()]
        #errors[d] = [Z_mass.getAsymErrorLo(), sigma.getAsymErrorLo()]

    if plot:
        frame.Draw()
        c1.Update()
        ROOT.gStyle.SetGridColor(ROOT.kGray+1)

        cmsTex=ROOT.TLatex()
        cmsTex.SetTextFont(42)
        cmsTex.SetTextSize(0.025)
        cmsTex.SetNDC()
        #if CMSONLY:
        cmsTex.SetTextSize(0.035)
        cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

        cmsTex.SetTextSize(0.025)
        cmsTex.SetLineWidth(2)
        cmsTex.SetTextFont(42)
        i=0
        stats = ROOT.TPaveText(0.65, 0.65, 0.88, 0.88, "br ARC NDC")
        stats.AddText("M(Z) = {} ({})".format(round(results[0],3), round(errors[0],3)))
        stats.AddText("res(M) = {} ({})".format(round(results[1],3), round(errors[1],3)))
        stats.AddText("chi2/dof = {}".format(round(results[2],2)))

        stats.Draw("SAME")
        c1.Update()
        c1.SaveAs("{}.png".format(plot))
        c1.SaveAs("{}.pdf".format(plot))
    
    return results, errors





def roofit_masscb(hist_mc, hist_dt, plot=False, fitf='bwxcb', tag='', massrange = [80, 102]):
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

    hist_mc.GetXaxis().SetRangeUser(massrange[0], massrange[1])
    hist_dt.GetXaxis().SetRangeUser(massrange[0], massrange[1])

    results, errors, plots = {}, {}, {}
    if not hist_dt:
        if 'mc' in tag:
            hist_mc.Scale(1000) # scaled to fb
    else:
        hist_mc.Scale(hist_dt.Integral()/hist_mc.Integral())
    #hist_mc.Sumw2(True)
    #hist_dt.Sumw2(True)
    
    x = ROOT.RooRealVar("x", "m_vis (GeV)", massrange[0], massrange[1])
    x.setBins(10000,"cache")
    x.setMin("cache",0)
    x.setMax("cache",500)

    if plot:
        c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
        c1.SetGridx()
        c1.SetGridy()
        c1.GetFrame().SetFillColor( 21 )
        c1.GetFrame().SetBorderMode(-1 )
        c1.GetFrame().SetBorderSize( 5 )
        frame = x.frame()
        frame.SetTitle('')
    
    roohist_mc = ROOT.RooDataHist('roomc', 'mc hist', ROOT.RooArgSet(x), hist_mc)

    Z_mass = ROOT.RooRealVar("Z_mass", "Z_mass", 91.1876, 60, 120)
    Z_width = ROOT.RooRealVar("Z_width", "Z_widthan", 2.4952, 0, 10)
    Z_mass.setConstant(True)
    Z_width.setConstant(True)
    mean = ROOT.RooRealVar("mean", "mean", 0, -1, 1)
    #cb_mean.setConstant(True)
    sigma = ROOT.RooRealVar("sigma", "sigma", 1., 0, 3)
    n_sigma = 1 # since sigma_total = sigma_left + sigma_right
    n_L = ROOT.RooRealVar("n_L", "n_L", 5, 0, 1000)
    n_R = ROOT.RooRealVar("n_R", "n_R", 5, 0, 1000)
    alpha_L = ROOT.RooRealVar("alpha_L", "alpha_L", 2, 1.8, 10)
    alpha_R = ROOT.RooRealVar("alpha_R", "alpha_R", 2.5, 2, 10)
    #alpha_L.setConstant(True)
    #alpha_R.setConstant(True)

    bw = ROOT.RooBreitWigner("bw", "BreitWigner", x, Z_mass, Z_width)
    cb = ROOT.RooCrystalBall("cb", "CrystalBall", x, mean, sigma,
                            sigma, alpha_L, n_L, alpha_R, n_R)

    func = ROOT.RooFFTConvPdf("func", "func", x, bw, cb)

    roohist_dt = ROOT.RooDataHist('dt', 'dt', ROOT.RooArgSet(x), hist_dt)

    #print("="*50)
    fitResult = func.fitTo(roohist_dt)#, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintEvalErrors(-1))
    #ws.Import(func)

    if plot:
        roohist_mc.plotOn(frame, ROOT.RooFit.DrawOption("B"), ROOT.RooFit.FillStyle(0), ROOT.RooFit.FillColor(ROOT.kBlue))
        func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
        #chi2 = frame.chiSquare(6)

        roohist_dt.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack))
        #func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack))
        chi2 = frame.chiSquare(6) 

        # print(d)
        results = [mean.getVal(), sigma.getVal(), n_L.getVal(), n_R.getVal(), alpha_L.getVal(), alpha_R.getVal()]
        errors = [mean.getAsymErrorHi(), sigma.getAsymErrorHi()]
        #errors[d] = [Z_mass.getAsymErrorLo(), sigma.getAsymErrorLo()]

    if plot:
        frame.Draw()
        c1.Update()
        ROOT.gStyle.SetGridColor(ROOT.kGray+1)

        cmsTex=ROOT.TLatex()
        cmsTex.SetTextFont(42)
        cmsTex.SetTextSize(0.025)
        cmsTex.SetNDC()
        #if CMSONLY:
        cmsTex.SetTextSize(0.035)
        cmsTex.DrawLatex(0.1,0.92,'#bf{CMS} #it{Preliminary}')

        cmsTex.SetTextSize(0.025)
        cmsTex.SetLineWidth(2)
        cmsTex.SetTextFont(42)
        i=0
        stats = ROOT.TPaveText(0.65, 0.65, 0.88, 0.88, "br ARC NDC")
        stats.AddText("M(Z) = {} ({})".format(round(results[0],3), round(errors[0],3)))
        stats.AddText("res(M) = {} ({})".format(round(results[1],3), round(errors[1],3)))
        stats.AddText("chi2/dof = {}".format(round(results[2],2)))

        stats.Draw("SAME")
        c1.Update()
        c1.SaveAs("{}.png".format(plot))
        c1.SaveAs("{}.pdf".format(plot))
        c1.Delete()
    
    print(results)
    
    return results, errors