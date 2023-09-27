import ROOT
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.optimize as opt

def get_canvas(num):
    height = int(np.sqrt(num))
    if num/height == num//height: wide = num//height    
    else: wide = num//height +1
    return height, wide

def plot_canvas(plots, title, filename, style, line=None):
    height, wide = get_canvas(len(plots))
    c = ROOT.TCanvas("c", title, wide*800, height*600)
    c.Divide(wide, height)
    i=1
    for plot in plots:
        c.cd(i)
        legend = ROOT.TLegend(0.6,0.7,0.9,0.9)
        hist = plots[plot]
        hist.SetStats(0)
        hist.Draw(style)
        legend.AddEntry(plot+""+filename.split(".")[0], plot, "l")
        if line: line.Draw("SAME")
        # legend.Draw("SAME")
        i+=1
    c.SaveAs(filename)


def plot_ratio(plots, title, outfile, text=['','',''], xrange=[80,102]):
    plots['mc'].Scale(plots['dt'].Integral()/plots['mc'].Integral())
    c = ROOT.TCanvas("c", title, 800, 700)
    ROOT.gROOT.SetBatch(1)
    pad1 = ROOT.TPad("pad1", "pad1", 0,0.3,1,1)
    pad1.SetBottomMargin(0)
    pad1.Draw()
    pad1.cd()
    plots['dt'].SetMarkerStyle(21)
    plots['dt'].SetMarkerSize(.8)
    plots['mc'].SetStats(0)
    plots['mc'].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    plots['mc'].GetYaxis().SetRangeUser(0.0001, 1.2* plots['dt'].GetMaximum())
    plots['dt'].GetXaxis().SetRangeUser(xrange[0], xrange[1])
    plots['mc'].Draw("same hist")
    plots['dt'].DrawCopy("same")
    plots['mc'].SetLineWidth(2)
    plots['mc'].SetTitle("dilepton mass")

    legend = ROOT.TLegend(0.12, 0.7, 0.3, 0.88)
    legend.AddEntry('mc', 'DY #rightarrow #mu#mu')
    legend.AddEntry('dt', 'Muon Data', "lep")
    legend.SetBorderSize(0)
    legend.Draw('same')

    print(plots['dt'], plots['mc'])
    #test_ad = plots['dt'].AndersonDarlingTest(plots['mc'], "D")
    test_chi2 = plots['dt'].Chi2Test(plots['mc'], "WW CHI2/NDF")
    #test_ks = plots['dt'].KolmogorovTest(plots['mc'], "WW")
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
    plots['dt'].SetStats(0)
    plots['dt'].Divide(plots['mc'])
    #plots['dt'].SetMarkerStyle(20)
    plots['dt'].Draw("ep")
    plots['dt'].GetYaxis().SetRangeUser(0.8,1.2)
    plots['dt'].GetXaxis().SetLabelSize(0.07)
    plots['dt'].GetYaxis().SetLabelSize(0.07)
    plots['dt'].GetYaxis().SetTitle("dt/mc")
    plots['dt'].GetYaxis().SetTitleOffset(0.4)
    plots['dt'].GetXaxis().SetTitleOffset(1)
    plots['dt'].GetYaxis().SetTitleSize(0.07)
    plots['dt'].GetXaxis().SetTitleSize(0.07)
    plots['dt'].SetTitle("")
    line = ROOT.TLine(xrange[0], 1, xrange[1], 1)
    line.SetLineWidth(2)
    line.Draw("same")
    c.cd()
    c.SaveAs(outfile+'.pdf')    
    c.SaveAs(outfile+".png")
    return test_chi2


def plot_same(plots, rcolors, title, outfile, text={}, ratio=False):
    c = ROOT.TCanvas("c", title, 800, 600)
    legend = ROOT.TLegend(0.7,0.75,0.9,0.9)
    ROOT.gROOT.SetBatch(1)

    if ratio:

        rp = ROOT.TRatioPlot(plots['dt'], plots['mc'])

        rp.GetXaxis().SetRangeUser(80,102)
        rp.Draw()
        rp.GetLowerRefGraph().SetMinimum(0.8)
        rp.GetLowerRefGraph().SetMaximum(1.2)
        rp.GetLowerRefYaxis().SetTitle("Ratio")
        
        plots['mc'].SetTitle('mc')
        plots['dt'].SetTitle('data')

        cmsTex=ROOT.TLatex()
        cmsTex.SetTextFont(42)
        cmsTex.SetTextSize(0.025)
        cmsTex.SetNDC()
        #if CMSONLY:
        cmsTex.SetTextSize(0.035)
        cmsTex.DrawLatex(0.11,0.91,'#bf{CMS} #it{Preliminary}')
        rp.GetUpperPad().BuildLegend()
        rp.GetUpperPad().cd()
        rp.SetH1DrawOpt("PE")
        plots['dt'].SetMarkerStyle(20)
        plots['dt'].Draw("PE SAME")
        rp.GetUpperPad().Update()


    else:
        for plot in plots:
            hist = plots[plot]
            hist.SetLineColor(rcolors[plot])
            hist.SetStats(0)
            hist.Draw("SAME")
            legend.AddEntry(plot, plot, "l")

        legend.Draw("SAME")
        cmsTex=ROOT.TLatex()
        cmsTex.SetTextFont(42)
        cmsTex.SetTextSize(0.025)
        cmsTex.SetNDC()
        #if CMSONLY:
        cmsTex.SetTextSize(0.035)
        cmsTex.DrawLatex(0.11,0.91,'#bf{CMS} #it{Preliminary}')


    cmsTex.SetTextSize(0.025)
    cmsTex.SetLineWidth(2)
    cmsTex.SetTextFont(42)
    i = 0
    for t in text:
        cmsTex.DrawLatex(0.75,0.7 - i/10., text[t][0])
        cmsTex.DrawLatex(0.75,0.65 - i/10., text[t][1])

        i+=1
    c.Update
    c.SaveAs(outfile)    
    c.SaveAs(outfile.split(".pdf")[0]+".png")


def plot2d(matrix, outfile, title, x, xbins, y, ybins, cmin, cmax, xticks, yticks):
    nbinsx, nbinsy = len(xbins)-1, len(ybins)-1

    fig, ax = plt.subplots()
    #plt.margins(x=0, y=0)
    # plt.imshow(np.transpose(matrix))
    plt.imshow(matrix)
    plt.colorbar()
    plt.clim(cmin, cmax)
    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    ax.set_xticks(np.arange(len(xticks)))
    ax.set_xticklabels(np.round(xticks,2))    
    ax.set_yticks(np.arange(len(yticks)))
    ax.set_yticklabels(np.round(np.flip(yticks),2))

    for i in range(nbinsx):
        for j in range(nbinsy):
            test = ax.text(i, j, round(matrix[j][i],3), ha='center', va='center', color='w')
    print(outfile)

    plt.savefig(outfile+'.pdf')
    plt.savefig(outfile+'.png')
    plt.clf()


def roofit_mass(hist_mc, hist_dt, plot=False, fitf='bwxcb', pdf=False, e=0, p=0):
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

    results, errors, plots = {}, {}, {}
    if not hist_dt:
        hist_mc.Scale(1000) # scaled to fb
    else:
        hist_mc.Scale(hist_dt.Integral()/hist_mc.Integral())
    #hist_mc.Sumw2(True)
    #hist_dt.Sumw2(True)
    
    hists = {'mc': hist_mc, 'dt': hist_dt}

    x = ROOT.RooRealVar("x", "m_vis (GeV)", 85, 97)
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

    ws = ROOT.RooWorkspace('ws')

    for d in ['dt', 'mc']:

        if not hists[d]:
            continue

        if fitf=='vgt':
            Z_mass = ROOT.RooRealVar("Z_mass", "Z_mass", 91.1876, 60, 120)
            Z_width = ROOT.RooRealVar("Z_width", "Z_widthan", 2.4952, 2.4952, 2.4952)
            Z_width.setConstant(True)
            sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 0, 10)
            n_sigma=1

            func = ROOT.RooVoigtian("vgt_mc", "Voigt", x, Z_mass, Z_width, sigma)

        elif fitf=='bwxcb':
            Z_mass = ROOT.RooRealVar("Z_mass", "Z_mass", 91.1876, 60, 120)
            Z_width = ROOT.RooRealVar("Z_width", "Z_widthan", 2.4952, 0, 10)
            Z_mass.setConstant(True)
            Z_width.setConstant(True)
            mean = ROOT.RooRealVar("mean", "mean", 0, -10, 10)
            #cb_mean.setConstant(True)
            sigma = ROOT.RooRealVar("sigma", "sigma", 1.42, 0, 10)
            n_sigma = 1 # since sigma_total = sigma_left + sigma_right
            n_L = ROOT.RooRealVar("n_L", "n_L", 5, 0, 1000)
            n_R = ROOT.RooRealVar("n_R", "n_R", 5, 0, 1000)
            alpha_L = ROOT.RooRealVar("alpha_L", "alpha_L", 1.3, 1, 5)
            alpha_R = ROOT.RooRealVar("alpha_R", "alpha_R", 1.2, 1, 5)
            #alpha_L.setConstant(True)
            #alpha_R.setConstant(True)

            bw = ROOT.RooBreitWigner("bw", "BreitWigner", x, Z_mass, Z_width)
            cb = ROOT.RooCrystalBall("cb", "CrystalBall", x, mean, sigma,
                                    sigma, alpha_L, n_L, alpha_R, n_R)

            func = ROOT.RooFFTConvPdf("func"+d, "func"+d, x, bw, cb)

        elif fitf=='template':
            roohist_mc = ROOT.RooDataHist('mc', 'mc hist', ROOT.RooArgSet(x), hist_mc)

            histpdf = ROOT.RooHistPdf("histpdf", "histpdf", x, roohist_mc, 1) # 1 is order of interpolation
            # definition of gaussian distribution
            mean = ROOT.RooRealVar("mean", "mean", 0, -100, 100) # not really Z mass but for the sake of not having too many exceptions in the code
            sigma = ROOT.RooRealVar("sigma", "sigma", 1.5, 0, 5)
            n_sigma=1
            gaus = ROOT.RooGaussian("gaus", "gaus", x, mean, sigma)

            func = ROOT.RooFFTConvPdf("func", "func", x, histpdf, gaus)

        else:
            # definition of gaussian distribution
            mean_g = ROOT.RooRealVar("mean_g", "mean_g", 0, 0, 10)
            sigma = ROOT.RooRealVar("sigma_g", "sigma_g", 1, 0, 5)
            n_sigma=1
            gaus = ROOT.RooGaussian("gaus", "gaus", x, mean_g, sigma)

            # definition of dscb
            Z_mass = ROOT.RooRealVar("mean_cb", "mean_cb", 90, 0, 100)
            sigma_cb = ROOT.RooRealVar("sigma_cb", "sigma_cb", 1, 0, 10)
            alpha_L = ROOT.RooRealVar("alpha_L", "alpha_L", 1, 0, 100)
            alpha_R = ROOT.RooRealVar("alpha_R", "alpha_R", 1, 0, 100)
            n_L = ROOT.RooRealVar("n_L", "n_L", 5, 0, 100)
            n_R = ROOT.RooRealVar("n_R", "n_R", 5, 0, 100)
            cb = ROOT.RooCrystalBall("cb", "CrystalBall", x, Z_mass, sigma_cb, sigma_cb, alpha_L, n_L, alpha_R, n_R)

            mean_g.setVal(0)
            mean_g.setConstant(True)
            sigma_cb.setVal(2.4952/2.)
            sigma_cb.setConstant(True)

            func = ROOT.RooFFTConvPdf("func", "func", x, cb, gaus)
        

        roohist = ROOT.RooDataHist(d, d, ROOT.RooArgSet(x), hists[d])
        fitResult = func.fitTo(roohist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintEvalErrors(-1))
        ws.Import(func)

        if plot:
            if d=='mc':
                roohist.plotOn(frame, ROOT.RooFit.DrawOption("B"), ROOT.RooFit.FillStyle(0), ROOT.RooFit.FillColor(ROOT.kBlue))
                func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue))
                chi2 = frame.chiSquare(6)

            if d=='dt':
                roohist.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack))
                func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack))
                chi2 = frame.chiSquare(6) 

        # print(d)
        results[d] = [mean.getVal(), n_sigma*sigma.getVal(), chi2]
        errors[d] = [mean.getAsymErrorHi(), n_sigma*sigma.getAsymErrorHi()]
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
        if hist_dt:
            stats.AddText("M(Z) = {} ({})".format(round(results['dt'][0],3), round(errors[d][0],3)))
            stats.AddText("res(M) = {} ({})".format(round(results['dt'][1],3), round(errors[d][1],3)))
            stats.AddText("chi2/dof = {}".format(round(results['dt'][2],2)))
            stats.AddLine(.0,.5,1.,.5)
        if hist_mc:
            stats.AddText("M(Z) = {} ({})".format(round(results['mc'][0],3), round(errors[d][0],3)))
            stats.GetListOfLines().Last().SetTextColor(ROOT.kBlue)
            stats.AddText("res(M) = {} ({})".format(round(results['mc'][1],3), round(errors[d][1],3)))
            stats.GetListOfLines().Last().SetTextColor(ROOT.kBlue)
            stats.AddText("chi2/dof = {}".format(round(results['mc'][2],2)))
            stats.GetListOfLines().Last().SetTextColor(ROOT.kBlue)

        stats.Draw("SAME")
        c1.Update()
        c1.SaveAs("{}.png".format(plot))
        c1.SaveAs("{}.pdf".format(plot))

    if pdf:
        return ws, results['mc'][0]
    
    return results, errors, e, p


def func(x, m, x0):
     return np.sqrt((m * x)**2 + x0**2)


def plot_resol(smear, ress_mc, errs_mc, res_mc, ress_dt, errs_dt, res_dt, name):

    plt.errorbar(smear, ress_mc, errs_mc, label='mc', marker='+', linestyle="")
    eps=0.0000000001
    par_mc, pcov = opt.curve_fit(
        func, 
        smear, 
        ress_mc, 
        #bounds=((0, min(ress_mc[0], res_mc)), (np.inf, max(ress_mc[0], res_mc)))
        )
    
    plt.errorbar(smear, ress_dt, errs_dt, label='dt', marker='+', linestyle="")
    par_dt, pcov = opt.curve_fit(
        func, 
        smear, 
        ress_dt, 
        #bounds=((0, min(ress_dt[0], res_dt)), (np.inf, max(ress_dt[0], res_dt)))
        )
    x = np.linspace(0, 4, 1000)
    y_mc = func(x, *par_mc)
    y_dt = func(x, *par_dt)

    if y_dt[0] < y_mc[0]:
        x_sf = 1
    else:
        x_sf = x[np.where(y_mc > y_dt[0])[0][0]]
    print(x_sf, y_dt[0])

    plt.plot(x, y_mc, label='fit to mc')
    plt.plot(x, y_dt, label='fit to data')
    plt.plot(x, np.ones_like(x)*y_dt[0], label='data resolution')
    plt.plot(x, np.ones_like(x)*y_mc[0], label='mc resolution')
    xmin, xmax, ymin, ymax = plt.axis()
    plt.plot([x_sf, x_sf], [ymin,ymax], label='additional smearing', color="grey", linestyle="dashed")
    plt.legend()
    plt.text(2, ymax + 0.04 * (ymax - ymin), 'smearing: {}'.format(round(x_sf,3)))
    plt.ylim(ymin,ymax)
    # plt.text(0.1, 2.7, 'gradient: {}'.format(round(par_dt[0],3)))
    plt.grid()
    plt.xlabel('smearing factor of pt')
    plt.ylabel('additional smearing of m_vis')
    plt.savefig(f'plots/smeared_fits/{name}.png')
    plt.savefig(f'plots/smeared_fits/{name}.pdf')
    plt.clf()
    return x_sf, y_mc[0], y_dt[0]


def usedir(directory, overwrite = False):
    if os.path.exists(directory) and overwrite:
        print('Files in {} may be overwritten'.format(directory))
        return True

    elif os.path.exists(directory) and not overwrite:
        ow = input("Directory {} already exists. Files may be overwritten. Continue? (y/n) ".format(directory))
        if ow == 'y':
            print("Files in {} will be overwritten.".format(directory))
            return True
        else:
            print("Files in {} will not be overwritten.".format(directory))
            return False
    else:
        os.makedirs(directory)        
        print("Created new directory: ", directory)
        return True
