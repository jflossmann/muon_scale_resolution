import ROOT
import numpy as np
from array import array
from tqdm import tqdm
from multiprocessing import Pool, RLock
import os

#Dictionary for names and variables of histograms
plotdict = {
    'MC': {
        'h_EtaPhiVsMz_neg_MC': ['eta_1', 'phi_1', 'mass_Z'],
        'h_EtaPhiVsMz_pos_MC': ['eta_2', 'phi_2', 'mass_Z'],
        'h_EtaPhiVsMz_neg_mean_roccor_MC': ['eta_1', 'phi_1', 'mass_Z_mean_roccor'],
        'h_EtaPhiVsMz_pos_mean_roccor_MC': ['eta_2', 'phi_2', 'mass_Z_mean_roccor'],
        'h_EtaPhiVsMz_neg_median_roccor_MC': ['eta_1', 'phi_1', 'mass_Z_median_roccor'],
        'h_EtaPhiVsMz_pos_median_roccor_MC': ['eta_2', 'phi_2', 'mass_Z_median_roccor'],
        'h_EtaPhiVsMz_neg_peak_roccor_MC': ['eta_1', 'phi_1', 'mass_Z_peak_roccor'],
        'h_EtaPhiVsMz_pos_peak_roccor_MC': ['eta_2', 'phi_2', 'mass_Z_peak_roccor'],
        'h_EtaPhiVsMz_neg_GEN': ['geneta_1', 'genphi_1', 'genmass_Z'],
        'h_EtaPhiVsMz_pos_GEN': ['geneta_2', 'genphi_2', 'genmass_Z'],
    },
    'DATA': {
        'h_EtaPhiVsMz_neg_DATA': ['eta_1', 'phi_1', 'mass_Z'],
        'h_EtaPhiVsMz_pos_DATA': ['eta_2', 'phi_2', 'mass_Z'],
        'h_EtaPhiVsMz_neg_mean_roccor_DATA': ['eta_1', 'phi_1', 'mass_Z_mean_roccor'],
        'h_EtaPhiVsMz_pos_mean_roccor_DATA': ['eta_2', 'phi_2', 'mass_Z_mean_roccor'],
        'h_EtaPhiVsMz_neg_median_roccor_DATA': ['eta_1', 'phi_1', 'mass_Z_median_roccor'],
        'h_EtaPhiVsMz_pos_median_roccor_DATA': ['eta_2', 'phi_2', 'mass_Z_median_roccor'],
        'h_EtaPhiVsMz_neg_peak_roccor_DATA': ['eta_1', 'phi_1', 'mass_Z_peak_roccor'],
        'h_EtaPhiVsMz_pos_peak_roccor_DATA': ['eta_2', 'phi_2', 'mass_Z_peak_roccor'],        
    }
}


def roofit_mass(hist, e, p):
    """fits a Voigt function to a given histogram and returns result and errors"""
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

    results, errors, plots = {}, {}, {}

    #definition of mass-variable
    x = ROOT.RooRealVar("x", "m_vis (GeV)", 85, 97)
    x.setBins(10000,"cache")
    x.setMin("cache",0)
    x.setMax("cache",500)

    #definition of fit parameters
    Z_mass = ROOT.RooRealVar("Z_mass", "Z_mass", 91.1876, 60, 120)
    Z_width = ROOT.RooRealVar("Z_width", "Z_width", 2.4952, 2.4952, 2.4952)
    Z_width.setConstant(True)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 0, 10)
    n_sigma=1

    #creation of voigt-model
    func = ROOT.RooVoigtian("vgt_mc", "Voigt", x, Z_mass, Z_width, sigma)

    #creation of a RooDataHist from histogram
    roohist = ROOT.RooDataHist('h', '', ROOT.RooArgSet(x), hist)
    #perform fit
    fitResult = func.fitTo(roohist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintEvalErrors(-1))
    
    results = [Z_mass.getVal(), n_sigma*sigma.getVal()]
    errors = [Z_mass.getAsymErrorHi(), n_sigma*sigma.getAsymErrorHi()]

    return results, errors, e, p



def hist_zmass(ntuples, eta_bins, phi_bins, mass_bins, hdir)->None:
    """produces 3D and 1D histograms from given ntuples"""
    hists_3d = []
    hists_1d = []
    for s in ntuples:
        rdf = ROOT.RDataFrame("Events", ntuples[s])

        for n in plotdict[s]:
            #create 3D histogram from ntuple
            hists_3d.append(
                rdf.Histo3D(
                    (
                        n, "",
                        len(eta_bins)-1, array('d', eta_bins),
                        len(phi_bins)-1, array('d', phi_bins),
                        len(mass_bins)-1, array('d', mass_bins),
                    ),
                plotdict[s][n][0], plotdict[s][n][1], plotdict[s][n][2],
                'zPtWeight' 
                )
            )

            # extract bin by bin mass histograms from 3d
            for e in range(len(eta_bins)-1):
                for p in range(len(phi_bins)-1):
                    h = hists_3d[-1].ProjectionZ(f"{n}_eta{e}_phi{p}", e+1, e+1, p+1, p+1)
                    hists_1d.append(h)
                    
    tf = ROOT.TFile(f'{hdir}zmass.root', 'recreate')
    for h in hists_3d + hists_1d:
        h.Write()
    tf.Close()


def fit_zmass(eta_bins, phi_bins, hdir):
    """performs fit for z-mass histogram. produces a new histogram for the fit results"""
    hists = []

    tf = ROOT.TFile(f'{hdir}zmass.root', 'read')
    for s in plotdict:
        for n in plotdict[s]:
            print(n)
            h_2d = ROOT.TH2D(
                f"{n}_fitresult", "", 
                len(eta_bins)-1, array('d', eta_bins), 
                len(phi_bins)-1, array('d', phi_bins)
            )
            arguments = []
            for e in range(len(eta_bins)-1):
                for p in range(len(phi_bins)-1):
                    h = tf.Get(f"{n}_eta{e}_phi{p}")
                    h.SetDirectory(ROOT.nullptr)
                    arguments.append((h, e, p))
            pool = Pool(8, initargs=(RLock(),), initializer=tqdm.set_lock)
            for results in tqdm(pool.imap_unordered(job_wrapper, arguments)):
                fit, fit_err, e, p = results            
                h_2d.SetBinContent(e+1, p+1, fit[0])
                h_2d.SetBinError(e+1, p+1, fit[1])
            
            h_2d.SetDirectory(ROOT.nullptr)
            hists.append(h_2d)

    tf.Close()

    tf = ROOT.TFile(f'{hdir}zmass_fitresults.root', 'recreate')
    for h in hists:
        h.Write()
    tf.Close()


def job_wrapper(args):
    return roofit_mass(*args)


def plot_zmass(eta_bins, phi_bins, hdir):
    """produces plots for fit results and saves them in "plots/zmass_fits/" """
    ROOT.gROOT.SetBatch(1)
    tf = ROOT.TFile(f'{hdir}zmass_fitresults.root', 'read')
    os.makedirs(f'plots/zmass_fits/', exist_ok=True)

    for np in ['neg', 'pos']:
        for roccor in ['', '_mean_roccor', '_median_roccor', '_peak_roccor']:
            h_mc = tf.Get(f'h_EtaPhiVsMz_{np+roccor}_MC_fitresult')
            h_gen = tf.Get(f'h_EtaPhiVsMz_{np}_GEN_fitresult')
            h_data = tf.Get(f'h_EtaPhiVsMz_{np+roccor}_DATA_fitresult')

            h_mc.SetDirectory(ROOT.nullptr)
            h_gen.SetDirectory(ROOT.nullptr)
            h_data.SetDirectory(ROOT.nullptr)

            h1_mc = h_mc.ProjectionX()
            h1_gen = h_gen.ProjectionX()
            h1_data = h_data.ProjectionX()

            for e in range(len(eta_bins)-1):
                h1_mc.SetBinContent(e+1, h1_mc.GetBinContent(e+1)/(len(phi_bins)-1))
                h1_mc.SetBinError(e+1, h1_mc.GetBinError(e+1)/(len(phi_bins)-1))
                h1_gen.SetBinContent(e+1, h1_gen.GetBinContent(e+1)/(len(phi_bins)-1))
                h1_gen.SetBinError(e+1, h1_gen.GetBinError(e+1)/(len(phi_bins)-1))
                h1_data.SetBinContent(e+1, h1_data.GetBinContent(e+1)/(len(phi_bins)-1))
                h1_data.SetBinError(e+1, h1_data.GetBinError(e+1)/(len(phi_bins)-1))

            ROOT.gStyle.SetOptStat(0)

            c = ROOT.TCanvas()

            h1_mc.Draw()
            h1_mc.SetLineColor(ROOT.kBlue)
            h1_mc.SetLineWidth(2)

            h1_gen.Draw("SAME")
            h1_gen.SetLineColor(ROOT.kGreen)
            h1_gen.SetLineWidth(2)

            h1_data.Draw("same")
            h1_data.SetLineColor(ROOT.kBlack)
            h1_data.SetLineWidth(2)

            tl = ROOT.TLegend()
            tl.AddEntry(h1_mc, 'MC')
            tl.AddEntry(h1_gen, 'GEN')
            tl.AddEntry(h1_data, 'Data')
            tl.Draw("same")

            c.SaveAs(f'plots/zmass_fits/zmass_fits_{np+roccor}.png')
    
