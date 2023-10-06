import ROOT
from array import array
from tqdm import tqdm
from multiprocessing import Pool, RLock
import os

#define different statistics to base corrections on
modes = ['mean', 'median', 'peak']


def fit_oneOverpT(hist, plot, e, p):
    """fits a Crystall Ball fuction to 1/pt histogram"""
    
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.gROOT.SetBatch(1)

    #make a variable object for independent variable
    x = ROOT.RooRealVar("x", "oneOverPt (1/GeV)", 1/200, 1/20)
    x.setBins(10000,"cache")
    x.setMin("cache",0)
    x.setMax("cache",500)

    #make canvas and frame for plot
    c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
    c1.SetGridx()
    c1.SetGridy()
    c1.GetFrame().SetFillColor( 21 )
    c1.GetFrame().SetBorderMode(-1 )
    c1.GetFrame().SetBorderSize( 5 )
    frame = x.frame()
    frame.SetTitle('')

    #create parameter objects for the CB-function
    mean = ROOT.RooRealVar("mean", "mean", 0, 1/50, 1./40)
    #cb_mean.setConstant(True)
    sigma_L = ROOT.RooRealVar("sigma_L", "sigma", 0.0007, 0, 0.001)
    sigma_R = ROOT.RooRealVar("sigma_R", "sigma", 0.001, 0, 0.001)
    n_sigma = 1 # since sigma_total = sigma_left + sigma_right
    n_L = ROOT.RooRealVar("n_L", "n_L", 5, 0, 1000)
    n_R = ROOT.RooRealVar("n_R", "n_R", 5, 0, 1000)
    alpha_L = ROOT.RooRealVar("alpha_L", "alpha_L", 0.3, 0, 2)
    alpha_R = ROOT.RooRealVar("alpha_R", "alpha_R", 0.3, 0, 2)

    #create CB function
    func =  ROOT.RooCrystalBall("cb", "CrystalBall", x, mean, sigma_L,
                            sigma_R, alpha_L, n_L, alpha_R, n_R)

    #convert histogram to RooDataHist-Object
    roohist = ROOT.RooDataHist('h', '', ROOT.RooArgSet(x), hist)
    #perform fit
    fitResult = func.fitTo(roohist, ROOT.RooFit.AsymptoticError(True), ROOT.RooFit.PrintEvalErrors(-1))

    #draw histogram and fit on frame
    roohist.plotOn(frame, ROOT.RooFit.DrawOption("B"), ROOT.RooFit.FillStyle(0), ROOT.RooFit.FillColor(ROOT.kBlue))
    func.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue))

    #draw frame and safe plot
    frame.Draw()
    c1.Update()
    ROOT.gStyle.SetGridColor(ROOT.kGray+1)
    c1.SaveAs("{}.png".format(plot))
    c1.SaveAs("{}.pdf".format(plot))
    
    print(sigma_L.getVal(), sigma_R.getVal())
    
    return mean.getVal(), e, p
    

def job_wrapper(args):
    return fit_oneOverpT(*args)


def hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='')->None:
    """create histograms for inverse transversal momentum from ntuple data"""
    
    os.makedirs(pdir+'fits/', exist_ok=True)
    hists = []
    ntuples_tmp = dict(ntuples)
    ntuples_tmp["GEN"] = ntuples_tmp["MC"]
    negpos = ["neg", "pos"]

    #iterate over data sources
    for s in ntuples_tmp:
        os.makedirs(pdir+f'fits/{s+corr}/', exist_ok=True)
        gen = ""
        roccor = corr
        if s == "GEN":
            gen = "gen"
            roccor = ''
        #create RDataframe to acess data
        rdf = ROOT.RDataFrame("Events", ntuples_tmp[s])

        for np in range(2):
            #define new column for 1/pt
            rdf = rdf.Define(f"oneOverPt_{np+1}", f"1./{gen}pt_{np+1}{roccor}")
            #create 3D histogram
            h_3d = rdf.Histo3D(
                (
                    f"h_oneOverPt_{s}_{negpos[np]}", "", 
                    len(eta_bins)-1, array('f', eta_bins), 
                    len(phi_bins)-1, array('f', phi_bins), 
                    len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
                ),
                f"{gen}eta_{np+1}",
                f"{gen}phi_{np+1}",
                f"oneOverPt_{np+1}",
                "zPtWeight" # TODO: improve method. averaging over bins not precise enough
            )
            #create histogram for median and peak
            h_median = ROOT.TH2D(
                f"h_oneOverPt_{s}_median_{negpos[np]}", "",
                len(eta_bins)-1, array('d', eta_bins), 
                len(phi_bins)-1, array('d', phi_bins)
            )

            h_peak = ROOT.TH2D(
                f"h_oneOverPt_{s}_peak_{negpos[np]}", "",
                len(eta_bins)-1, array('d', eta_bins), 
                len(phi_bins)-1, array('d', phi_bins)
            )

            quantile = array('d', [0.5])
            median = array('d', [0])

            arguments = []
            for e in range(len(eta_bins)-1):
                for p in range(len(phi_bins)-1):
                    h_ep = h_3d.ProjectionZ(f"h_tmp_{e}_{p}", e+1, e+1, p+1, p+1)

                    arguments.append((h_ep, f'{pdir}fits/{s+corr}/{e}_{p}', e, p))

                    h_ep.GetQuantiles(1, median, quantile)
                    h_median.SetBinContent(e+1, p+1, median[0])

            pool = Pool(8, initargs=(RLock(),), initializer=tqdm.set_lock)
            for results in tqdm(pool.imap_unordered(job_wrapper, arguments)):
                fit, e, p = results            
                h_peak.SetBinContent(e+1, p+1, fit)
                h_peak.SetBinError(e+1, p+1, fit)

            hists += [h_3d, h_3d.Project3DProfile("yx"), h_median, h_peak]
    #save
    tf = ROOT.TFile(f"{hdir}oneOverPt{corr}.root","RECREATE")
    for h in hists:
        h.Write()
    tf.Close()



def get_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir)->None: 
    """extract scale corrections from ntuple data"""

    #charges
    negpos = ["neg", "pos"]
    
    # get 3D histograms from TFile
    tf = ROOT.TFile(f"{hdir}oneOverPt.root", "READ")
    oneOverPt_hists = {}

    #iterate over data, mc and gen
    for s in list(ntuples.keys())+["GEN"]:
        for np in negpos:
            #read mean-histogram
            oneOverPt_hists[f"{s}_mean_{np}"] = tf.Get(f'h_oneOverPt_{s}_{np}_pyx')
            oneOverPt_hists[f"{s}_mean_{np}"].SetDirectory(ROOT.nullptr)
            #read median histogram
            oneOverPt_hists[f"{s}_median_{np}"] = tf.Get(f'h_oneOverPt_{s}_median_{np}')
            oneOverPt_hists[f"{s}_median_{np}"].SetDirectory(ROOT.nullptr)
            #rean peak histogram
            oneOverPt_hists[f"{s}_peak_{np}"] = tf.Get(f'h_oneOverPt_{s}_peak_{np}')
            oneOverPt_hists[f"{s}_peak_{np}"].SetDirectory(ROOT.nullptr)
    tf.Close()

    # define correction factor C from paper as root 3d histogram
    C, Dm, Da, M, A = {}, {}, {}, {}, {}
    for mode in modes:
        for s in ntuples:
            #print(s)
            #make 3D-histograms for C, Dm, Da, M and A
            C[mode+s] = ROOT.TH3D(
                f"C_{mode}_{s}", "",
                len(charge_bins)-1, array('f', charge_bins),
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            Dm[mode+s] = ROOT.TH2D(
                f"Dm_{mode}_{s}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            Da[mode+s] = ROOT.TH2D(
                f"Da_{mode}_{s}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            M[mode+s] = ROOT.TH2D(
                f"M_{mode}_{s}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            A[mode+s] = ROOT.TH2D(
                f"A_{mode}_{s}", "",
                len(eta_bins)-1, array('f', eta_bins),
                len(phi_bins)-1, array('f', phi_bins)
            )
            #iterate over eta, phi and charge
            for eta in range(len(eta_bins)-1):
                for phi in range(len(phi_bins)-1):
                    for charge in range(len(charge_bins)-1):
                        mean_gen = oneOverPt_hists[f"GEN_{mode}_{negpos[charge]}"].GetBinContent(eta+1, phi+1)
                        mean = oneOverPt_hists[f"{s}_{mode}_{negpos[charge]}"].GetBinContent(eta+1, phi+1)

                        #set value of C[mode+s] for current bin
                        C[mode+s].SetBinContent(
                            charge+1,
                            eta+1,
                            phi+1,
                            mean_gen - mean
                        )
                        #print(mean_gen, mean)
                    #calculate Dm, Da, M, A for current bin
                    Dm[mode+s].SetBinContent(
                        eta+1,
                        phi+1,
                        (C[mode+s].GetBinContent(2, eta+1, phi+1) + C[mode+s].GetBinContent(1, eta+1, phi+1)) / 2.
                    )
                    Da[mode+s].SetBinContent(
                        eta+1,
                        phi+1,
                        (C[mode+s].GetBinContent(2, eta+1, phi+1) - C[mode+s].GetBinContent(1, eta+1, phi+1)) / 2.
                    )
                    M[mode+s].SetBinContent(
                        eta+1,
                        phi+1,
                        1 + 2*Dm[mode+s].GetBinContent(eta+1, phi+1) / (
                            oneOverPt_hists[f"{s}_{mode}_{negpos[0]}"].GetBinContent(eta+1, phi+1) + 
                            oneOverPt_hists[f"{s}_{mode}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                        )
                    )

                    A[mode+s].SetBinContent(
                        eta+1,
                        phi+1,
                        Da[mode+s].GetBinContent(eta+1, phi+1) - (
                            Dm[mode+s].GetBinContent(eta+1, phi+1)*(
                            oneOverPt_hists[f"{s}_{mode}_{negpos[0]}"].GetBinContent(eta+1, phi+1) -
                            oneOverPt_hists[f"{s}_{mode}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                        ) / (
                            oneOverPt_hists[f"{s}_{mode}_{negpos[0]}"].GetBinContent(eta+1, phi+1) +
                            oneOverPt_hists[f"{s}_{mode}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                        )
                        )
                    )
    #safe histograms as Tfile
    tf = ROOT.TFile(f"{hdir}C.root", "RECREATE")
    for mode in modes:
        for s in ntuples:
            C[mode+s].Write()
            Dm[mode+s].Write()
            Da[mode+s].Write()
            M[mode+s].Write()
            A[mode+s].Write()
    tf.Close()



def apply_scale_corrections(ntuples, eta_bins, phi_bins, charge_bins, hdir):
    #open Tfile with scale corrections
    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}C.root", "READ");')
    for mode in modes:
        #read histograms for curent mode
        ROOT.gROOT.ProcessLine(f'TH2D* M_{mode}_DATA = (TH2D*)tf->Get("M_{mode}_DATA");')
        ROOT.gROOT.ProcessLine(f'TH2D* M_{mode}_MC = (TH2D*)tf->Get("M_{mode}_MC");')
        ROOT.gROOT.ProcessLine(f'TH2D* A_{mode}_DATA = (TH2D*)tf->Get("A_{mode}_DATA");')
        ROOT.gROOT.ProcessLine(f'TH2D* A_{mode}_MC = (TH2D*)tf->Get("A_{mode}_MC");')   
    
    for s in ntuples:
        #make RDataFrame object for current dataset
        rdf = ROOT.RDataFrame("Events", ntuples[s])
        quants = list(rdf.GetColumnNames())
        #calculate reciprocal of momentum for both muons
        rdf = rdf.Define("oneOverPt_1", "1./pt_1")
        rdf = rdf.Define("oneOverPt_2", "1./pt_2")

        
        for mode in modes:
            #calculate the corrected reciprocal pt
            rdf = rdf.Define(
                f"oneOverPt_1_{mode}_roccor",
                f"""oneOverPt_1 * M_{mode}_{s}->GetBinContent(M_{mode}_{s}->FindBin(eta_1, phi_1)) + 
                charge_1 * A_{mode}_{s}->GetBinContent(A_{mode}_{s}->FindBin(eta_1, phi_1))
                """
            )
            rdf = rdf.Define(f"pt_1_{mode}_roccor", f"1./oneOverPt_1_{mode}_roccor")
            rdf = rdf.Define(
                f"oneOverPt_2_{mode}_roccor",
                f"""oneOverPt_2 * M_{mode}_{s}->GetBinContent(M_{mode}_{s}->FindBin(eta_2, phi_2)) + 
                charge_2 * A_{mode}_{s}->GetBinContent(A_{mode}_{s}->FindBin(eta_2, phi_2))
                """
            )
            rdf = rdf.Define(f"pt_2_{mode}_roccor", f"1./oneOverPt_2_{mode}_roccor")

            #calculate corrected 4-momenta and corrected Z-mass
            rdf = rdf.Define(f"p4_1_{mode}", f"ROOT::Math::PtEtaPhiMVector(pt_1_{mode}_roccor, eta_1, phi_1, mass_1)")
            rdf = rdf.Define(f"p4_2_{mode}", f"ROOT::Math::PtEtaPhiMVector(pt_2_{mode}_roccor, eta_2, phi_2, mass_2)")
            rdf = rdf.Define(f"p4_Z_{mode}", f"p4_1_{mode} + p4_2_{mode}")
            rdf = rdf.Define(f"mass_Z_{mode}_roccor", f"p4_Z_{mode}.M()")

            #add calculated quantities to quantities
            quants += [f"pt_1_{mode}_roccor", f"pt_2_{mode}_roccor", f"mass_Z_{mode}_roccor"]
            
        #save everything in new ntuple
        rdf.Snapshot("Events", ntuples[s].replace('.root', '_corr.root'), quants)
