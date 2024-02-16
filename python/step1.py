import ROOT
from array import array
from tqdm import tqdm
import os

def hist_oneOverpT(ntuples, oneOverPt_bins, eta_bins, phi_bins, hdir, pdir, corr='')->None:
    """create histograms for inverse transversal momentum from ntuples"""

    hists = {}
    hists["tosave"] = []
    negpos = ["neg", "pos"]

    #iterate over data sources
    for typ in ntuples:
        hists[typ+'0'] = []
        hists[typ+'1'] = []

        for sample in ntuples[typ]:
            #print(sample)
            gen, sgen = "", ""
            roccor = corr
            if typ == "GEN":
                sgen = "smearedgen"
                gen = 'gen'
                roccor = ''
            #create RDataframe to acess data
            rdf = ROOT.RDataFrame("Events", ntuples[typ][sample])
            #print(ntuples[typ][sample])
            rdf = rdf.Define("weight", "zPtWeight*genWeight/sumwWeight*xsec*sf_id*sf_iso")

            for np in range(2):
                #define new column for 1/pt
                rdf = rdf.Define(f"oneOverPt_{np+1}", f"1./{sgen}pt_{np+1}{roccor}") # TODO check sgen vs gen
                #create 3D histogram
                h_3d = rdf.Histo3D(
                    (
                        f"h_oneOverPt_{sample}_{negpos[np]}", "", 
                        len(eta_bins)-1, array('f', eta_bins), 
                        len(phi_bins)-1, array('f', phi_bins), 
                        len(oneOverPt_bins)-1, array('f', oneOverPt_bins)
                    ),
                    f"{gen}eta_{np+1}",
                    f"{gen}phi_{np+1}",
                    f"oneOverPt_{np+1}",
                    "weight"
                )
                hists[typ+str(np)] += [h_3d]

    # get data-bkg histogram
    for np in range(2):
        # get histograms of Signal MC and Data and GEN
        h_mc = hists['SIG'+str(np)][0].Clone(f"h_oneOverPt_MC_{negpos[np]}")
        h_dt = hists['DATA'+str(np)][0].Clone(f"h_oneOverPt_DATA_{negpos[np]}")

        # get histogram of combined backgrounds (not normalized to data)
        h_sum_bkg = hists['BKG'+str(np)][0].Clone(f"h_oneOverPt_BKG_{negpos[np]}")
        for i in range(len(hists['BKG'+str(np)])-1):
            h_tmp = hists['BKG'+str(np)][i+1].Clone("h_tmp")
            h_sum_bkg.Add(h_tmp)
        
        # calculate combined MC contribution to extract MC->Data normalization SF
        h_mc.Add(h_sum_bkg)
        sf = h_dt.Integral() / h_mc.Integral()

        # scale backgrounds to data
        h_sum_bkg.Scale(sf)
        h_dt.Add(h_sum_bkg, -1)

        hists["tosave"] += [h_dt, hists['SIG'+str(np)][0], hists['GEN'+str(np)][0]]

        # get median
        for h in hists["tosave"][-3:]:
            hists["tosave"] += [h.Project3DProfile("yx")]            
        
    #save
    tf = ROOT.TFile(f"{hdir}step1_oneOverPt{corr}.root","RECREATE")
    for h in hists["tosave"]:
        h.Write()
    tf.Close()


def get_scale_corrections(samples, eta_bins, phi_bins, charge_bins, hdir)->None: 
    """extract scale corrections from ntuple data"""

    #charges
    negpos = ["neg", "pos"]
    
    # get 3D histograms from TFile
    tf = ROOT.TFile(f"{hdir}step1_oneOverPt.root", "READ")
    oneOverPt_hists = {}

    #iterate over data, mc and gen
    for typ in samples:
        for np in negpos:
            #read mean-histogram
            oneOverPt_hists[f"{typ}_{np}"] = tf.Get(f'h_oneOverPt_{typ}_{np}_pyx')
            oneOverPt_hists[f"{typ}_{np}"].SetDirectory(ROOT.nullptr)
    tf.Close()

    # define correction factor C from paper as root 3d histogram
    C, Dm, Da, M, A = {}, {}, {}, {}, {}
    for typ in samples[:-1]:
        #print(s)
        #make 3D-histograms for C, Dm, Da, M and A
        C[typ] = ROOT.TH3D(
            f"C_{typ}", "",
            len(charge_bins)-1, array('f', charge_bins),
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        Dm[typ] = ROOT.TH2D(
            f"Dm_{typ}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        Da[typ] = ROOT.TH2D(
            f"Da_{typ}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        M[typ] = ROOT.TH2D(
            f"M_{typ}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        A[typ] = ROOT.TH2D(
            f"A_{typ}", "",
            len(eta_bins)-1, array('f', eta_bins),
            len(phi_bins)-1, array('f', phi_bins)
        )
        #iterate over eta, phi and charge
        for eta in range(len(eta_bins)-1):
            for phi in range(len(phi_bins)-1):
                for charge in range(len(charge_bins)-1):
                    mean_gen = oneOverPt_hists[f"GEN_{negpos[charge]}"].GetBinContent(eta+1, phi+1)
                    mean = oneOverPt_hists[f"{typ}_{negpos[charge]}"].GetBinContent(eta+1, phi+1)

                    #set value of C[mode+s] for current bin
                    C[typ].SetBinContent(
                        charge+1,
                        eta+1,
                        phi+1,
                        mean_gen - mean
                    )
                    #print(mean_gen, mean)
                #calculate Dm, Da, M, A for current bin
                Dm[typ].SetBinContent(
                    eta+1,
                    phi+1,
                    (C[typ].GetBinContent(2, eta+1, phi+1) + C[typ].GetBinContent(1, eta+1, phi+1)) / 2.
                )
                Da[typ].SetBinContent(
                    eta+1,
                    phi+1,
                    (C[typ].GetBinContent(2, eta+1, phi+1) - C[typ].GetBinContent(1, eta+1, phi+1)) / 2.
                )
                M[typ].SetBinContent(
                    eta+1,
                    phi+1,
                    1 + 2*Dm[typ].GetBinContent(eta+1, phi+1) / (
                        oneOverPt_hists[f"{typ}_{negpos[0]}"].GetBinContent(eta+1, phi+1) + 
                        oneOverPt_hists[f"{typ}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                    )
                )

                A[typ].SetBinContent(
                    eta+1,
                    phi+1,
                    Da[typ].GetBinContent(eta+1, phi+1) - (
                        Dm[typ].GetBinContent(eta+1, phi+1)*(
                        oneOverPt_hists[f"{typ}_{negpos[0]}"].GetBinContent(eta+1, phi+1) -
                        oneOverPt_hists[f"{typ}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                    ) / (
                        oneOverPt_hists[f"{typ}_{negpos[0]}"].GetBinContent(eta+1, phi+1) +
                        oneOverPt_hists[f"{typ}_{negpos[1]}"].GetBinContent(eta+1, phi+1)
                    )
                    )
                )

    #safe histograms as Tfile
    tf = ROOT.TFile(f"{hdir}step1_C.root", "RECREATE")
    for typ in samples[:-1]:
        C[typ].Write()
        Dm[typ].Write()
        Da[typ].Write()
        M[typ].Write()
        A[typ].Write()
    tf.Close()