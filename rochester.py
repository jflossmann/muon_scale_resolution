from Higgs_Mass_setup import DrawResolution
import ROOT
import array

#TString USER_DIR = "/afs/cern.ch/work/f/ferrico/private/Rochester/CMSSW_12_4_3/src/";
#TString USER_DIR = "/afs/cern.ch/work/x/xzuo/UF_NanoAODTool_10619p2/src/Roch_test/";
USER_DIR="/afs/cern.ch/work/x/xzuo/UF_NanoAODTool_10619p2/src/Roch_test/"

print(USER_DIR)

def rochester():
# Initialisieren von ROOT und setzen Sie den Batch-Modus
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)

    Anno_2018 = True
    lumi = 86.6855
    if Anno_2018:
        lumi = 373

    reweight = 1

    _file0 = ROOT.TFile()  # Sie müssen den Dateinamen hier hinzufügen, z.B., "my_file.root"
    tree = ROOT.TTree()
    nentries = 0
    rand = ROOT.TRandom3()
    u1 = 0

    crossSection = 13730  # pb (aus DAS)
    weight = 1

    samples = ["MC", "DATA", "GEN"]

    pt_bins = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 100, 140, 200]

    eta_bins = [-2.4, 2.4]

    phi_bins = [-3.2, 3.2]

    qOverPt_bins = [i / 200.0 for i in range(21)]

    mass_bins = [60 + i for i in range(61)]

    h_Zboson_mass = [ROOT.TH1F() for _ in range(3)]
    h_Zboson_pT = [ROOT.TH1F() for _ in range(3)]
    h_lepton_qOverPt = [ROOT.TH1F() for _ in range(3)]
    h_lepton_qOverPt_scan = [[[ROOT.TH1F() for _ in range(40)] for _ in range(40)] for _ in range(3)]
    h_lepton_pt = [ROOT.TH1F() for _ in range(3)]
    h_lepton_eta = [ROOT.TH1F() for _ in range(3)]
    h_lepton_phi = [ROOT.TH1F() for _ in range(3)]

    h_Zboson_mass_vs_phi = [[[ROOT.TH2F() for _ in range(2)] for _ in range(3)] for _ in range(3)]
    h_Zboson_mass_vs_eta = [[[ROOT.TH2F() for _ in range(2)] for _ in range(3)] for _ in range(3)]

    h_Zboson_mass_EtaPhi = [ROOT.TH3F() for _ in range(3)]
    h_lepton_qOverPt_scan_EtaPhi = [[[ROOT.TH3F() for _ in range(2)] for _ in range(3)] for _ in range(3)]

    h_pTRes = ROOT.TH2F("h_pTRes", "h_pTRes", len(pt_bins) - 1, array('f', pt_bins), 100, -0.1, 0.1)  # Sie müssen das 'array'-Modul importieren

    ############################################################################################

    for i in range(len(samples)):
        histo_name = "h_Zboson_mass_" + samples[i]
        h_Zboson_mass[i] = ROOT.TH1F(histo_name, histo_name, 80, 81, 101)
        histo_name = "h_Zboson_pT_" + samples[i]
        h_Zboson_pT[i] = ROOT.TH1F(histo_name, histo_name, 80, 0.0, 160.0)
        histo_name = "h_lepton_qOverPt_" + samples[i]
        h_lepton_qOverPt[i] = ROOT.TH1F(histo_name, histo_name, 50, -0.1, 0.1)
        histo_name = "h_lepton_pt_" + samples[i]
        h_lepton_pt[i] = ROOT.TH1F(histo_name, histo_name, 40, 0.0, 200.0)
        histo_name = "h_lepton_eta_" + samples[i]
        h_lepton_eta[i] = ROOT.TH1F(histo_name, histo_name, 24, -2.4, 2.4)
        histo_name = "h_lepton_phi_" + samples[i]
        h_lepton_phi[i] = ROOT.TH1F(histo_name, histo_name, 1, -3.2, 3.2)

        histo_name = "h_Zboson_mass_EtaPhi_" + samples[i]
        h_Zboson_mass_EtaPhi[i] = ROOT.TH3F(histo_name, histo_name, len(eta_bins) - 1, array('f', eta_bins), len(phi_bins) - 1, array('f', phi_bins), len(mass_bins) - 1, array('f', mass_bins))

        for q in range(2):
            if q == 0:
                histo_name = "h_lepton_qOverPt_scan_EtaPhi_%s_positive" % samples[i]
            else:
                histo_name = "h_lepton_qOverPt_scan_EtaPhi_%s_negative" % samples[i]
            h_lepton_qOverPt_scan_EtaPhi[i][q] = ROOT.TH3F(histo_name, histo_name, len(eta_bins) - 1, array('f', eta_bins), len(phi_bins) - 1, array('f', phi_bins), len(qOverPt_bins) - 1, array('f', qOverPt_bins))

            if q == 0:
                histo_name = "h_Zboson_mass_vs_phi_%s_positive" % samples[i]
            else:
                histo_name = "h_Zboson_mass_vs_phi_%s_negative" % samples[i]
            h_Zboson_mass_vs_phi[i][q] = ROOT.TH2F(histo_name, histo_name, len(phi_bins) - 1, array('f', phi_bins), 80, 81, 101)

            if q == 0:
                histo_name = "h_Zboson_mass_vs_eta_%s_positive" % samples[i]
            else:
                histo_name = "h_Zboson_mass_vs_eta_%s_negative" % samples[i]
            h_Zboson_mass_vs_eta[i][q] = ROOT.TH2F(histo_name, histo_name, len(eta_bins) - 1, array('f', eta_bins), 80, 81, 101)

    for j in range(40):
        for h in range(40):
            for i in range(len(samples)):
                histo_name = "h_lepton_qOverPt_%s_%d_%d" % (samples[i], j, h)
                h_lepton_qOverPt_scan[i][j][h] = ROOT.TH1F(histo_name, histo_name, 50, -0.1, 0.1)

    muon_pt = ROOT.std.vector('double')()
    muon_eta = ROOT.std.vector('double')()
    muon_phi = ROOT.std.vector('double')()
    muon_mass = ROOT.std.vector('double')()
    muon_charge = ROOT.std.vector('double')()
    muon_isLoose = ROOT.std.vector('double')()
    muon_isMedium = ROOT.std.vector('double')()
    muon_isTight = ROOT.std.vector('double')()
    muon_trackIso = ROOT.std.vector('double')()

    GEN_pt = ROOT.std.vector('double')()
    GEN_eta = ROOT.std.vector('double')()
    GEN_phi = ROOT.std.vector('double')()
    GEN_mass = ROOT.std.vector('double')()
    GEN_id = ROOT.std.vector('double')()
    GEN_status = ROOT.std.vector('double')()
    GEN_motherId = ROOT.std.vector('double')()

    hlt_pt = ROOT.std.vector('double')()
    hlt_eta = ROOT.std.vector('double')()
    hlt_phi = ROOT.std.vector('double')()

    passedTrig = False

    for s in range(2):
        print("step1", USER_DIR + samples[s])
        if Anno_2018:
            _file0 = ROOT.TFile(USER_DIR + samples[s] + ".root")

        if _file0:
            tree = _file0.Get("Ana/passedEvents")
        else:
            print("File not found")
            return

        if not tree:
            print("tree not found")
            return

        h_num_eventi = _file0.Get("Ana/sumWeights")
        num_event = h_num_eventi.Integral()

        tree.SetBranchAddress("muon_pt", muon_pt)
        tree.SetBranchAddress("muon_eta", muon_eta)
        tree.SetBranchAddress("muon_phi", muon_phi)
        tree.SetBranchAddress("muon_mass", muon_mass)
        tree.SetBranchAddress("muon_charge", muon_charge)
        tree.SetBranchAddress("muon_isLoose", muon_isLoose)
        tree.SetBranchAddress("muon_isMedium", muon_isMedium)
        tree.SetBranchAddress("muon_isTight", muon_isTight)
        tree.SetBranchAddress("muon_trackIso", muon_trackIso)

        tree.SetBranchAddress("GEN_pt", GEN_pt)
        tree.SetBranchAddress("GEN_eta", GEN_eta)
        tree.SetBranchAddress("GEN_phi", GEN_phi)
        tree.SetBranchAddress("GEN_mass", GEN_mass)
        tree.SetBranchAddress("GEN_id", GEN_id)
        tree.SetBranchAddress("GEN_status", GEN_status)
        tree.SetBranchAddress("GEN_motherId", GEN_motherId)

        tree.SetBranchAddress("hlt_pt", hlt_pt)
        tree.SetBranchAddress("hlt_eta", hlt_eta)
        tree.SetBranchAddress("hlt_phi", hlt_phi)

        tree.SetBranchAddress("passedTrig", passedTrig)

        nentries = tree.GetEntries()

        print(nentries)

        for entry in range(nentries):
            tree.GetEntry(entry)
            if entry == 0:
                print(weight)

            if entry % 250000 == 0:
                print(entry, "/", nentries)

            # if not passedTrig:
            #     continue

            if len(muon_pt) < 2:
                continue

            muon_hlt_match = []
            muon_GEN_match = []

            for j in range(len(muon_pt)):
                if len(muon_hlt_match) > 1:
                    continue

                deltaR = 1000
                if muon_pt[j] < 10:
                    continue
                if abs(muon_eta[j]) > 2.4:
                    continue
                if muon_trackIso[j] / muon_pt[j] > 0.1:
                    continue
                if not muon_isMedium[j]:
                    continue
                # if not muon_isMedium[j]:
                #     continue
                # if not muon_isLoose[j]:
                #     continue

                for i in range(len(hlt_pt)):
                    deltaEta = hlt_eta[i] - muon_eta[j]
                    deltaPhi = hlt_phi[i] - muon_phi[j]
                    tmp_deltaR = (deltaEta ** 2 + deltaPhi ** 2) ** 0.5
                    if tmp_deltaR < deltaR:
                        deltaR = tmp_deltaR
                if deltaR < 0.35:
                    muon_hlt_match.append(j)

                deltaR = 1000
                k = -1
                for i in range(len(GEN_pt)):
                    if GEN_status[i] != 1:
                        continue
                    if abs(GEN_id[i]) != 13:
                        continue
                    if GEN_motherId[i] != 23:
                        continue

                    deltaEta = GEN_eta[i] - muon_eta[j]
                    deltaPhi = GEN_phi[i] - muon_phi[j]
                    tmp_deltaR = (deltaEta ** 2 + deltaPhi ** 2) ** 0.5
                    if tmp_deltaR < deltaR:
                        deltaR = tmp_deltaR
                        if deltaR < 0.1:
                            k = i
                if k != -1:
                    muon_GEN_match.append(k)

            if len(muon_hlt_match) < 2:
                continue
            if muon_pt[muon_hlt_match[0]] < 25:
                continue
            if muon_charge[muon_hlt_match[0]] + muon_charge[muon_hlt_match[1]] != 0:
                continue

            lep0 = ROOT.TLorentzVector()
            lep1 = ROOT.TLorentzVector()
            lep0.SetPtEtaPhiM(muon_pt[muon_hlt_match[0]], muon_eta[muon_hlt_match[0]], muon_phi[muon_hlt_match[0]], muon_mass[muon_hlt_match[0]])
            lep1.SetPtEtaPhiM(muon_pt[muon_hlt_match[1]], muon_eta[muon_hlt_match[1]], muon_phi[muon_hlt_match[1]], muon_mass[muon_hlt_match[1]])
            Z = lep0 + lep1

            if 81 < Z.M() < 101:
                if s != 1:
                    weight = lumi * crossSection / num_event
                else:
                    weight = 1

                h_Zboson_pT[s].Fill(Z.Pt(), weight)

                if s == 0:
                    if len(muon_GEN_match) < 2:
                        continue
                    if GEN_id[muon_GEN_match[0]] + GEN_id[muon_GEN_match[1]] != 0:
                        continue

                    lep0 = ROOT.TLorentzVector()
                    lep1 = ROOT.TLorentzVector()
                    lep0.SetPtEtaPhiM(GEN_pt[muon_GEN_match[0]], GEN_eta[muon_GEN_match[0]], GEN_phi[muon_GEN_match[0]], GEN_mass[muon_GEN_match[0]])
                    lep1.SetPtEtaPhiM(GEN_pt[muon_GEN_match[1]], GEN_eta[muon_GEN_match[1]], GEN_phi[muon_GEN_match[1]], GEN_mass[muon_GEN_match[1]])
                    Z = lep0 + lep1

                    if muon_charge[muon_hlt_match[0]] * GEN_id[muon_GEN_match[0]] < 0:
                        h_pTRes.Fill(muon_pt[muon_hlt_match[0]], (GEN_pt[muon_GEN_match[0]] - muon_pt[muon_hlt_match[0]]) / GEN_pt[muon_GEN_match[0]])
                        h_pTRes.Fill(muon_pt[muon_hlt_match[1]], (GEN_pt[muon_GEN_match[1]] - muon_pt[muon_hlt_match[1]]) / GEN_pt[muon_GEN_match[1]])
                    else:
                        h_pTRes.Fill(muon_pt[muon_hlt_match[0]], (GEN_pt[muon_GEN_match[1]] - muon_pt[muon_hlt_match[0]]) / GEN_pt[muon_GEN_match[1]])
                        h_pTRes.Fill(muon_pt[muon_hlt_match[1]], (GEN_pt[muon_GEN_match[0]] - muon_pt[muon_hlt_match[1]]) / GEN_pt[muon_GEN_match[0]])

    
    save = USER_DIR + "rochester_v1.pdf"

    h_resGEN = ROOT.TH1F("h_resGEN", "h_resGEN", pt_bins.size()-1, pt_bins[0])
        

    for xxx in range(1, h_pTRes.GetNbinsX() + 1):
        name = "Res_pt_%.0f_%.0f" % (pt_bins[xxx - 1], pt_bins[xxx])

        h_tmp_MC = h_pTRes.ProjectionY(name, xxx, xxx)
        c1 = ROOT.TCanvas(name, name, 600, 600)
        
        var = ROOT.RooRealVar("var", "var", -0.20, 0.20)
        MeanA = ROOT.RooRealVar("MeanA", "MeanA", 0, -0.01, 0.01)
        Sigma_GaussA = ROOT.RooRealVar("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.0001, 0.5)
        GaussA = ROOT.RooGaussian("GaussA", "GaussA", var, MeanA, Sigma_GaussA)
        histo = ROOT.RooDataHist("prova_MC", "prova_MC", var, h_tmp_MC)
        xframe = var.frame(ROOT.RooFit.Title(name))
        histo.plotOn(xframe)
        GaussA.fitTo(histo, ROOT.RooFit.Range(-0.03, 0.03))
        GaussA.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed + 2))
        GaussA.paramOn(xframe, ROOT.RooFit.Layout(0.1, 0.4, 0.7))
        xframe.Draw()

        h_resGEN.SetBinContent(xxx, Sigma_GaussA.getVal())

        if xxx == 1:
            c1.Print(save + "[")

        c1.Print(save)

    h_ZBoson_pt_correction = h_Zboson_pT[1].Clone()

    # DrawResolution(h_Zboson_pT[0], h_Zboson_pT[1], h_Zboson_pT[2], samples, "Zboson_pT", save, "pT_{2l} (GeV)")

    h_ZBoson_pt_correction.Divide(h_Zboson_pT[0])
    c_ZBoson_pt_correction = ROOT.TCanvas("ZBoson_pt_correction", "ZBoson_pt_correction", 600, 600)
    h_ZBoson_pt_correction.Draw()
    c_ZBoson_pt_correction.Print(save)

    ZbosonPt = 1

    for s in range(3):
        h_Zboson_pT[s].Reset()

    for s in range(2):
        print(USER_DIR + samples[s])
        if Anno_2018:
            _file0 = ROOT.TFile(USER_DIR + samples[s] + ".root")

        if _file0:
            tree = _file0.Get("Ana/passedEvents")
        else:
            print("File not found")
            return

        if not tree:
            print("tree not found")
            return

        h_num_eventi = ROOT.TH1F(_file0.Get("Ana/sumWeights"))
        num_event = h_num_eventi.Integral()
        print("Num event =", num_event)

        tree.SetBranchAddress("muon_pt", muon_pt)
        tree.SetBranchAddress("muon_eta", muon_eta)
        tree.SetBranchAddress("muon_phi", muon_phi)
        tree.SetBranchAddress("muon_mass", muon_mass)
        tree.SetBranchAddress("muon_charge", muon_charge)
        tree.SetBranchAddress("muon_isLoose", muon_isLoose)
        tree.SetBranchAddress("muon_isMedium", muon_isMedium)
        tree.SetBranchAddress("muon_isTight", muon_isTight)
        tree.SetBranchAddress("muon_trackIso", muon_trackIso)

        tree.SetBranchAddress("GEN_pt", GEN_pt)
        tree.SetBranchAddress("GEN_eta", GEN_eta)
        tree.SetBranchAddress("GEN_phi", GEN_phi)
        tree.SetBranchAddress("GEN_mass", GEN_mass)
        tree.SetBranchAddress("GEN_id", GEN_id)
        tree.SetBranchAddress("GEN_status", GEN_status)
        tree.SetBranchAddress("GEN_motherId", GEN_motherId)

        tree.SetBranchAddress("hlt_pt", hlt_pt)
        tree.SetBranchAddress("hlt_eta", hlt_eta)
        tree.SetBranchAddress("hlt_phi", hlt_phi)

        tree.SetBranchAddress("passedTrig", passedTrig)

        nentries = tree.GetEntries()
        print(nentries)

        for entry in range(nentries):
            # for entry in range(2500000):
            # for entry in range(1500000):
            # for entry in range(500000):
            # for entry in range(100000):

            tree.GetEntry(entry)
            if entry == 0:
                print(weight)

            if entry % 250000 == 0:
                print(entry, "/", nentries)

            if len(muon_pt) < 2:
                continue

            muon_hlt_match = []
            muon_GEN_match = []

            for j in range(len(muon_pt)):
                if len(muon_hlt_match) > 1:
                    continue

                deltaR = 1000
                if muon_pt[j] < 10:
                    continue
                if abs(muon_eta[j]) > 2.4:
                    continue
                if muon_trackIso[j] / muon_pt[j] > 0.1:
                    continue
                if not muon_isMedium:
                    continue
                # if not muon_isMedium:
                #     continue
                # if not muon_isLoose:
                #     continue

                for i in range(len(hlt_pt)):
                    deltaEta = hlt_eta[i] - muon_eta[j]
                    deltaPhi = hlt_phi[i] - muon_phi[j]
                    tmp_deltaR = (deltaEta * deltaEta + deltaPhi * deltaPhi) ** 0.5
                    if tmp_deltaR < deltaR:
                        deltaR = tmp_deltaR

                if deltaR < 0.35:
                    muon_hlt_match.append(j)

                if s == 0:
                    deltaR = 1000
                    k = -1
                    for i in range(len(GEN_pt)):
                        if GEN_status[i] != 1:
                            continue
                        if abs(GEN_id[i]) != 13:
                            continue
                        if GEN_motherId[i] != 23:
                            continue

                        deltaEta = GEN_eta[i] - muon_eta[j]
                        deltaPhi = GEN_phi[i] - muon_phi[j]
                        tmp_deltaR = (deltaEta * deltaEta + deltaPhi * deltaPhi) ** 0.5
                        if tmp_deltaR < deltaR:
                            deltaR = tmp_deltaR
                            if deltaR < 0.1:
                                k = i

                    if k != -1:
                        muon_GEN_match.append(k)

            if len(muon_hlt_match) < 2:
                continue
            if muon_charge[muon_hlt_match[0]] + muon_charge[muon_hlt_match[1]] != 0:
                continue
            if muon_pt[muon_hlt_match[0]] < 25:
                continue

            lep0 = ROOT.TLorentzVector()
            lep1 = ROOT.TLorentzVector()
            lep0.SetPtEtaPhiM(muon_pt[muon_hlt_match[0]], muon_eta[muon_hlt_match[0]], muon_phi[muon_hlt_match[0]], muon_mass[muon_hlt_match[0]])
            lep1.SetPtEtaPhiM(muon_pt[muon_hlt_match[1]], muon_eta[muon_hlt_match[1]], muon_phi[muon_hlt_match[1]], muon_mass[muon_hlt_match[1]])
            Z = lep0 + lep1

            if 81 < Z.M() < 101:
                # if 60 < Z.M() < 120:

                if s == 0:
                    ZbosonPt = h_ZBoson_pt_correction.GetBinContent(h_ZBoson_pt_correction.FindBin(Z.Pt()))
                    weight = lumi * crossSection / num_event

                    reweight = ZbosonPt * weight
                else:
                    weight = 1
                    reweight = 1

                h_Zboson_mass[s].Fill(Z.M(), reweight)
                h_Zboson_pT[s].Fill(Z.Pt(), reweight)

                h_lepton_qOverPt[s].Fill(muon_charge[muon_hlt_match[0]] / muon_pt[muon_hlt_match[0]], reweight)
                h_lepton_qOverPt[s].Fill(muon_charge[muon_hlt_match[1]] / muon_pt[muon_hlt_match[1]], reweight)

                h_lepton_pt[s].Fill(muon_pt[muon_hlt_match[0]], reweight)
                h_lepton_pt[s].Fill(muon_pt[muon_hlt_match[1]], reweight)

                h_lepton_eta[s].Fill(muon_eta[muon_hlt_match[0]], reweight)
                h_lepton_eta[s].Fill(muon_eta[muon_hlt_match[1]], reweight)

                h_lepton_phi[s].Fill(muon_phi[muon_hlt_match[0]], reweight)
                h_lepton_phi[s].Fill(muon_phi[muon_hlt_match[1]], reweight)

                if muon_charge[muon_hlt_match[0]] > 0:
                    h_lepton_qOverPt_scan_EtaPhi[s][0].Fill(muon_eta[muon_hlt_match[0]], muon_phi[muon_hlt_match[0]], 1 / muon_pt[muon_hlt_match[0]], reweight)
                    h_lepton_qOverPt_scan_EtaPhi[s][1].Fill(muon_eta[muon_hlt_match[1]], muon_phi[muon_hlt_match[1]], 1 / muon_pt[muon_hlt_match[1]], reweight)

                    h_Zboson_mass_vs_eta[s][0].Fill(muon_eta[muon_hlt_match[0]], Z.M(), reweight)
                    h_Zboson_mass_vs_eta[s][1].Fill(muon_eta[muon_hlt_match[1]], Z.M(), reweight)
                    h_Zboson_mass_vs_phi[s][0].Fill(muon_phi[muon_hlt_match[0]], Z.M(), reweight)
                    h_Zboson_mass_vs_phi[s][1].Fill(muon_phi[muon_hlt_match[1]], Z.M(), reweight)

                else:
                    h_lepton_qOverPt_scan_EtaPhi[s][0].Fill(muon_eta[muon_hlt_match[1]], muon_phi[muon_hlt_match[1]], 1 / muon_pt[muon_hlt_match[1]], reweight)
                    h_lepton_qOverPt_scan_EtaPhi[s][1].Fill(muon_eta[muon_hlt_match[0]], muon_phi[muon_hlt_match[0]], 1 / muon_pt[muon_hlt_match[0]], reweight)

                    h_Zboson_mass_vs_eta[s][0].Fill(muon_eta[muon_hlt_match[1]], Z.M(), reweight)
                    h_Zboson_mass_vs_eta[s][1].Fill(muon_eta[muon_hlt_match[0]], Z.M(), reweight)
                    h_Zboson_mass_vs_phi[s][0].Fill(muon_phi[muon_hlt_match[1]], Z.M(), reweight)
                    h_Zboson_mass_vs_phi[s][1].Fill(muon_phi[muon_hlt_match[0]], Z.M(), reweight)

    DrawResolution(h_Zboson_mass[0], h_Zboson_mass[1], h_Zboson_mass[2], samples, "Zboson_mass", save, "mass_{2l} (GeV)")
    DrawResolution(h_Zboson_pT[0], h_Zboson_pT[1], h_Zboson_pT[2], samples, "Zboson_pT", save, "pT_{2l} (GeV)")
    DrawResolution(h_lepton_pt[0], h_lepton_pt[1], h_lepton_pt[2], samples, "h_lepton_pt", save, "p_{T} (GeV)")
    DrawResolution(h_lepton_eta[0], h_lepton_eta[1], h_lepton_eta[2], samples, "h_lepton_eta", save, "#eta")
    DrawResolution(h_lepton_phi[0], h_lepton_phi[1], h_lepton_phi[2], samples, "h_lepton_phi", save, "#phi")
    DrawResolution(h_lepton_qOverPt[0], h_lepton_qOverPt[1], h_lepton_qOverPt[2], samples, "h_lepton_qOverPt", save, "q/p_{T} (GeV)")

    # Create TMultiGraph objects
    Mass_vs_phi_pos = ROOT.TMultiGraph()
    Mass_vs_phi_neg = ROOT.TMultiGraph()
    Mass_vs_eta_pos = ROOT.TMultiGraph()
    Mass_vs_eta_neg = ROOT.TMultiGraph()

    for s in range(3):
        mean_pos = []
        mean_neg = []
        eta_bins_centered = []
        phi_bins_centered = []

        for xxx in range(1, h_Zboson_mass_vs_eta[s][0].GetNbinsX() + 1):
            name = "Mass_eta_{:.1f}_{:.1f}_pos".format(eta_bins[xxx - 1], eta_bins[xxx])
            
            h_tmp_MC = ROOT.TH1F(h_Zboson_mass_vs_eta[s][0].ProjectionY(name, xxx, xxx))
            massZ_P = ROOT.RooRealVar("massZ", "massZ", 81, 101)
            histo = ROOT.RooDataHist("mass Z", "mass Z", massZ_P, h_tmp_MC)
            
            PDGmZ = ROOT.RooRealVar("PDGmZ", "PDGmZ", 91.19, 86, 96)
            PDGwZ = ROOT.RooRealVar("PDGwZ", "PDGwZ", 2.5, 1, 4)
            PDGBW = ROOT.RooBreitWigner("PDGBW", "PDGBW", massZ_P, PDGmZ, PDGwZ)
            PDGBW_2 = ROOT.RooBreitWigner("PDGBW_2", "PDGBW_2", massZ_P, PDGmZ, PDGwZ)

            CB_mean = ROOT.RooRealVar("CB_mean", "CB_mean", 0, -5, 5)
            Sigma = ROOT.RooRealVar("Sigma", "Sigma", 1, 0.1, 20)
            AlphaL = ROOT.RooRealVar("AlphaL", "AlphaL", 1, 0.1, 30)
            ExpL = ROOT.RooRealVar("ExpL", "ExpL", 1, 0.1, 30)
            AlphaR = ROOT.RooRealVar("AlphaR", "AlphaR", 1, 0.1, 30)
            ExpR = ROOT.RooRealVar("ExpR", "ExpR", 1, 0.1, 50)
            DSCB = ROOT.RooDoubleCB("DSCB", "DSCB", massZ_P, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR)

            model = ROOT.RooFFTConvPdf("CW", "CW", massZ_P, PDGBW, DSCB)
            model_2 = ROOT.RooFFTConvPdf("CW_2", "CW_2", massZ_P, PDGBW, DSCB)

            if s != 20:
                model.fitTo(histo, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))
            else:
                PDGBW.fitTo(histo, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))

            mean_pos.append(PDGmZ.getVal())
            eta_bins_centered.append(eta_bins[xxx - 1] + abs((eta_bins[1] - eta_bins[0]) / 2))

            c_MC = ROOT.TCanvas(name, name, 900, 700)
            c_MC.SetFrameFillColor(0)
            xframe = massZ_P.frame("massZ")
            xframe = massZ_P.frame(ROOT.Title(name))
            histo.plotOn(xframe)
            if s != 20:
                model.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kBlue))
                model.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            else:
                PDGBW.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))
                PDGBW.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            xframe.Draw()
            c_MC.Print(save + name)

            name = "Mass_eta_{:.1f}_{:.1f}_neg".format(eta_bins[xxx - 1], eta_bins[xxx])
            h_tmp_MC_2 = ROOT.TH1F(h_Zboson_mass_vs_eta[s][1].ProjectionY(name, xxx, xxx))
            histo_neg = ROOT.RooDataHist("mass Z", "mass Z", massZ_P, h_tmp_MC_2)
            if s != 20:
                model_2.fitTo(histo_neg, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))
            else:
                PDGBW_2.fitTo(histo_neg, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))
            mean_neg.append(PDGmZ.getVal())

            c_MC = ROOT.TCanvas(name, name, 900, 700)
            c_MC.SetFrameFillColor(0)
            xframe = massZ_P.frame("massZ")
            xframe = massZ_P.frame(ROOT.Title(name))
            histo_neg.plotOn(xframe)
            if s != 20:
                model_2.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kBlue))
                model_2.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            else:
                PDGBW_2.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))
                PDGBW_2.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            xframe.Draw()
            c_MC.Print(save + name)

        gr_new = ROOT.TGraph(len(eta_bins) - 1, array("d", eta_bins_centered), array("d", mean_pos))
        if s == 1:
            gr_new.SetMarkerColor(ROOT.kRed)
        if s == 2:
            gr_new.SetMarkerColor(ROOT.kGreen)
        Mass_vs_eta_pos.Add(gr_new)

        gr_new = ROOT.TGraph(len(eta_bins) - 1, array("d", eta_bins_centered), array("d", mean_neg))
        if s == 1:
            gr_new.SetMarkerColor(ROOT.kRed)
        if s == 2:
            gr_new.SetMarkerColor(ROOT.kGreen)
        Mass_vs_eta_neg.Add(gr_new)

        for xxx in range(1, h_Zboson_mass_vs_phi[s][0].GetNbinsX() + 1):
            name = "Mass_phi_{:.1f}_{:.1f}_pos".format(phi_bins[xxx - 1], phi_bins[xxx])
            
            h_tmp_MC = ROOT.TH1F(h_Zboson_mass_vs_phi[s][0].ProjectionY(name, xxx, xxx))
            massZ_P = ROOT.RooRealVar("massZ", "massZ", 81, 101)
            histo = ROOT.RooDataHist("mass Z", "mass Z", massZ_P, h_tmp_MC)
            
            PDGmZ = ROOT.RooRealVar("PDGmZ", "PDGmZ", 91.19, 86, 96)
            PDGwZ = ROOT.RooRealVar("PDGwZ", "PDGwZ", 2.5, 1, 4)
            PDGBW = ROOT.RooBreitWigner("PDGBW", "PDGBW", massZ_P, PDGmZ, PDGwZ)
            PDGBW_2 = ROOT.RooBreitWigner("PDGBW_2", "PDGBW_2", massZ_P, PDGmZ, PDGwZ)

            CB_mean = ROOT.RooRealVar("CB_mean", "CB_mean", 0, -5, 5)
            Sigma = ROOT.RooRealVar("Sigma", "Sigma", 1, 0.1, 20)
            AlphaL = ROOT.RooRealVar("AlphaL", "AlphaL", 1, 0.1, 30)
            ExpL = ROOT.RooRealVar("ExpL", "ExpL", 1, 0.1, 30)
            AlphaR = ROOT.RooRealVar("AlphaR", "AlphaR", 1, 0.1, 30)
            ExpR = ROOT.RooRealVar("ExpR", "ExpR", 1, 0.1, 50)
            DSCB = ROOT.RooDoubleCB("DSCB", "DSCB", massZ_P, CB_mean, Sigma, AlphaL, ExpL, AlphaR, ExpR)

            model = ROOT.RooFFTConvPdf("CW", "CW", massZ_P, PDGBW, DSCB)
            model_2 = ROOT.RooFFTConvPdf("CW_2", "CW_2", massZ_P, PDGBW, DSCB)

            if s != 20:
                model.fitTo(histo, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))
            else:
                PDGBW.fitTo(histo, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))

            mean_pos.append(PDGmZ.getVal())
            phi_bins_centered.append(phi_bins[xxx - 1] + abs((phi_bins[1] - phi_bins[0]) / 2))

            c_MC = ROOT.TCanvas(name, name, 900, 700)
            c_MC.SetFrameFillColor(0)
            xframe = massZ_P.frame("massZ")
            xframe = massZ_P.frame(ROOT.Title(name))
            histo.plotOn(xframe)
            if s != 20:
                model.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kBlue))
                model.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            else:
                PDGBW.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))
                PDGBW.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            xframe.Draw()
            c_MC.Print(save + name)

            name = "Mass_phi_{:.1f}_{:.1f}_neg".format(phi_bins[xxx - 1], phi_bins[xxx])
            h_tmp_MC_2 = ROOT.TH1F(h_Zboson_mass_vs_phi[s][1].ProjectionY(name, xxx, xxx))
            histo_neg = ROOT.RooDataHist("mass Z", "mass Z", massZ_P, h_tmp_MC_2)
            if s != 20:
                model_2.fitTo(histo_neg, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))
            else:
                PDGBW_2.fitTo(histo_neg, ROOT.Range(86, 96), ROOT.Save(True), ROOT.SumW2Error(True), ROOT.Verbose(False),
                            ROOT.PrintLevel(-1), ROOT.Warnings(False), ROOT.NumCPU(12), ROOT.Timer(True))
            mean_neg.append(PDGmZ.getVal())

            c_MC = ROOT.TCanvas(name, name, 900, 700)
            c_MC.SetFrameFillColor(0)
            xframe = massZ_P.frame("massZ")
            xframe = massZ_P.frame(ROOT.Title(name))
            histo_neg.plotOn(xframe)
            if s != 20:
                model_2.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kBlue))
                model_2.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            else:
                PDGBW_2.plotOn(xframe, ROOT.RooFit.LineColor(ROOT.kRed))
                PDGBW_2.paramOn(xframe, ROOT.RooFit.Layout(0.13, 0.5, 0.80))
            xframe.Draw()
            c_MC.Print(save + name)

        gr_new = ROOT.TGraph(len(phi_bins) - 1, array("d", phi_bins_centered), array("d", mean_pos))
        if s == 1:
            gr_new.SetMarkerColor(ROOT.kRed)
        if s == 2:
            gr_new.SetMarkerColor(ROOT.kGreen)
        Mass_vs_phi_pos.Add(gr_new)

        gr_new = ROOT.TGraph(len(phi_bins) - 1, array("d", phi_bins_centered), array("d", mean_neg))
        if s == 1:
            gr_new.SetMarkerColor(ROOT.kRed)
        if s == 2:
            gr_new.SetMarkerColor(ROOT.kGreen)
        Mass_vs_phi_neg.Add(gr_new)
    
    for xxx in range(1, h_Zboson_mass_vs_eta[0][0].GetNbinsX() + 1):
        name = "Mass_eta_{:.1f}_{:.1f}_pos_MC".format(eta_bins[xxx - 1], eta_bins[xxx])
        h_tmp_MC = ROOT.TH1F(h_Zboson_mass_vs_eta[0][0].ProjectionY(name, xxx, xxx))
        name = "Mass_eta_{:.1f}_{:.1f}_pos_DATA".format(eta_bins[xxx - 1], eta_bins[xxx])
        h_tmp_DATA = ROOT.TH1F(h_Zboson_mass_vs_eta[1][0].ProjectionY(name, xxx, xxx))
        name = "Mass_eta_{:.1f}_{:.1f}_pos".format(eta_bins[xxx - 1], eta_bins[xxx])
        h_tmp_GEN = ROOT.TH1F(h_Zboson_mass_vs_eta[2][0].ProjectionY(name, xxx, xxx))

        DrawResolution(h_tmp_MC, h_tmp_DATA, h_tmp_GEN, samples, name, save, "mass_{2l} (GeV)")

        name = "Mass_eta_{:.1f}_{:.1f}_neg_MC".format(eta_bins[xxx - 1], eta_bins[xxx])
        h_tmp_MC_2 = ROOT.TH1F(h_Zboson_mass_vs_eta[0][1].ProjectionY(name, xxx, xxx))
        name = "Mass_eta_{:.1f}_{:.1f}_neg_DATA".format(eta_bins[xxx - 1], eta_bins[xxx])
        h_tmp_DATA_2 = ROOT.TH1F(h_Zboson_mass_vs_eta[1][1].ProjectionY(name, xxx, xxx))
        name = "Mass_eta_{:.1f}_{:.1f}_neg".format(eta_bins[xxx - 1], eta_bins[xxx])
        h_tmp_GEN_2 = ROOT.TH1F(h_Zboson_mass_vs_eta[2][1].ProjectionY(name, xxx, xxx))
        name = "Mass_eta_{:.1f}_{:.1f}_neg".format(eta_bins[xxx - 1], eta_bins[xxx])

        DrawResolution(h_tmp_MC_2, h_tmp_DATA_2, h_tmp_GEN_2, samples, name, save, "mass_{2l} (GeV)")

    for xxx in range(1, h_Zboson_mass_vs_phi[0][0].GetNbinsX() + 1):
        name = "Mass_phi_{:.1f}_{:.1f}_pos_MC".format(phi_bins[xxx - 1], phi_bins[xxx])
        h_tmp_MC = ROOT.TH1F(h_Zboson_mass_vs_phi[0][0].ProjectionY(name, xxx, xxx))
        name = "Mass_phi_{:.1f}_{:.1f}_pos_DATA".format(phi_bins[xxx - 1], phi_bins[xxx])
        h_tmp_DATA = ROOT.TH1F(h_Zboson_mass_vs_phi[1][0].ProjectionY(name, xxx, xxx))
        name = "Mass_phi_{:.1f}_{:.1f}_pos".format(phi_bins[xxx - 1], phi_bins[xxx])
        h_tmp_GEN = ROOT.TH1F(h_Zboson_mass_vs_phi[2][0].ProjectionY(name, xxx, xxx))

        DrawResolution(h_tmp_MC, h_tmp_DATA, h_tmp_GEN, samples, name, save, "mass_{2l} (GeV)")

        name = "Mass_phi_{:.1f}_{:.1f}_neg_MC".format(phi_bins[xxx - 1], phi_bins[xxx])
        h_tmp_MC_2 = ROOT.TH1F(h_Zboson_mass_vs_phi[0][1].ProjectionY(name, xxx, xxx))
        name = "Mass_phi_{:.1f}_{:.1f}_neg_DATA".format(phi_bins[xxx - 1], phi_bins[xxx])
        h_tmp_DATA_2 = ROOT.TH1F(h_Zboson_mass_vs_phi[1][1].ProjectionY(name, xxx, xxx))
        name = "Mass_phi_{:.1f}_{:.1f}_neg".format(phi_bins[xxx - 1], phi_bins[xxx])
        h_tmp_GEN_2 = ROOT.TH1F(h_Zboson_mass_vs_phi[2][1].ProjectionY(name, xxx, xxx))
        name = "Mass_phi_{:.1f}_{:.1f}_neg".format(phi_bins[xxx - 1], phi_bins[xxx])

        DrawResolution(h_tmp_MC_2, h_tmp_DATA_2, h_tmp_GEN_2, samples, name, save, "mass_{2l} (GeV)")
    
    c_Mass_vs_eta_pos = ROOT.TCanvas("c_Mass_vs_eta_pos", "c_Mass_vs_eta_pos", 600, 600)
    Mass_vs_eta_pos.SetTitle("c_Mass_vs_eta_pos")
    Mass_vs_eta_pos.Draw("AP*")
    c_Mass_vs_eta_pos.Print(save)

    c_Mass_vs_eta_neg = ROOT.TCanvas("c_Mass_vs_eta_neg", "c_Mass_vs_eta_neg", 600, 600)
    Mass_vs_eta_neg.SetTitle("c_Mass_vs_eta_neg")
    Mass_vs_eta_neg.Draw("AP*")
    c_Mass_vs_eta_neg.Print(save)

    c_Mass_vs_phi_pos = ROOT.TCanvas("c_Mass_vs_phi_pos", "c_Mass_vs_phi_pos", 600, 600)
    Mass_vs_phi_pos.SetTitle("c_Mass_vs_phi_pos")
    Mass_vs_phi_pos.Draw("AP*")
    c_Mass_vs_phi_pos.Print(save)

    c_Mass_vs_phi_neg = ROOT.TCanvas("c_Mass_vs_phi_neg", "c_Mass_vs_phi_neg", 600, 600)
    Mass_vs_phi_neg.SetTitle("c_Mass_vs_phi_neg")
    Mass_vs_phi_neg.Draw("AP*")
    c_Mass_vs_phi_neg.Print(save)
    # c_Mass_vs_phi_neg.Print(save + "]")

    for s in range(3):
        h_Zboson_mass[s].Reset()
        h_Zboson_pT[s].Reset()
        h_lepton_pt[s].Reset()
        h_lepton_eta[s].Reset()
        h_lepton_phi[s].Reset()
        h_lepton_qOverPt[s].Reset()
        h_Zboson_mass_vs_phi[s][0].Reset()
        h_Zboson_mass_vs_phi[s][1].Reset()
        h_Zboson_mass_vs_eta[s][0].Reset()
        h_Zboson_mass_vs_eta[s][1].Reset()

    Corr_DATA_pos = 0.0
    Corr_MC_pos = 0.0
    Corr_DATA_neg = 0.0
    Corr_MC_neg = 0.0
    D_m_MC = 0.0
    D_a_MC = 0.0
    D_m_DATA = 0.0
    D_a_DATA = 0.0

    h_Multiplicative_MC = ROOT.TH2F("Multiplicative_MC", "Multiplicative_MC", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_Additive_MC = ROOT.TH2F("Additive_MC", "Additive_MC", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_Multiplicative_DATA = ROOT.TH2F("Multiplicative_DATA", "Multiplicative_DATA", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_Additive_DATA = ROOT.TH2F("Additive_DATA", "Additive_DATA", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_denominator_MC_pos = ROOT.TH2F("h_denominator_MC_pos", "h_denominator_MC_pos", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_denominator_MC_neg = ROOT.TH2F("h_denominator_MC_neg", "h_denominator_MC_neg", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_denominator_DATA_pos = ROOT.TH2F("h_denominator_DATA_pos", "h_denominator_DATA_pos", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))
    h_denominator_DATA_neg = ROOT.TH2F("h_denominator_DATA_neg", "h_denominator_DATA_neg", len(eta_bins) - 1, array('d', eta_bins), len(phi_bins) - 1, array('d', phi_bins))

    for xxx in range(1, h_lepton_qOverPt_scan_EtaPhi[0][0].GetNbinsX() + 1):
        for yyy in range(1, h_lepton_qOverPt_scan_EtaPhi[0][0].GetNbinsY() + 1):
            name = "1overPt_eta_{:.1f}_{:.1f}_phi_{:.1f}_{:.1f}".format(
                eta_bins[xxx - 1], eta_bins[xxx], phi_bins[yyy - 1], phi_bins[yyy]
            )
            print("========")
            print(name)

            h_tmp_MC = h_lepton_qOverPt_scan_EtaPhi[0][0].ProjectionZ(
                name + "_MC+", xxx, xxx, yyy, yyy
            )
            h_MC = ROOT.TH1F(name + "_MC", name + "_MC", 100, -0.1, 0.1)
            for binx in range(1, h_tmp_MC.GetNbinsX() + 1):
                h_MC.SetBinContent(binx, h_tmp_MC.GetBinContent(binx))
            h_tmp_DATA = h_lepton_qOverPt_scan_EtaPhi[1][0].ProjectionZ(
                name + "_DATA+", xxx, xxx, yyy, yyy
            )
            h_DATA = ROOT.TH1F(name + "_DATA", name + "_DATA", 100, -0.1, 0.1)
            for binx in range(1, h_tmp_DATA.GetNbinsX() + 1):
                h_DATA.SetBinContent(binx, h_tmp_DATA.GetBinContent(binx))
            h_tmp_GEN = h_lepton_qOverPt_scan_EtaPhi[2][0].ProjectionZ(
                name + "_GEN+", xxx, xxx, yyy, yyy
            )
            h_GEN = ROOT.TH1F(name + "_GEN", name + "_GEN", 100, -0.1, 0.1)
            for binx in range(1, h_tmp_GEN.GetNbinsX() + 1):
                h_GEN.SetBinContent(binx, h_tmp_GEN.GetBinContent(binx))

            h_denominator_MC_pos.SetBinContent(xxx, yyy, h_MC.GetMean())
            h_denominator_DATA_pos.SetBinContent(xxx, yyy, h_DATA.GetMean())

            print("Before\nDenominator+ =", h_MC.GetMean(), h_DATA.GetMean())

            Corr_MC_pos = h_GEN.GetMean() - h_MC.GetMean()
            Corr_DATA_pos = h_GEN.GetMean() - h_DATA.GetMean()
            print("Mean+  =", h_GEN.GetMean(), h_MC.GetMean(), h_DATA.GetMean())

            h_MC.Reset()
            h_DATA.Reset()
            h_GEN.Reset()

            h_tmp_MC = h_lepton_qOverPt_scan_EtaPhi[0][1].ProjectionZ(
                name + "_MC-", xxx, xxx, yyy, yyy
            )
            h_MC = ROOT.TH1F(name + "_MC", name + "_MC", 100, -0.1, 0.1)
            for binx in range(1, h_tmp_MC.GetNbinsX() + 1):
                h_MC.SetBinContent(binx, h_tmp_MC.GetBinContent(binx))
            h_tmp_DATA = h_lepton_qOverPt_scan_EtaPhi[1][1].ProjectionZ(
                name + "_DATA-", xxx, xxx, yyy, yyy
            )
            h_DATA = ROOT.TH1F(name + "_DATA", name + "_DATA", 100, -0.1, 0.1)
            for binx in range(1, h_tmp_DATA.GetNbinsX() + 1):
                h_DATA.SetBinContent(binx, h_tmp_DATA.GetBinContent(binx))
            h_tmp_GEN = h_lepton_qOverPt_scan_EtaPhi[2][1].ProjectionZ(
                name + "_GEN-", xxx, xxx, yyy, yyy
            )
            h_GEN = ROOT.TH1F(name + "_GEN", name + "_GEN", 100, -0.1, 0.1)
            for binx in range(1, h_tmp_GEN.GetNbinsX() + 1):
                h_GEN.SetBinContent(binx, h_tmp_GEN.GetBinContent(binx))

            h_denominator_MC_neg.SetBinContent(xxx, yyy, h_MC.GetMean())
            h_denominator_DATA_neg.SetBinContent(xxx, yyy, h_DATA.GetMean())

            print("Denominator- =", h_MC.GetMean(), h_DATA.GetMean())

            Corr_MC_neg = h_GEN.GetMean() - h_MC.GetMean()
            Corr_DATA_neg = h_GEN.GetMean() - h_DATA.GetMean()
            print("Mean-  =", h_GEN.GetMean(), h_MC.GetMean(), h_DATA.GetMean())

            print("C_MC =", Corr_MC_pos, Corr_MC_neg)
            print("C_DATA =", Corr_DATA_pos, Corr_DATA_neg)

            D_m_MC = (Corr_MC_pos + Corr_MC_neg) / 2.0
            D_a_MC = (Corr_MC_pos - Corr_MC_neg) / 2.0
            D_m_DATA = (Corr_DATA_pos + Corr_DATA_neg) / 2.0
            D_a_DATA = (Corr_DATA_pos - Corr_DATA_neg) / 2.0
            print("D_m_MC =", D_m_MC)
            print("D_m_DATA =", D_m_DATA)
            print("D_a_MC =", D_a_MC)
            print("D_a_DATA =", D_a_DATA)

            M = D_m_MC
            h_Additive_MC.SetBinContent(xxx, yyy, D_a_MC)
            h_Multiplicative_MC.SetBinContent(xxx, yyy, M)
            print("MC: M =", M, "A =", D_a_MC)

            M = D_m_DATA
            h_Additive_DATA.SetBinContent(xxx, yyy, D_a_DATA)
            h_Multiplicative_DATA.SetBinContent(xxx, yyy, M)
            print("DATA: M =", M, "A =", D_a_DATA)

    c_Additive_MC = ROOT.TCanvas("c_Additive_MC", "c_Additive_MC", 600, 600)
    c_Additive_MC.SetRightMargin(0.2)
    h_Additive_MC.Draw("COLZ")
    c_Additive_MC.Print(save)
    c_Additive_MC.Close()

    c_Additive_DATA = ROOT.TCanvas("c_Additive_DATA", "c_Additive_DATA", 600, 600)
    c_Additive_DATA.SetRightMargin(0.2)
    h_Additive_DATA.Draw("COLZ")
    c_Additive_DATA.Print(save)
    c_Additive_DATA.Close()

    c_Multiplicative_MC = ROOT.TCanvas("c_Multiplicative_MC", "c_Multiplicative_MC", 600, 600)
    c_Multiplicative_MC.SetRightMargin(0.2)
    h_Multiplicative_MC.Draw("COLZ")
    c_Multiplicative_MC.Print(save)
    c_Multiplicative_MC.Close()

    c_Multiplicative_DATA = ROOT.TCanvas("c_Multiplicative_DATA", "c_Multiplicative_DATA", 600, 600)
    c_Multiplicative_DATA.SetRightMargin(0.2)
    h_Multiplicative_DATA.Draw("COLZ")
    c_Multiplicative_DATA.Print(save)
    c_Multiplicative_DATA.Close()

    # Clear histograms
    for s in range(3):
        for q in range(2):
            h_lepton_qOverPt_scan_EtaPhi[s][q].Reset()

    # Create histograms
    h_MC_mass = ROOT.TH2F("h_MC_mass", "h_MC_mass", 40, 81, 101, 40, 81, 101)
    h_DATA_mass = ROOT.TH2F("h_DATA_mass", "h_DATA_mass", 40, 81, 101, 40, 81, 101)

    h_pTVariance_MC_pos = ROOT.TH1F("h_pTVariance_MC_pos", "h_pTVariance_MC_pos", 100, -0.025, 0.025)
    h_pTVariance_MC_neg = ROOT.TH1F("h_pTVariance_MC_neg", "h_pTVariance_MC_neg", 100, -0.025, 0.025)
    h_pTVariance_DATA_pos = ROOT.TH1F("h_pTVariance_DATA_pos", "h_pTVariance_DATA_pos", 100, -0.1, 0.1)
    h_pTVariance_DATA_neg = ROOT.TH1F("h_pTVariance_DATA_neg", "h_pTVariance_DATA_neg", 100, -0.1, 0.1)

    for s in range(2):
        if Anno_2018:
            _file0 = ROOT.TFile(USER_DIR + samples[s] + ".root")

        if _file0:
            tree = _file0.Get("Ana/passedEvents")
        else:
            print("File not found")
            return

        if not tree:
            print("Tree not found")
            return

        h_num_eventi = _file0.Get("Ana/sumWeights")
        num_event = h_num_eventi.Integral()

        tree.SetBranchAddress("muon_pt", ROOT.AddressOf(muon_pt))
        tree.SetBranchAddress("muon_eta", ROOT.AddressOf(muon_eta))
        tree.SetBranchAddress("muon_phi", ROOT.AddressOf(muon_phi))
        tree.SetBranchAddress("muon_mass", ROOT.AddressOf(muon_mass))
        tree.SetBranchAddress("muon_charge", ROOT.AddressOf(muon_charge))
        tree.SetBranchAddress("muon_isLoose", ROOT.AddressOf(muon_isLoose))
        tree.SetBranchAddress("muon_isMedium", ROOT.AddressOf(muon_isMedium))
        tree.SetBranchAddress("muon_isTight", ROOT.AddressOf(muon_isTight))
        tree.SetBranchAddress("muon_trackIso", ROOT.AddressOf(muon_trackIso))

        tree.SetBranchAddress("GEN_pt", ROOT.AddressOf(GEN_pt))
        tree.SetBranchAddress("GEN_eta", ROOT.AddressOf(GEN_eta))
        tree.SetBranchAddress("GEN_phi", ROOT.AddressOf(GEN_phi))
        tree.SetBranchAddress("GEN_mass", ROOT.AddressOf(GEN_mass))
        tree.SetBranchAddress("GEN_id", ROOT.AddressOf(GEN_id))
        tree.SetBranchAddress("GEN_status", ROOT.AddressOf(GEN_status))
        tree.SetBranchAddress("GEN_motherId", ROOT.AddressOf(GEN_motherId))

        tree.SetBranchAddress("hlt_pt", ROOT.AddressOf(hlt_pt))
        tree.SetBranchAddress("hlt_eta", ROOT.AddressOf(hlt_eta))
        tree.SetBranchAddress("hlt_phi", ROOT.AddressOf(hlt_phi))

        tree.SetBranchAddress("passedTrig", ROOT.AddressOf(passedTrig))

        nentries = tree.GetEntries()

        print(nentries)

        
#zeile 1090 im C++ File




