import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

# Erstelle ein leeres Histogramm
h_blank = ROOT.TH1F("h_blank", "h_blank", 10, -0.5, 9.5)


h_num_eventi = ROOT.TH1F()
num_event = float()
nentries = ROOT.Long64_t()

hMuScaleFac = ROOT.TH2F()
hMuScaleFacUnc = ROOT.TH2F()

h_mu_B = ROOT.TH1F()
h_mu_E = ROOT.TH1F()
h_el_B = ROOT.TH1F()
h_el_E = ROOT.TH1F()
x_pTaxis = ROOT.TAxis()

Run, Event, LumiSect = ROOT.ULong64_t(), ROOT.ULong64_t(), ROOT.ULong64_t()

pT1_FromMuonBestTrack = ROOT.Double(0)
pT2_FromMuonBestTrack = ROOT.Double(0)
eta1_FromMuonBestTrack = ROOT.Double(0)
eta2_FromMuonBestTrack = ROOT.Double(0)
phi1_FromMuonBestTrack = ROOT.Double(0)
phi2_FromMuonBestTrack = ROOT.Double(0)
muonPV_x1 = ROOT.Double(0)
muonPV_x2 = ROOT.Double(0)
muonPV_y1 = ROOT.Double(0)
muonPV_y2 = ROOT.Double(0)
muonPV_z1 = ROOT.Double(0)
muonPV_z2 = ROOT.Double(0)
m1 = ROOT.Double(0)
pT1 = ROOT.Double(0)
eta1 = ROOT.Double(0)
phi1 = ROOT.Double(0)
pterr1old = ROOT.Double(0)
Iso1 = ROOT.Double(0)
pT_FSR1 = ROOT.Double(0)
eta_FSR1 = ROOT.Double(0)
phi_FSR1 = ROOT.Double(0)
m_FSR1 = ROOT.Double(0)
Id1 = ROOT.Int(0)
d0BS1 = ROOT.Double(0)
d0BS2 = ROOT.Double(0)
d0PV1 = ROOT.Double(0)
d0PV2 = ROOT.Double(0)
dzPV1 = ROOT.Double(0)
dzPV2 = ROOT.Double(0)
Sip1 = ROOT.Double(0)
Sip2 = ROOT.Double(0)
genLep_pt1 = ROOT.Double(0)
genLep_pt2 = ROOT.Double(0)
pterr1 = ROOT.Double(0)
pterr1_single = ROOT.Double(0)
pterr1_VX = ROOT.Double(0)
pterr1_VX_BS = ROOT.Double(0)
genLep_eta1 = ROOT.Double(0)
genLep_eta2 = ROOT.Double(0)
genLep_phi1 = ROOT.Double(0)
genLep_phi2 = ROOT.Double(0)
m2 = ROOT.Double(0)
pT2 = ROOT.Double(0)
eta2 = ROOT.Double(0)
phi2 = ROOT.Double(0)
pterr2old = ROOT.Double(0)
Iso2 = ROOT.Double(0)
pT_FSR2 = ROOT.Double(0)
eta_FSR2 = ROOT.Double(0)
phi_FSR2 = ROOT.Double(0)
m_FSR2 = ROOT.Double(0)
Id2 = ROOT.Int(0)
pterr2 = ROOT.Double(0)
pterr2_single = ROOT.Double(0)
pterr2_VX = ROOT.Double(0)
pterr2_VX_BS = ROOT.Double(0)
weight = ROOT.Float(0)
k_ggZZ = ROOT.Float(0)
k_qqZZ_ewk = ROOT.Float(0)
k_qqZZ_qcd_M = ROOT.Float(0)
lep1_ecalDriven = ROOT.Int(0)
lep2_ecalDriven = ROOT.Int(0)
massZ = ROOT.Double(0)
massZErr = ROOT.Double(0)
massZ_FSR = ROOT.Double(0)
massZErr_FSR = ROOT.Double(0)
genzm = ROOT.Double(0)
GENmass2l = ROOT.Double(0)
GENmass4l = ROOT.Float(0)
nZXCRFailedLeptons = ROOT.Int(0)

genWeight = ROOT.Float(0)
pileupWeight = ROOT.Float(0)
prefiringWeight = ROOT.Float(0)
dataMCWeight = ROOT.Float(0)
pileupWeightUp = ROOT.Float(0)
pileupWeightDn = ROOT.Float(0)

massZ_single_BS = ROOT.Double(0)
massZ_single_BS_FSR = ROOT.Double(0)
massZ_single_BS_BS = ROOT.Double(0)
massZ_single_BS_BS_FSR = ROOT.Double(0)
massZ_vtx = ROOT.Double(0)
massZ_vtx_FSR = ROOT.Double(0)
massZ_vtx_BS = ROOT.Double(0)
massZ_vtx_BS_FSR = ROOT.Double(0)

massZErr_single_BS = ROOT.Double(0)
massZErr_single_BS_FSR = ROOT.Double(0)
massZErr_vtx = ROOT.Double(0)
massZErr_vtx_FSR = ROOT.Double(0)
massZErr_vtx_BS = ROOT.Double(0)
massZErr_vtx_BS_FSR = ROOT.Double(0)

single_BS_pT1 = ROOT.Double(0)
single_BS_pT2 = ROOT.Double(0)
single_BS_eta1 = ROOT.Double(0)
single_BS_eta2 = ROOT.Double(0)
single_BS_m1 = ROOT.Double(0)
single_BS_m2 = ROOT.Double(0)
single_BS_phi1 = ROOT.Double(0)
single_BS_phi2 = ROOT.Double(0)

vtx_pT1 = ROOT.Double(0)
vtx_pT2 = ROOT.Double(0)
vtx_eta1 = ROOT.Double(0)
vtx_eta2 = ROOT.Double(0)
vtx_m1 = ROOT.Double(0)
vtx_m2 = ROOT.Double(0)
vtx_phi1 = ROOT.Double(0)
vtx_phi2 = ROOT.Double(0)

vtx_BS_pT1 = ROOT.Double(0)
vtx_BS_pT2 = ROOT.Double(0)
vtx_BS_eta1 = ROOT.Double(0)
vtx_BS_eta2 = ROOT.Double(0)
vtx_BS_m1 = ROOT.Double(0)
vtx_BS_m2 = ROOT.Double(0)
vtx_BS_phi1 = ROOT.Double(0)
vtx_BS_phi2 = ROOT.Double(0)

single_BS_pT_FSR1 = ROOT.Double(0)
single_BS_pT_FSR2 = ROOT.Double(0)
single_BS_eta_FSR1 = ROOT.Double(0)
single_BS_eta_FSR2 = ROOT.Double(0)
single_BS_m_FSR1 = ROOT.Double(0)
single_BS_m_FSR2 = ROOT.Double(0)
single_BS_phi_FSR1 = ROOT.Double(0)
single_BS_phi_FSR2 = ROOT.Double(0)

vtx_pT_FSR1 = ROOT.Double(0)
vtx_pT_FSR2 = ROOT.Double(0)
vtx_eta_FSR1 = ROOT.Double(0)
vtx_eta_FSR2 = ROOT.Double(0)
vtx_m_FSR1 = ROOT.Double(0)
vtx_m_FSR2 = ROOT.Double(0)
vtx_phi_FSR1 = ROOT.Double(0)
vtx_phi_FSR2 = ROOT.Double(0)

vtx_BS_pT_FSR1 = ROOT.Double(0)
vtx_BS_pT_FSR2 = ROOT.Double(0)
vtx_BS_eta_FSR1 = ROOT.Double(0)
vtx_BS_eta_FSR2 = ROOT.Double(0)
vtx_BS_m_FSR1 = ROOT.Double(0)
vtx_BS_m_FSR2 = ROOT.Double(0)
vtx_BS_phi_FSR1 = ROOT.Double(0)
vtx_BS_phi_FSR2 = ROOT.Double(0)

d0BS_vtx_BS1 = ROOT.Double(0)
d0BS_vtx_BS2 = ROOT.Double(0)

pT1_genFromReco = ROOT.Double(0)
pT2_genFromReco = ROOT.Double(0)

Tracker1 = ROOT.Double(0)
Tracker2 = ROOT.Double(0)

Tight1 = ROOT.Int(0)
Tight2 = ROOT.Int(0)

massZ2 = ROOT.Float(0)

# Vektoren
lep_id = ROOT.std.vector('int')()
lep_tightId = ROOT.std.vector('float')()
lep_RelIso = ROOT.std.vector('float')()
GENlep_hasFSR = ROOT.std.vector('float')()
GENlep_status = ROOT.std.vector('float')()
GENlep_MomId = ROOT.std.vector('float')()
GENlep_MomMomId = ROOT.std.vector('float')()
GENlep_id = ROOT.std.vector('float')()
GENlep_mass = ROOT.std.vector('float')()
GENlep_pt = ROOT.std.vector('float')()
pho_pt = ROOT.std.vector('float')()
fsrPhotons_pt = ROOT.std.vector('float')()
lep_pt = ROOT.std.vector('float')()
lep_pt_UnS = ROOT.std.vector('float')()
lep_pt_genFromReco = ROOT.std.vector('float')()
lep_pt_FromMuonBestTrack = ROOT.std.vector('float')()
lep_trackerLayersWithMeasurement = ROOT.std.vector('float')()
lepFSR_pt = ROOT.std.vector('float')()
vtxLep_pt = ROOT.std.vector('float')()
vtxLep_BS_pt_NoRoch = ROOT.std.vector('float')()
vtxLep_BS_pt = ROOT.std.vector('float')()
vtxLepFSR_pt = ROOT.std.vector('float')()
vtxLepFSR_BS_pt = ROOT.std.vector('float')()
GENlep_eta = ROOT.std.vector('float')()
lep_eta = ROOT.std.vector('float')()
pho_eta = ROOT.std.vector('float')()
fsrPhotons_eta = ROOT.std.vector('float')()
lepFSR_eta = ROOT.std.vector('float')()
vtxLep_eta = ROOT.std.vector('float')()
vtxLep_BS_eta = ROOT.std.vector('float')()
vtxLepFSR_eta = ROOT.std.vector('float')()
vtxLepFSR_BS_eta = ROOT.std.vector('float')()
GENlep_phi = ROOT.std.vector('float')()
lep_phi = ROOT.std.vector('float')()
pho_phi = ROOT.std.vector('float')()
fsrPhotons_phi = ROOT.std.vector('float')()
lepFSR_phi = ROOT.std.vector('float')()
vtxLep_phi = ROOT.std.vector('float')()
vtxLep_BS_phi = ROOT.std.vector('float')()
vtxLepFSR_phi = ROOT.std.vector('float')()
vtxLepFSR_BS_phi = ROOT.std.vector('float')()
lep_mass = ROOT.std.vector('float')()
lepFSR_mass = ROOT.std.vector('float')()
vtxLep_mass = ROOT.std.vector('float')()
vtxLep_BS_mass = ROOT.std.vector('float')()
vtxLepFSR_mass = ROOT.std.vector('float')()
vtxLepFSR_BS_mass = ROOT.std.vector('float')()
lep_numberOfValidPixelHits = ROOT.std.vector

lep_dataMC = ROOT.std.vector('float')()
lep_dataMCErr = ROOT.std.vector('float')()
lep_normalizedChi2 = ROOT.std.vector('float')()
lep_numberOfValidMuonHits = ROOT.std.vector('float')()
lep_numberOfMatchedStations = ROOT.std.vector('float')()
lep_Sip = ROOT.std.vector('float')()
lep_d0PV = ROOT.std.vector('float')()
lep_dzPV = ROOT.std.vector('float')()

commonPV_x = ROOT.std.vector('float')()
commonPV_y = ROOT.std.vector('float')()
commonPV_z = ROOT.std.vector('float')()
commonBS_x = ROOT.std.vector('double')()
commonBS_y = ROOT.std.vector('double')()
commonBS_z = ROOT.std.vector('double')()

PV_x = ROOT.Float(0)
PV_y = ROOT.Float(0)
PV_z = ROOT.Float(0)

BS_x = ROOT.Float(0)
BS_y = ROOT.Float(0)
BS_z = ROOT.Float(0)

BS_xErr = ROOT.Float(0)
BS_yErr = ROOT.Float(0)
BS_zErr = ROOT.Float(0)

BeamWidth_x = ROOT.Float(0)
BeamWidth_y = ROOT.Float(0)

pTL1 = ROOT.Float(0)
pTL2 = ROOT.Float(0)
pTL3 = ROOT.Float(0)
pTL4 = ROOT.Float(0)

etaL1 = ROOT.Float(0)
etaL2 = ROOT.Float(0)
etaL3 = ROOT.Float(0)
etaL4 = ROOT.Float(0)

phiL1 = ROOT.Float(0)
phiL2 = ROOT.Float(0)
phiL3 = ROOT.Float(0)
phiL4 = ROOT.Float(0)

mL1 = ROOT.Float(0)
mL2 = ROOT.Float(0)
mL3 = ROOT.Float(0)
mL4 = ROOT.Float(0)

met = ROOT.Float(0)

GENMH = ROOT.Float(0)

finalState = ROOT.Int(0)
passedFullSelection = ROOT.Bool(0)
is2P2F = ROOT.Int(0)
is3P1F = ROOT.Int(0)
isMCzz = ROOT.Int(0)
eventWeightFR = ROOT.Float(0)
eventWeightFR_up = ROOT.Float(0)
eventWeightFR_down = ROOT.Float(0)
fr3 = ROOT.Float(0)
fr2 = ROOT.Float(0)

mass4l = ROOT.Float(0)
mass4lREFIT = ROOT.Float(0)
mass4mu = ROOT.Float(0)
mass2e2mu = ROOT.Float(0)
mass4l_vtx = ROOT.Float(0)
mass4l_vtxFSR = ROOT.Float(0)
mass4lREFIT_vtx = ROOT.Float(0)
mass4l_vtx_BS = ROOT.Float(0)
mass4l_vtxFSR_BS = ROOT.Float(0)
mass4lREFIT_vtx_BS = ROOT.Float(0)

mass4lErr = ROOT.Float(0)
mass4lErrREFIT = ROOT.Float(0)
mass4lErr_vtx = ROOT.Float(0)
mass4lErr_vtx_BS = ROOT.Float(0)
mass4lErrREFIT_vtx = ROOT.Float(0)
mass4lErrREFIT_vtx_BS = ROOT.Float(0)

D_bkg_kin = ROOT.Float(0)
D_bkg_kin_vtx_BS = ROOT.Float(0)

massZ1 = ROOT.Float(0)
massZ1REFIT = ROOT.Float(0)

GENmassZ1 = ROOT.Float(0)

# Vektoren
lep_pterr = ROOT.std.vector('float')()
vtxLep_ptError = ROOT.std.vector('float')()
vtxLep_BS_ptError = ROOT.std.vector('float')()

nVtx = ROOT.Int(0)
nInt = ROOT.Int(0)
NVtx = ROOT.Int(0)
NInt = ROOT.Int(0)

directory = ROOT.TString("")
execute = ROOT.TString("")
filename = ROOT.TString("")
save_nome = ROOT.TString("")

c1 = ROOT.TCanvas()


prova_2 = []
prova_4mu = []
prova_4e = []
prova_2e2mu = []
prova_2mu2e = []

pT_bins = []
pT_bins_err = []
pT_bins_center = []
eta_bins = []
eta_bins_name = []
eta_bins_count = []
phi_bins = []

MC_mean_vs_gamma_DSCB = []
MC_meanErr_vs_gamma_DSCB = []
DATA_mean_vs_gamma_DSCB = []
DATA_meanErr_vs_gamma_DSCB = []

MC_mean_vs_gamma_CB = []
MC_meanErr_vs_gamma_CB = []
DATA_mean_vs_gamma_CB = []
DATA_meanErr_vs_gamma_CB = []

MC_mean_vs_pt_DSCB = []
MC_meanErr_vs_pt_DSCB = []
DATA_mean_vs_pt_DSCB = []
DATA_meanErr_vs_pt_DSCB = []

MC_mean_vs_pt_CB = []
MC_meanErr_vs_pt_CB = []
DATA_mean_vs_pt_CB = []
DATA_meanErr_vs_pt_CB = []

scale_DATA_MC_DSCB = []
scaleErr_DATA_MC_DSCB = []
scale_DATA_MC_CB = []
scaleErr_DATA_MC_CB = []

spreadScale_DSCB = []
spreadScaleErr_DSCB = []
spreadScale_CB = []
spreadErrScale_CB = []

spreadMean_DSCB = []
spreadMeanErr_DSCB = []
spreadMean_CB = []
spreadErrMean_CB = []

MC_mean = ROOT.Float(0)
DATA_mean = ROOT.Float(0)
MC_meanErr = ROOT.Float(0)
DATA_meanErr = ROOT.Float(0)

hist_lut = ROOT.TH2F()

DATA = ROOT.Bool(0)

fs = ROOT.TString("")
mode = ROOT.TString("")
histo_name = ROOT.TString("")
fit_param_latex = ROOT.TString("")
canvas_name = ROOT.TString("")

chi_dof = ROOT.Int(0)

lineRef = ROOT.TLine(60, 0, 120, 0)
lineRef_H = ROOT.TLine(105, 0, 140, 0)

_file0 = ROOT.TFile()
tree = ROOT.TTree()

lep_1 = ROOT.TLorentzVector()
lep_2 = ROOT.TLorentzVector()
lep_1err = ROOT.TLorentzVector()
lep_2err = ROOT.TLorentzVector()
ZPrime = ROOT.TLorentzVector()

scale_up = ROOT.TLine()
scale_down = ROOT.TLine()

mass_type = []

gamma_bins = []
gammaErr_bins = []

massPoint = []
massPointName = []
ProdMode = []
ProdMode_File = []
year = []
category = []
decayMode = []
luminosity = []

numberEvent_ggH = []
numberEvent_VBF = []
numberEvent_WplusH = []
numberEvent_WminusH = []
numberEvent_ZH = []
numberEvent_ttH = []
numberEvent = []
numberEvent_2018 = []

crossSection_ggH = []
crossSection_VBF = []
crossSection_WplusH = []
crossSection_WminusH = []
crossSection_ZH = []
crossSection_ttH = []
crossSection_bbH = []
crossSection_tHq = []
crossSection = []

color = ROOT.Color_t(0)
nome_file = ROOT.TString("")

x_min = ROOT.Float(60)
x_max = ROOT.Float(120)

massZ_min = ROOT.Float(60)
massZ_max = ROOT.Float(120)
massZErr_min = ROOT.Float(0.2)
massZErr_max = ROOT.Float(7.2)

BW_mean_PDG = ROOT.Float(91.1876)
BW_sigma_PDG = ROOT.Float(2.4952)

BW_mean_min = ROOT.Float(86)
BW_mean_max = ROOT.Float(96)
BW_mean_DATA_min = ROOT.Float(86)
BW_mean_DATA_max = ROOT.Float(96)
BW_mean_MC_min = ROOT.Float(86)
BW_mean_MC_max = ROOT.Float(96)

CB_mean_min = ROOT.Float(-5)
CB_mean_max = ROOT.Float(5)
CB_sigma_min = ROOT.Float(0)
CB_sigma_max = ROOT.Float(30)
CB_alpha_min = ROOT.Float(0)
CB_alpha_max = ROOT.Float(30)
CB_exp_min = ROOT.Float(0)
CB_exp_max = ROOT.Float(30)

DSCB_mean_min = ROOT.Float(-5)
DSCB_mean_max = ROOT.Float(5)
DSCB_sigma_min = ROOT.Float(0)
DSCB_sigma_max = ROOT.Float(5)
DSCB_alphaL_min = ROOT.Float(0)
DSCB_alphaL_max = ROOT.Float(20)
DSCB_expL_min = ROOT.Float(0)
DSCB_expL_max = ROOT.Float(20)
DSCB_alphaR_min = ROOT.Float(0)
DSCB_alphaR_max = ROOT.Float(10)
DSCB_expR_min = ROOT.Float(0)
DSCB_expR_max = ROOT.Float(10)

tau_min = ROOT.Float(-1)
tau_max = ROOT.Float(1)

fsig_min = ROOT.Float(0.01)
fsig_max = ROOT.Float(1)



def DrawResolution(h1, h2, h3, histos_name, nome_canvas, save, x_name):
    h1.SetStats(0)
    h2.SetStats(0)
    h3.SetStats(0)

    max_val = 1
    max_val = max(max_val, h1.GetMaximum())
    max_val = max(max_val, h2.GetMaximum())
    max_val = max(max_val, h3.GetMaximum())
    
    max_val *= 1.1

    canvas = ROOT.TCanvas(nome_canvas, nome_canvas, 750, 750)
    canvas.cd()
    pad11 = ROOT.TPad("pad1", "pad1", 0, 0.15, 1, 1)
    pad11.SetGrid()
    pad11.SetBottomMargin(0.1)
    pad11.Draw()
    pad11.cd()

    h1.SetTitle(nome_canvas)
    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(0.75)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.GetYaxis().SetRangeUser(0, max_val)
    h1.GetXaxis().SetTitle(x_name)

    h2.SetTitle(nome_canvas)
    h2.SetLineColor(ROOT.kRed)
    h2.SetMarkerStyle(20)
    h2.SetMarkerSize(0.75)
    h2.SetMarkerColor(ROOT.kRed)
    h2.GetYaxis().SetRangeUser(0, max_val)
    h2.GetXaxis().SetTitle(x_name)

    h3.SetTitle(nome_canvas)
    h3.SetLineColor(ROOT.kGreen)
    h3.SetMarkerStyle(20)
    h3.SetMarkerSize(0.75)
    h3.SetMarkerColor(ROOT.kGreen)
    h3.GetYaxis().SetRangeUser(0, max_val)
    h3.GetXaxis().SetTitle(x_name)

    ratio_2 = ROOT.TH1F("ratio_2", "Ratio 2", h2.GetNbinsX(), h2.GetXaxis().GetXmin(), h2.GetXaxis().GetXmax())
    ratio_2.Divide(h2, h1)
    ratio_2.GetYaxis().SetRangeUser(0.75, 1.25)
    ratio_2.SetTitle("")
    ratio_2.SetStats(0)
    ratio_2.GetYaxis().SetTitleSize(0.05)
    ratio_2.GetYaxis().SetLabelSize(0.14)
    label_Yaxis = f"{histos_name[1]} / {histos_name[0]}"
    ratio_2.GetYaxis().SetTitle(label_Yaxis)
    ratio_2.GetXaxis().SetTitle(x_name)

    ratio_3 = ROOT.TH1F("ratio_3", "Ratio 3", h3.GetNbinsX(), h3.GetXaxis().GetXmin(), h3.GetXaxis().GetXmax())
    ratio_3.Divide(h3, h1)
    ratio_3.GetYaxis().SetRangeUser(0.5, 1.05)
    ratio_3.SetTitle("")
    ratio_3.SetStats(0)
    ratio_3.GetYaxis().SetTitleSize(0.05)
    ratio_3.GetYaxis().SetLabelSize(0.14)
    ratio_3.GetYaxis().SetTitle(label_Yaxis)
    ratio_3.GetXaxis().SetTitle(x_name)

    legend = ROOT.TLegend(0.65, 0.15, 0.9, 0.35)
    legend.AddEntry(h1, histos_name[0])
    legend.AddEntry(h2, histos_name[1])
    legend.AddEntry(h3, histos_name[2])

    canvas.Update()
    h1.Draw("E")
    h2.Draw("same E")
    h3.Draw("same E")
    legend.Draw()

    line = ROOT.TLine(pad11.GetUxmin(), 1, pad11.GetUxmax(), 1)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)
    line.Draw()

    canvas.Update()

    canvas.cd()
    pad22 = ROOT.TPad("pad2", "pad2", 0, 0.001, 1, 0.15)
    pad22.Draw()
    pad22.cd()

    ratio_2.Draw("E")
    ratio_3.Draw("same E")

    line2 = ROOT.TLine(pad22.GetUxmin(), 1, pad22.GetUxmax(), 1)
    line2.SetLineColor(ROOT.kBlack)
    line2.SetLineWidth(1)
    line2.Draw()

    canvas.Update()

    canvas.Print(save)