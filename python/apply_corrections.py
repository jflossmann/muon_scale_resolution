import ROOT


def get_rdf(path):
    chain = ROOT.TChain("Events")

    chain.AddFriend("Events", path.split('/')[-1])

    #rdf = rdf.Define("weight", "zPtWeight*genWeight/sumwWeight*xsec*sf_id*sf_iso*bs_weight")

    #print(rdf.Mean("weight").GetValue())
    return chain


def step1(rdf, hdir, typ):
    if typ in ['GEN', 'BKG']:
        typ = 'SIG'

    #open Tfile with scale corrections
    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step1_C.root", "READ");')
    #read histograms for curent mode
    ROOT.gROOT.ProcessLine(f'TH2D* M_DATA = (TH2D*)tf->Get("M_DATA");')
    ROOT.gROOT.ProcessLine(f'TH2D* M_SIG = (TH2D*)tf->Get("M_SIG");')
    ROOT.gROOT.ProcessLine(f'TH2D* M_BKG = (TH2D*)tf->Get("M_SIG");')
    ROOT.gROOT.ProcessLine(f'TH2D* M_GEN = (TH2D*)tf->Get("M_SIG");')
    ROOT.gROOT.ProcessLine(f'TH2D* A_DATA = (TH2D*)tf->Get("A_DATA");')
    ROOT.gROOT.ProcessLine(f'TH2D* A_SIG = (TH2D*)tf->Get("A_SIG");')  
    ROOT.gROOT.ProcessLine(f'TH2D* A_BKG = (TH2D*)tf->Get("A_SIG");')
    ROOT.gROOT.ProcessLine(f'TH2D* A_GEN = (TH2D*)tf->Get("A_SIG");')

    rdf = rdf.Define(
        "M_1",
        f"M_{typ}->GetBinContent(M_{typ}->FindBin(eta_1, phi_1))"
    )
    rdf = rdf.Define(
        "M_2",
        f"M_{typ}->GetBinContent(M_{typ}->FindBin(eta_2, phi_2))"
    )
    rdf = rdf.Define(
        "A_1",
        f"A_{typ}->GetBinContent(A_{typ}->FindBin(eta_1, phi_1))"
    )
    rdf = rdf.Define(
        "A_2",
        f"A_{typ}->GetBinContent(A_{typ}->FindBin(eta_2, phi_2))"
    )
    rdf = rdf.Define(f"pt_1_roccor", f"1./(1./pt_1 * M_1 + charge_1 * A_1)")
    rdf = rdf.Define(f"pt_2_roccor", f"1./(1./pt_2 * M_2 + charge_2 * A_2)")

    #calculate corrected 4-momenta and corrected Z-mass
    rdf = rdf.Define(f"p4_1", f"ROOT::Math::PtEtaPhiMVector(pt_1_roccor, eta_1, phi_1, mass_1)")
    rdf = rdf.Define(f"p4_2", f"ROOT::Math::PtEtaPhiMVector(pt_2_roccor, eta_2, phi_2, mass_2)")
    rdf = rdf.Define(f"p4_Z", f"p4_1 + p4_2")
    rdf = rdf.Define(f"mass_Z_roccor", f"p4_Z.M()")

    return rdf



def step2(df, hdir, typ):
    if typ == 'GEN':
        # application of mean and sigma to create pull distributions
        ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step2_fitresults.root", "read");')
        ROOT.gROOT.ProcessLine(f'TH3D* h_results_cb = (TH3D*)tf->Get("h_results_cb");')
        ROOT.gROOT.ProcessLine(f'TH3D* h_results_poly = (TH3D*)tf->Get("h_results_poly");')
        #ROOT.gROOT.ProcessLine('tf->Close()')

        df = df.Define("R_1", "genpt_1/pt_1_roccor")
        df = df.Define("R_2", "genpt_2/pt_2_roccor")

        df = df.Filter(
            "abs(eta_1) < 2.4 && abs(eta_2) < 2.4 &&\
            nTrkLayers_1 > 6.5 && nTrkLayers_1 < 17.5 && \
            nTrkLayers_2 > 6.5 && nTrkLayers_2 < 17.5"
        )

        df = df.Define(
            "genpt_1_smeared",
            'double pt1;\
            Int_t etabin1 = h_results_cb->GetXaxis()->FindBin(abs(eta_1));\
            Int_t nlbin1 = h_results_cb->GetYaxis()->FindBin(nTrkLayers_1);\
            double sig_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 2);\
            double sig_poly1_a = h_results_poly->GetBinContent(etabin1, nlbin1, 1);\
            double sig_poly1_b = h_results_poly->GetBinContent(etabin1, nlbin1, 2);\
            double sig_poly1_c = h_results_poly->GetBinContent(etabin1, nlbin1, 3);\
            double sig_poly1 = sig_poly1_a + sig_poly1_b * genpt_1 + sig_poly1_c * genpt_1*genpt_1;\
            double sig1 = sig_cb1 * sig_poly1;\
            if (sig1 < 0) sig1 = 0;\
            pt1 = genpt_1 * ( 1 + sig1 * (float)(gaus()));\
            return pt1;'
        )

        df = df.Define(
            "genpt_2_smeared",
            "double pt2;\
            Int_t etabin2 = h_results_cb->GetXaxis()->FindBin(abs(eta_2));\
            Int_t nlbin2 = h_results_cb->GetYaxis()->FindBin(nTrkLayers_2);\
            double sig_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 2);\
            double sig_poly2_a = h_results_poly->GetBinContent(etabin2, nlbin2, 1);\
            double sig_poly2_b = h_results_poly->GetBinContent(etabin2, nlbin2, 2);\
            double sig_poly2_c = h_results_poly->GetBinContent(etabin2, nlbin2, 3);\
            double sig_poly2 = sig_poly2_a + sig_poly2_b * genpt_2 + sig_poly2_c * genpt_2*genpt_2;\
            double sig2 = sig_cb2 * sig_poly2;\
            if (sig2 < 0) sig2 = 0;\
            pt2 = genpt_2 * ( 1 + sig2 * (float)(gaus()));\
            return pt2;"
        )

        df = df.Define(
            "genmass_Z_smeared",
            "sqrt(2 * genpt_1_smeared * genpt_2_smeared * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)));"
        )

    return df



def step3(df, hdir, typ):
    if typ in ['GEN', 'BKG']:
        typ = 'SIG'

    ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step3_it_{typ}.root", "READ");')
    ROOT.gROOT.ProcessLine(f'TH2D* h_kappa_{typ} = (TH2D*)tf->Get("M_{typ}");')
    ROOT.gROOT.ProcessLine(f'h_kappa_{typ}->SetDirectory(nullptr);')
    ROOT.gROOT.ProcessLine(f'TH2D* h_lambd_{typ} = (TH2D*)tf->Get("A_{typ}");')
    ROOT.gROOT.ProcessLine(f'h_lambd_{typ}->SetDirectory(nullptr);')
    ROOT.gROOT.ProcessLine(f'tf->Close();')

    df = df.Define(
        f"pt_1_roccor_it",
        f"double pt;\
        pt = 1./ (h_kappa_{typ}->GetBinContent( h_kappa_{typ}->GetXaxis()->FindBin(eta_1) , h_kappa_{typ}->GetYaxis()->FindBin(phi_1) ) / pt_1 - \
        h_lambd_{typ}->GetBinContent( h_lambd_{typ}->GetXaxis()->FindBin(eta_1) , h_lambd_{typ}->GetYaxis()->FindBin(phi_1) ));\
        return pt;"
    )
    df = df.Define(
        f"pt_2_roccor_it",
        f"1./ (h_kappa_{typ}->GetBinContent( h_kappa_{typ}->GetXaxis()->FindBin(eta_2) , h_kappa_{typ}->GetYaxis()->FindBin(phi_2) ) / pt_2 + \
        h_lambd_{typ}->GetBinContent( h_lambd_{typ}->GetXaxis()->FindBin(eta_2) , h_lambd_{typ}->GetYaxis()->FindBin(phi_2) ))"
    )
    df = df.Define(
        f"mass_Z_roccor_it",
        f"sqrt( 2 * pt_1_roccor_it * pt_2_roccor_it * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)) )"
    )
    
    return df


def step4(df, hdir, typ):          
    if typ == 'GEN':
        ROOT.gROOT.ProcessLine(f'TFile* tf = TFile::Open("{hdir}step4_k.root", "read");')

        for dtsg in ['DATA', 'SIG']:

            ROOT.gROOT.ProcessLine(f'TH1D* h_k_{dtsg} = (TH1D*)tf->Get("h_k_{dtsg}");')

            df = df.Define(
                f"genpt_1_smeared_{dtsg}", 
                f"genpt_1 * (1 + h_k_{dtsg}->GetBinContent(h_k_{dtsg}->FindBin(abs(eta_1))) * (genpt_1_smeared/genpt_1 - 1))"
            )
            df = df.Define(
                f"genpt_2_smeared_{dtsg}", 
                f"genpt_2 * (1 + h_k_{dtsg}->GetBinContent(h_k_{dtsg}->FindBin(abs(eta_2))) * (genpt_2_smeared/genpt_2 - 1))"
            )
            df = df.Define(
                f"genmass_Z_smeared_{dtsg}",
                f"sqrt(2 * genpt_1_smeared_{dtsg} * genpt_2_smeared_{dtsg} * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)));"
            )

        # correct mc directly (without gen information)
        df = df.Define(
            "pt_1_smeared",
            f"double pt; \
            double k_sig = h_k_SIG->GetBinContent(h_k_SIG->FindBin(abs(eta_1))); \
            double k_data = h_k_DATA->GetBinContent(h_k_DATA->FindBin(abs(eta_1))); \
            pt = pt_1 * (1 + sqrt(k_data*k_data - k_sig*k_sig) * (genpt_1_smeared/genpt_1 -1)); \
            return pt;"
        )
        df = df.Define(
            "pt_2_smeared",
            f"double pt; \
            double k_sig = h_k_SIG->GetBinContent(h_k_SIG->FindBin(abs(eta_2))); \
            double k_data = h_k_DATA->GetBinContent(h_k_DATA->FindBin(abs(eta_2))); \
            pt = pt_2 * (1 + sqrt(k_data*k_data - k_sig*k_sig) * (genpt_2_smeared/genpt_2 -1)); \
            return pt;"
        )
        df = df.Define(
            f"mass_Z_smeared",
            f"sqrt(2 * pt_1_smeared * pt_2_smeared * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)));"
        )

    return df