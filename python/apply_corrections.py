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



ROOT.gInterpreter.Declare(
    """
    #include <boost/math/special_functions/erf.hpp>
    struct CrystalBall{
        double pi=3.14159;
        double sqrtPiOver2=sqrt(pi/2.0);
        double sqrt2=sqrt(2.0);
        double m;
        double s;
        double a;
        double n;
        double B;
        double C;
        double D;
        double N;
        double NA;
        double Ns;
        double NC;
        double F;
        double G;
        double k;
        double cdfMa;
        double cdfPa;
    CrystalBall():m(0),s(1),a(10),n(10){
        init();
    }
    CrystalBall(double mean, double sigma, double alpha, double n)
        :m(mean),s(sigma),a(alpha),n(n){
        init();
    }
    void init(){
        double fa = fabs(a);
        double ex = exp(-fa*fa/2);
        double A  = pow(n/fa, n) * ex;
        double C1 = n/fa/(n-1) * ex; 
        double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
        B = n/fa-fa;
        C = (D1+2*C1)/C1;   
        D = (D1+2*C1)/2;   
        N = 1.0/s/(D1+2*C1); 
        k = 1.0/(n-1);  
        NA = N*A;       
        Ns = N*s;       
        NC = Ns*C1;     
        F = 1-fa*fa/n; 
        G = s*n/fa;    
        cdfMa = cdf(m-a*s);
        cdfPa = cdf(m+a*s);
    }
    double pdf(double x) const{ 
    double d=(x-m)/s;
        if(d<-a) return NA*pow(B-d, -n);
        if(d>a) return NA*pow(B+d, -n);
        return N*exp(-d*d/2);
    }
    double pdf(double x, double ks, double dm) const{ 
        double d=(x-m-dm)/(s*ks);
        if(d<-a) return NA/ks*pow(B-d, -n);
        if(d>a) return NA/ks*pow(B+d, -n);
        return N/ks*exp(-d*d/2);

    }
    double cdf(double x) const{
        double d = (x-m)/s;
        if(d<-a) return NC / pow(F-s*d/G, n-1);
        if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
        return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
    }
    double invcdf(double u) const{
        if(u<cdfMa) return m + G*(F - pow(NC/u, k));
        if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
        return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
    }
    };
    """
)

ROOT.gInterpreter.Declare(
    """
    double cb_rndm(double mean, double sigma, double alpha, double n) {
        CrystalBall cb(mean, sigma, alpha, n);
        TRandom3 rnd(time(0));
        double rndm = gRandom->Rndm();
        //cout<< rndm << "  " << cb.invcdf(rndm) << std::endl;
        return cb.invcdf(rndm);
    }
    """
)



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
            double mean_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 1);\
            double sig_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 2);\
            double n_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 3);\
            double alpha_cb1 = h_results_cb->GetBinContent(etabin1, nlbin1, 4);\
            double rndm_cb1 = (double) cb_rndm(mean_cb1, sig_cb1, alpha_cb1, n_cb1);\
            double sig_poly1_a = h_results_poly->GetBinContent(etabin1, nlbin1, 1);\
            double sig_poly1_b = h_results_poly->GetBinContent(etabin1, nlbin1, 2);\
            double sig_poly1_c = h_results_poly->GetBinContent(etabin1, nlbin1, 3);\
            double sig_poly1 = sig_poly1_a + sig_poly1_b * genpt_1 + sig_poly1_c * genpt_1*genpt_1;\
            if (sig_poly1 < 0) sig_poly1 = 0;\
            pt1 = 1. / (1./genpt_1 * ( 1 + sig_poly1 * rndm_cb1));\
            return pt1;'
        )

        df = df.Define(
            "genpt_2_smeared",
            "double pt2;\
            Int_t etabin2 = h_results_cb->GetXaxis()->FindBin(abs(eta_2));\
            Int_t nlbin2 = h_results_cb->GetYaxis()->FindBin(nTrkLayers_2);\
            double mean_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 1);\
            double sig_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 2);\
            double n_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 3);\
            double alpha_cb2 = h_results_cb->GetBinContent(etabin2, nlbin2, 4);\
            double rndm_cb2 = (double) cb_rndm(mean_cb2, sig_cb2, alpha_cb2, n_cb2);\
            double sig_poly2_a = h_results_poly->GetBinContent(etabin2, nlbin2, 1);\
            double sig_poly2_b = h_results_poly->GetBinContent(etabin2, nlbin2, 2);\
            double sig_poly2_c = h_results_poly->GetBinContent(etabin2, nlbin2, 3);\
            double sig_poly2 = sig_poly2_a + sig_poly2_b * genpt_2 + sig_poly2_c * genpt_2*genpt_2;\
            if (sig_poly2 < 0) sig_poly2 = 0;\
            pt2 = 1. / (1./genpt_2 * ( 1 + sig_poly2 * rndm_cb2));\
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
                f"1. / (1./genpt_1 * (1 + h_k_{dtsg}->GetBinContent(h_k_{dtsg}->FindBin(abs(eta_1))) * (genpt_1_smeared/genpt_1 - 1)))"
            )
            df = df.Define(
                f"genpt_2_smeared_{dtsg}", 
                f"1. / (1./genpt_2 * (1 + h_k_{dtsg}->GetBinContent(h_k_{dtsg}->FindBin(abs(eta_2))) * (genpt_2_smeared/genpt_2 - 1)))"
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
            pt = 1. / (1./pt_1_roccor_it * (1 + sqrt(k_data*k_data - k_sig*k_sig) * (genpt_1_smeared/genpt_1 -1))); \
            return pt;"
        )
        df = df.Define(
            "pt_2_smeared",
            f"double pt; \
            double k_sig = h_k_SIG->GetBinContent(h_k_SIG->FindBin(abs(eta_2))); \
            double k_data = h_k_DATA->GetBinContent(h_k_DATA->FindBin(abs(eta_2))); \
            pt = 1. / (1./pt_2_roccor_it * (1 + sqrt(k_data*k_data - k_sig*k_sig) * (genpt_2_smeared/genpt_2 -1))); \
            return pt;"
        )
        df = df.Define(
            f"mass_Z_smeared",
            f"sqrt(2 * pt_1_smeared * pt_2_smeared * (cosh(eta_1 - eta_2) - cos(phi_1 - phi_2)));"
        )

    return df