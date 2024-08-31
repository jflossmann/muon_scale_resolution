import ROOT
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.ticker import MultipleLocator
from python.apply_corrections import step1, step2, step3, step4, uncertainties
from time import time
import uproot
from multiprocessing import Pool, RLock
import mplhep as hep
from scipy.stats import norm

def job_wrapper(args):
    return calc_bs_ntuples(*args)

ROOT.gROOT.SetBatch()

def corr(typ):
    if typ=='DATA': #to compare some data to start and get used to it
        df = ROOT.RDataFrame("Events", f"/ceph/jdriesch/rochester_v0/2023_nlo/{typ}_zPt.root")
        df = df.Filter("abs(eta_1) < 2.4 && abs(eta_2) < 2.4")
        hdir_ref=f"/work/jflossmann/bachelor/muon_scale_resolution/hists/2023_nlo/"

        pt1step3,pt2step3,mstep4=[],[],[]
        pt_mean,pt_std,m_mean,m_std=[],[],[],[]
        ev=10
        for k in range(ev):
            dft=df.Range(k)
            dft_ref=step1(dft,hdir_ref,typ,'ref')
            dft_ref=step2(dft_ref,hdir_ref,typ,'ref')
            dft_ref=step3(dft_ref,hdir_ref,typ,'ref')
            dft_ref=step4(dft_ref,hdir_ref,typ,'ref')
            pt1_ref_mean=dft_ref.Mean(f"pt_1_step3_ref").GetValue()
            pt2_ref_mean=dft_ref.Mean(f"pt_2_step3_ref").GetValue()
            m_ref_mean=dft_ref.Mean(f"mass_Z_step4_ref").GetValue()
            for num in tqdm(range(0,100)):
                hdir=f"/work/jflossmann/bachelor/muon_scale_resolution/hists/2023_nlo/condor/stat/{num}/"
                dft=step1(dft,hdir,typ, num)
                dft=step2(dft,hdir,typ, num)
                dft=step3(dft,hdir,typ, num)
                dft=step4(dft,hdir,typ, num)
                pt1step3.append(dft.Mean(f"pt_1_step3_{num}").GetValue())
                pt2step3.append(dft.Mean(f"pt_2_step3_{num}").GetValue())
                mstep4.append(dft.Mean(f"mass_Z_step4_{num}").GetValue())
            pt_mean.append(np.mean(pt1step3)/pt1_ref_mean),pt_std.append(np.std(pt1step3)/pt1_ref_mean)
            pt_mean.append(np.mean(pt2step3)/pt2_ref_mean),pt_std.append(np.std(pt2step3)/pt2_ref_mean)
            m_mean.append(np.mean(mstep4)/m_ref_mean),m_std.append(np.std(mstep4)/m_ref_mean)
    
        x = np.arange(1,len(pt_mean)+1)
        x = np.repeat(x, 2)

        plt.errorbar(x[::2],pt_mean[::2],yerr=pt_std[::2],fmt='o',label=r'$\mu^-$')  # muon
        plt.errorbar(x[1::2],pt_mean[1::2],yerr=pt_std[1::2],fmt='o',color='red',label=r'$\mu^+$')  # antomuonmuon
        plt.legend()
        plt.set_xlabel('Nummer des Muons')
        plt.set_ylabel(r'$\frac{\mu(p_T)}{p_T}$')
        plt.xaxis.set_major_locator(MultipleLocator(1))
        plt.tight_layout()
        plt.savefig('corrDATA.pdf')
        plt.show()

    if typ=='SIG': #creating all the bootstrap samples
        bs_samples=1000

        arguments = [(num, typ) for num in range(800,bs_samples)]
        nthreads = 32

        # for a in tqdm(arguments[627:]):
        #     job_wrapper(a)
        
        pool = Pool(nthreads, initargs=(RLock(),), initializer=tqdm.set_lock)
        for _ in tqdm(
            pool.imap_unordered(job_wrapper, arguments),
            total=len(arguments),
            desc="Total progess",
            dynamic_ncols=True,
            leave=True,
            ):
            pass

        # calc_bs_ntuples(num, df, typ)


def calc_bs_ntuples(num, typ):

    df = ROOT.RDataFrame("Events", f"/ceph/jdriesch/rochester_v0/2023_nlo/DY_zPt.root")
    df = df.Filter("abs(eta_1) < 2.4 && abs(eta_2) < 2.4")

    dft=df.Range(0,100000)
    quants = []
    hdir=f"/work/jflossmann/bachelor/muon_scale_resolution/hists/2023_nlo/condor/stat/{num}/"
    dft=step1(dft, hdir, typ, f'_{num}')
    dft=step2(dft, hdir, typ, f'_{num}')
    dft=step3(dft, hdir, typ, f'_{num}')
    quants += [f"pt_1_step3_{num}_0", f"pt_2_step3_{num}_0"]
    for i in range(0,100):
        dft = dft.Define(f"pt_1_step3_{num}_{i}", f"pt_1_step3_{num}")
        dft = dft.Define(f"pt_2_step3_{num}_{i}", f"pt_2_step3_{num}")
        dft=step4(dft,hdir,typ, f'_{num}_{i}')
        quants += [f"pt_1_step4_{num}_{i}", f"pt_2_step4_{num}_{i}"]
        
    dft.Snapshot("Events", f"/ceph/jflossmann/bootstraps1000/2023_nlo/sig_processed_{num}.root", quants)
    print("Snapshot saved")
    return

def to_uproot(): #creating the refrence dataset to compare the uncertainties later
    typ='SIG'
    hdir_ref=f"/work/jflossmann/bachelor/muon_scale_resolution/hists/2023_nlo/"
    dft = ROOT.RDataFrame("Events", f"/ceph/jdriesch/rochester_v0/2023_nlo/DY_zPt.root")
    dft = dft.Filter("abs(eta_1) < 2.4 && abs(eta_2) < 2.4")
    dft=dft.Range(0,100000)
    dft_ref=step1(dft,hdir_ref,typ,'_ref')
    dft_ref=step2(dft_ref,hdir_ref,typ,'_ref')
    dft_ref=step3(dft_ref,hdir_ref,typ,'_ref')
    quants = [
        f"pt_1_step3_ref_0", f"pt_2_step3_ref_0",
        f"pt_1_step4_scaleup_ref_0", f"pt_2_step4_scaleup_ref_0",
        f"pt_1_step4_scaledn_ref_0", f"pt_2_step4_scaledn_ref_0",
                 "eta_1",                       "eta_2",
                 "phi_1",                       "phi_2",
                  "nTrkLayers_1", "nTrkLayers_2"
    ]
    for i in range(100):
        dft_ref=dft_ref.Define(f'pt_1_step3_ref_{i}', 'pt_1_step3_ref')
        dft_ref=dft_ref.Define(f'pt_2_step3_ref_{i}', 'pt_2_step3_ref')
        dft_ref=step4(dft_ref,hdir_ref,typ,f'_ref_{i}')
        dft_ref=uncertainties(dft_ref,hdir_ref,typ,f'_ref_{i}')
        quants += [
            f"pt_1_step4_ref_{i}", f"pt_2_step4_ref_{i}",
            f"pt_1_step4_resolup_ref_{i}", f"pt_2_step4_resolup_ref_{i}",
            f"pt_1_step4_resoldn_ref_{i}", f"pt_2_step4_resoldn_ref_{i}"
        ]
    dft_ref.Snapshot('Events','/ceph/jflossmann/bootstraps100/2023_nlo/sig_ref.root', quants)
    print("Snapshot saved")


def evaluate(): #all values or arrays coming from error propagation are labeled with unc at the start, whereas coming from the bootstrap method are labeled with boot

    input_columns = ["pt_1_step4_scaleup_ref_0", "pt_1_step4_scaledn_ref_0", "pt_1_step3_ref_0", "eta_1", "eta_2", "phi_1", "phi_2","nTrkLayers_1", "nTrkLayers_2", "pt_2_step4_scaleup_ref_0", "pt_2_step4_scaledn_ref_0", "pt_2_step3_ref_0"]
    input_columns += [f"pt_1_step4_ref_{i}" for i in range(100)]
    input_columns += [f"pt_2_step4_ref_{i}" for i in range(100)]
    input_columns += [f"pt_1_step4_resolup_ref_{i}" for i in range(100)]
    input_columns += [f"pt_2_step4_resolup_ref_{i}" for i in range(100)]
    input_columns += [f"pt_1_step4_resoldn_ref_{i}" for i in range(100)]
    input_columns += [f"pt_2_step4_resoldn_ref_{i}" for i in range(100)]

    with uproot.open('/ceph/jflossmann/bootstraps100/2023_nlo/sig_ref.root:Events') as f:
        pd_df_ref = f.arrays(input_columns, library='pd')

    prob1_resol,prob1_scale = [],[]
    prob2_resol,prob2_scale = [],[]
    eventcount = 10000
    bscount = 1000
    boot1_scale, boot1_res_std = [],[]
    boot2_scale, boot2_res_std = [],[]
    for i in tqdm(range(bscount)): 
        boot1_res_std.append([])
        boot2_res_std.append([])

        input_columns = [f"pt_1_step4_{i}_{j}"  for j in range(100)]
        input_columns += [f"pt_2_step4_{i}_{j}"  for j in range(100)]
        input_columns += [f"pt_1_step3_{i}_0"]
        input_columns += [f"pt_2_step3_{i}_0"]

        with uproot.open(f'/ceph/jflossmann/bootstraps1000/2023_nlo/sig_processed_{i}.root:Events') as f:
            pd_df = f.arrays(input_columns, library='pd')
    
        boot1_scale.append(pd_df[f'pt_1_step3_{i}_0'][:eventcount])
        boot2_scale.append(pd_df[f'pt_2_step3_{i}_0'][:eventcount])

        # pts2_bs.append(dft.Mean(f"pt_2_step3_{i}").GetValue()[0])
        row1_res,row2_res=[],[]
        for n in range(100):
            row1_res.append(pd_df[f'pt_1_step4_{i}_{n}'][:eventcount])
            row2_res.append(pd_df[f'pt_2_step4_{i}_{n}'][:eventcount])
        row1_res_tp = np.transpose(row1_res)
        row2_res_tp = np.transpose(row2_res)
        for n in range(eventcount):
            boot1_res_std[-1].append(np.std(row1_res_tp[n]))
            boot2_res_std[-1].append(np.std(row2_res_tp[n]))

    boot1_res_std_tp = np.transpose(boot1_res_std)
    boot2_res_std_tp = np.transpose(boot2_res_std)

    boot1_res_std_tp_f, boot2_res_std_tp_f = [],[]
    abseta1, abseta2 = [],[]
    testeta, testphi, testnl = [],[],[]
    eta_ref, nl_ref =[], []
    nbelow50 = []
    for n in range(eventcount):
        unc1_res_up,unc1_res_down,unc1_res_nom=[],[],[]
        unc2_res_up,unc2_res_down,unc2_res_nom=[],[],[]

        boot1_res_std_tp_f.append([])
        boot2_res_std_tp_f.append([])
        stds1_filter, stds2_filter = [],[]
        for std in boot1_res_std_tp[n]:
            if std > 1e-10:
                stds1_filter.append(std)
        for std in boot2_res_std_tp[n]:
            if std > 1e-10:
                stds2_filter.append(std)
        boot1_res_std_tp_f[n] = stds1_filter
        boot2_res_std_tp_f[n] = stds2_filter

        abseta1.append(np.abs(pd_df_ref["eta_1"][n]))
        abseta2.append(np.abs(pd_df_ref["eta_2"][n]))

        unc1_scale_nom_4 = pd_df_ref[f"pt_1_step4_ref_0"][n]
        unc2_scale_nom_4 = pd_df_ref[f"pt_2_step4_ref_0"][n]
        unc1_scale_up_4 = pd_df_ref[f"pt_1_step4_scaleup_ref_0"][n] 
        unc2_scale_up_4 = pd_df_ref[f"pt_2_step4_scaleup_ref_0"][n] 
        unc1_scale_dn_4 = pd_df_ref[f"pt_1_step4_scaledn_ref_0"][n]
        unc2_scale_dn_4 = pd_df_ref[f"pt_2_step4_scaledn_ref_0"][n]

        unc1_scale_nom = pd_df_ref[f"pt_1_step3_ref_0"][n] 
        unc1_scale_up = unc1_scale_nom + (unc1_scale_up_4 - unc1_scale_nom_4)
        unc1_scale_dn = unc1_scale_nom + (unc1_scale_dn_4 - unc1_scale_nom_4)

        unc2_scale_nom = pd_df_ref[f"pt_2_step3_ref_0"][n] 
        unc2_scale_up = unc2_scale_nom + (unc2_scale_up_4 - unc2_scale_nom_4)
        unc2_scale_dn = unc2_scale_nom + (unc2_scale_dn_4 - unc2_scale_nom_4)
        
        for i in range(100):
            unc1_res_nom.append(pd_df_ref[f"pt_1_step4_ref_{i}"][n])
            unc1_res_up.append(pd_df_ref[f"pt_1_step4_resolup_ref_{i}"][n])
            unc1_res_down.append(pd_df_ref[f"pt_1_step4_resoldn_ref_{i}"][n])

        for i in range(100):
            unc2_res_nom.append(pd_df_ref[f"pt_2_step4_ref_{i}"][n])
            unc2_res_up.append(pd_df_ref[f"pt_2_step4_resolup_ref_{i}"][n])
            unc2_res_down.append(pd_df_ref[f"pt_2_step4_resoldn_ref_{i}"][n])

        boot1_scale_n = [boot1_scale[i][n] for i in range(bscount)]
        boot2_scale_n = [boot2_scale[i][n] for i in range(bscount)]

        p1_resol=len([x for x in boot1_res_std_tp_f[n] if x>np.std(unc1_res_down) and x<np.std(unc1_res_up)])/len(boot1_res_std_tp_f[n])
        p1_scale=len([x for x in boot1_scale_n if x>unc1_scale_dn and x<unc1_scale_up])/len(boot1_scale_n)
        prob1_resol.append(p1_resol)
        prob1_scale.append(p1_scale)
        if p1_resol < 0.5:
            testeta.append(np.abs(pd_df_ref["eta_1"][n]))
            testphi.append(pd_df_ref["phi_1"][n])
            testnl.append(pd_df_ref["nTrkLayers_1"][n])
            nbelow50.append()
        eta_ref.append(np.abs(pd_df_ref["eta_1"][n]))
        nl_ref.append(pd_df_ref["nTrkLayers_1"][n])

        p2_resol=len([x for x in boot2_res_std_tp_f[n] if x>np.std(unc2_res_down) and x<np.std(unc2_res_up)])/len(boot2_res_std_tp_f[n])
        p2_scale=len([x for x in boot2_scale_n if x>unc2_scale_dn and x<unc2_scale_up])/len(boot2_scale_n)
        prob2_resol.append(p2_resol)
        prob2_scale.append(p2_scale)

        if n in nbelow50:
            hep.style.use("CMS")
            hep.cms.label(loc=0)
            count_res, bins_res, ignored =plt.hist(boot1_res_std_tp_f[n], bins=20, label = 'Data',color='blue')
            bin_centers_res = (bins_res[:-1] + bins_res[1:]) / 2
            errors_res=[]
            for i in count_res:
                if i == 0:
                   errors_res.append(0)
                else:
                    errors_res.append(1/np.sqrt(i)* i)
            plt.errorbar(bin_centers_res, count_res, yerr=errors_res, fmt='|', color='black')
            mu_res, sigma_res = np.mean(boot1_res_std_tp_f[n]), np.std(boot1_res_std_tp_f[n])
            x_res = np.linspace(min(bins_res), max(bins_res), 1000)
            g_res = norm.pdf(x_res, mu_res, sigma_res)
            g_res *= (count_res.max() / g_res.max())

            gaus_res = norm.pdf(bin_centers_res, mu_res, sigma_res)
            gaus_res *= (count_res.max() / gaus_res.max())
            chi2_res = np.sum(((count_res - gaus_res) ** 2) / gaus_res)  
            xndof_res= chi2_res/len(count_res)

            plt.plot(x_res, g_res, color = 'green', label = 'gaussian distribution')
            plt.axvline(np.std(unc1_res_up), color='red')
            plt.axvline(np.std(unc1_res_down), color='red', label = r'$\sigma^{up/down}$ ')
            plt.xlabel(r'$ \sigma_{res}$ in GeV')
            plt.ylabel(r'Counts')
            plt.legend(loc='upper left')
            plt.text(0.95, 0.95, f"$\chi^2 / ndof = {xndof_res:.2f}$", transform=plt.gca().transAxes, fontsize=12, va='top', ha='right')
            plt.text(0.95, 0.65, f"$\sigma = {sigma_res:.2f}$", transform=plt.gca().transAxes, fontsize=12, va='top', ha='right')
            plt.savefig(f'Histresbelow50/Histbelow50_res_{n}.pdf')
            plt.close()

            #hep.style.use("CMS")
            #hep.cms.label(loc=0)
            #count_scale, bins_scale, ignored =plt.hist(boot1_scale_n, label = 'Data', color='blue', bins=20)
            #bin_centers_scale = (bins_scale[:-1] + bins_scale[1:]) / 2
            #errors_scale=[]
            #for i in count_scale:
            #    if i == 0:
            #       errors_scale.append(0)
            #    else:
            #        errors_scale.append(1/np.sqrt(i)*i)
            #plt.errorbar(bin_centers_scale, count_scale, yerr=errors_scale, fmt='|', color='black')
            #mu_scale, sigma_scale = np.mean(boot1_scale_n), np.std(boot1_scale_n)
            #x_scale = np.linspace(min(bins_scale), max(bins_scale), 1000)
            #g_scale = norm.pdf(x_scale, mu_scale, sigma_scale)
            #g_scale *= (count_scale.max() / g_scale.max())

            #gaus_scale = norm.pdf(bin_centers_scale, mu_scale, sigma_scale) 
            #gaus_scale *= (count_scale.max() / gaus_scale.max())
            #chi2_scale = np.sum(((count_scale - gaus_scale) ** 2) / gaus_scale)  
            #xndof_scale= chi2_scale/len(count_scale)

            #plt.plot(x_scale, g_scale, color = 'green', label ='gaussian distribution')
            #plt.xlabel(r'$ p_T^{scale}$ in GeV')
            #plt.ylabel(r'Counts')
            #plt.axvline(unc1_scale_up,color='red', label = r'$p_T^{up/down}$ ')
            #plt.axvline(unc1_scale_dn,color='red')
            #plt.legend(loc='upper left')
            #plt.text(0.95, 0.95, f"$\chi^2 / ndof = {xndof_scale:.2f}$", transform=plt.gca().transAxes, fontsize=12, va='top', ha='right')
            #plt.savefig(f'HistogrammeUncVergleich1000/Hist_scale_{n}.pdf')
            #plt.close()
    
    hep.style.use("CMS")
    hep.cms.label(loc = 0)
    plt.hist(prob1_resol, color='red', histtype='step', label = 'Probability Resolution')
    plt.hist(prob1_scale, color='green', histtype='step', label = 'Probability Scale')
    plt.xlabel(r'$\frac{x_{in}}{x_{all}}$')
    plt.ylabel('Counts')
    plt.legend()
    plt.savefig('prob1000.pdf')
    plt.close()

    abseta_bins = np.linspace(0, 2.4, 13)
    prob_bins = np.linspace(0, 1, 13)
    nl_bins = np.linspace(8, 17, 10)

    counts, xedges, yedges, im = plt.hist2d(abseta1, prob1_resol, bins = (abseta_bins, prob_bins), cmap='viridis', density = True)
    norm_counts = counts / counts.sum(axis=1, keepdims=True)
    plt.imshow(norm_counts.T, origin='lower', aspect='auto', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='viridis')
    plt.colorbar(label=' normalized counts')
    plt.xlabel(r'|$\eta$|')
    plt.ylabel('Probability Resolution')
    plt.ylim(0,1)
    plt.savefig('prob1_binned_1000.pdf')
    plt.close()

    counts_ref, _, _, _ = plt.hist2d(eta_ref, nl_ref, cmap='viridis', bins = (abseta_bins, nl_bins))
    plt.close()
    counts1, xedges1, yedges1, _ = plt.hist2d(testeta, testnl, cmap= 'viridis', bins = (abseta_bins, nl_bins))
    plt.close()
    counts_ref_flat = counts_ref.flatten()
    counts_flat = counts1.flatten()
    counts_norm = []
    for i, count in enumerate(counts_ref_flat):
        if count == 0:
            counts_norm.append(0) 
        else:
            counts_norm.append(counts_flat[i] / count)
    counts_norm = np.array(counts_norm).reshape(counts1.shape)
    plt.imshow(counts_norm.T, origin='lower', extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]], aspect='auto', cmap='viridis')
    plt.colorbar(label=r'$\frac{events below 50%}{all events}$')
    plt.xlabel(r'|$\eta$|')
    plt.ylabel(r'$n_\mathrm{L}$')
    plt.savefig('eta-nl-norm1000_v1.pdf')
    plt.close()

# corr('SIG')
# to_uproot()
evaluate()