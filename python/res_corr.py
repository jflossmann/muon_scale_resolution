def apply_res_corr(ntuples_gen, hdir, pdir, do_plot, do_binwise_plot=False):

    pdir = pdir+'resolution/'
    #open fit_result file
    #open json, read as dict
    with open(f"{hdir}/fit_results_res.json", 'r') as openfile: fit_results = json.load(openfile)
    
    #get bins from
    abseta_bins=fit_results["bins"]["abseta"]
    nl_bins=fit_results["bins"]["nl"]

    #open data
    file=uproot.open(ntuples_gen)
    tree=file["Events"]
    variables=tree.keys()
    df=tree.arrays(variables, library="pd")

    df_new1=pd.DataFrame()

    #loop over bins, fill df_new1
    
    print("getting corrections for pt 1")
    for i in tqdm(range(len(abseta_bins)-1)):
        eta_1_filter=(abseta_bins[i]<df["eta_1"]) & (df["eta_1"]<=abseta_bins[i+1])
        df_e1=df[eta_1_filter]

        for j in range(len(nl_bins)-1):
            nl_1_filter=(nl_bins[j]<df_e1["nTrkLayers_1"]) & (df_e1["nTrkLayers_1"]<=nl_bins[j+1])
            df_en1=df_e1[nl_1_filter]

            if len(df_en1)>0:
                #get CB sigma for current bin:
                
                CB_sigma=fit_results["pull"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"]
                #get polynomial fit parameters
                poly_par=np.array([fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 0"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 2"]["value"]])
                pt1=df_en1["genpt_1"]
                #calc correction

                #calc std from parabula
                std1=pt1*pt1*poly_par[2]+pt1*poly_par[1]+poly_par[0]
                #smear pt
                
                df_en1["genpt_1_smeared"]=np.array(df_en1["pt_1"]) + np.random.normal(0, abs(std1), size=len(df_en1))#*CB_sigma

                df_new1=pd.concat([df_new1, df_en1])
                if do_plot and do_binwise_plot:
                    bins=50
                    rang=[0,200]
                    fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
                    plt.subplots_adjust(wspace=0, hspace=0)
                    ax0.set_xlim(left=rang[0],right=rang[1])
                    ax0.set_xticks([])
                    h_s=ax0.hist(df_en1["genpt_1_smeared"],bins=bins, range=rang, histtype="step",density=True, label="smeared")
                    h_r=ax0.hist(df_en1["pt_1"], bins=bins, histtype="step", range=rang, density=True, label="MC_reco")
                    ax0.legend()
                    ax0.set_xlabel("pt_µ1 (GeV/c)")
                    ax0.set_title("without CB-smearing")
                    ax1.set_ylim(bottom=0.85,top=1.15)
                    ax1.set_xlim(left=rang[0],right=rang[1])
                    ax1.plot(np.linspace(rang[0],rang[1],bins),h_s[0]/h_r[0],'.')
                    ax1.set_xlabel("M_µµ (GeV)")
                    ax1.set_ylabel("ratio")
                    ax1.grid(True)
                    plt.savefig(f"{pdir}binwise_pt_Distribution/abseta{i}_nL{j}_1.png")
                    plt.clf()
        

    df_new2=pd.DataFrame()
    df=df_new1
    print("getting corrections for pt 2")
    for i in tqdm(range(len(abseta_bins)-1)):
        eta_2_filter=(abseta_bins[i]<df["eta_2"]) & (df["eta_2"]<=abseta_bins[i+1])
        df_e2=df[eta_2_filter]

        for j in range(len(nl_bins)-1):
            nl_2_filter=(nl_bins[j]<df_e2["nTrkLayers_2"]) & (df_e2["nTrkLayers_2"]<=nl_bins[j+1])
            df_en2=df_e2[nl_2_filter]

            if len(df_en2)>0:
                #get CB sigma for current bin:
                CB_sigma=fit_results["pull"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"]
                #get polynomial fit parameters
                poly_par=np.array([fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 0"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 1"]["value"],
                                fit_results["poly"]["eta_"+str(i)]["nL_"+str(j)]["Parameter 2"]["value"]] )
                
                pt2=df_en2["genpt_2"]
                #calc correction
                #calc std from parabula
                std2=pt2*pt2*poly_par[2]+pt2*poly_par[1]+poly_par[0]

                #smear pt
                df_en2["genpt_2_smeared"]=np.array(df_en2["pt_2"]) + np.random.normal(0, abs(std2), size=len(df_en2))#*CB_sigma 
                df_new2=pd.concat([df_new2, df_en2])
                if do_plot and do_binwise_plot:
                    bins=50
                    rang=[0,200]
                    fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
                    ax0.set_xlim(left=rang[0],right=rang[1])
                    ax0.set_xticks([])
                    plt.subplots_adjust(wspace=0, hspace=0)
                    h_s=ax0.hist(df_en2["genpt_1_smeared"],bins=bins, range=rang, histtype="step",density=True, label="smeared")
                    h_r=ax0.hist(df_en2["pt_1"], bins=bins, histtype="step", range=rang, density=True, label="MC_reco")
                    ax0.legend()
                    ax0.set_xlabel("M_µµ (GeV)")
                    ax0.set_title("without CB-smearing")
                    
                    ax1.set_ylim(bottom=0.85,top=1.15)
                    ax1.set_xlim(left=rang[0],right=rang[1])
                    ax1.plot(np.linspace(rang[0],rang[1],bins),h_s[0]/h_r[0],'.')
                    ax1.set_xlabel("pt_µ1 (GeV/c)")
                    ax1.set_ylabel("ratio")
                    ax1.grid(True)
                    plt.savefig(f"{pdir}binwise_pt_Distribution/abseta{i}_nL{j}_2.png")
                    plt.clf()
    

    #calculate invarriant mass
    df_new2["mass_Z_smeared_Jost"]=np.sqrt( 2*df_new2.smearedgenpt_1*df_new2.smearedgenpt_2*(np.cosh(df_new2.eta_1-df_new2.eta_2)-np.cos(df_new2.phi_1-df_new2.phi_2)) )
    df_new2["mass_Z_smeared_Dori"]=np.sqrt( 2*df_new2.genpt_1_smeared*df_new2.genpt_2_smeared*(np.cosh(df_new2.eta_1-df_new2.eta_2)-np.cos(df_new2.phi_1-df_new2.phi_2)) )

    if do_plot:
        bins=100
        rang=[80,100]
        fig, (ax0,ax1)=plt.subplots(2,1, gridspec_kw={"height_ratios":[3,1] })
        plt.subplots_adjust(wspace=0, hspace=0)
        h_s=ax0.hist(df_new2["mass_Z_smeared_Dori"],bins=bins, range=rang, histtype="step",density=True, label="smeared")
        h_r=ax0.hist(df_new2["mass_Z"],bins=bins, histtype="step", range=rang, density=True, label="MC_reco")
        #ax0.hist(df_new2["mass_Z_smeared_Jost"],bins=bins, range=rang, histtype="step",density=True, label="MC_gen") 
        chi2=np.sum((h_r[0]-h_s[0])**2/h_s[0])
        ndf=bins
        ax0.annotate(f"χ²/NDF = {round(chi2,6)}/{ndf}",xy=(100,400), xycoords="figure pixels")
        ax0.set_xlim(left=rang[0],right=rang[1])
        ax0.set_xticks([])
        ax0.legend()
        ax0.set_xlabel("M_µµ (GeV)")
        ax0.set_title("without CB-smearing")
        ax1.set_ylim(bottom=0.85,top=1.15)
        ax1.set_xlim(left=rang[0],right=rang[1])
        ax1.plot(np.linspace(rang[0],rang[1],bins),h_s[0]/h_r[0])
        ax1.set_xlabel("M_µµ (GeV)")
        ax1.set_ylabel("ratio")
        ax1.grid(True)
        #plt.savefig(f"{pdir}Z_mass_comparison.pdf")
        plt.savefig(f"{pdir}Z_mass_comparison.png")
        plt.clf()

    #save data
    print(f"saving corrected gen to {ntuples_gen.replace('.root', '_corr.root')}")
    data={key: df_new2[key].values for key in df_new2.columns}
    rdf = ROOT.RDF.MakeNumpyDataFrame(data)
    rdf.Snapshot("Events", ntuples_gen.replace('.root', '_corr.root'))
    print("done")

