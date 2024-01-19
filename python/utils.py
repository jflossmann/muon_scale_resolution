import numpy as np

def datahist_minus_bkghist(variable, df_data, df_mc, variable_bkg, dfs_bkg, bins, Range, normalized=True):
    #make weighted hists
    df_data["weight"]=df_data["zPtWeight"]*df_data["genWeight"]*df_data["sumwWeight"]*df_data["xsec"]*df_data["sf_id"]*df_data["sf_iso"]
    hist_data=list(np.histogram(df_data[variable], weights=df_data["weight"], range=Range, bins=bins))

    df_mc["weight"]=df_mc["zPtWeight"]*df_mc["genWeight"]*df_mc["sumwWeight"]*df_mc["xsec"]*df_mc["sf_id"]*df_mc["sf_iso"]
    hist_mc=np.histogram(df_mc[variable], weights=df_mc["weight"], range=Range, bins=bins)[0]

    hists_bkg=np.zeros(len(hist_data[0]))

    for df in dfs_bkg:
        df["weight"]=df["zPtWeight"]*df["genWeight"]*df["sumwWeight"]*df["xsec"]*df["sf_id"]*df["sf_iso"]
        hists_bkg+=np.histogram(df[variable_bkg], weights=df["weight"], range=Range, bins=bins)[0]

    hist_mc+=hists_bkg


    binwidth=hist_data[1][1]-hist_data[1][0]
    sf=np.sum(hist_data[0]*binwidth)/np.sum(hist_mc*binwidth)

    hist_data[0]=hist_data[0]-hists_bkg*sf

    if normalized:
        hist_data[0]=hist_data[0]/np.sum(hist_data[0]*binwidth)

    bincenters=[]
    for i in range(len(hist_data[0])):
        bincenters.append(hist_data[1][i]+(hist_data[1][i+1]-hist_data[1][i])/2)

    return hist_data, bincenters
