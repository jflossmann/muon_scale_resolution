import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import crystalball
import uproot


def datahist_minus_bkghist(variable, df_data, df_mc, variable_bkg, dfs_bkg, bins, Range, normalized=True):
    #make weighted hists

    hist_data=list(np.histogram(df_data[variable], range=Range, bins=bins))

    hist_mc=np.histogram(df_mc[variable], weights=df_mc["weight"], range=Range, bins=bins)[0]
    hists_bkg=np.zeros(len(hist_data[0]))
    for df in dfs_bkg:
        
        hists_bkg+=np.histogram(df[variable_bkg], weights=df["weight"], range=Range, bins=bins)[0]

    hist_mc+=hists_bkg

    binwidth=hist_data[1][1]-hist_data[1][0]
    bincenters=[]
    for i in range(len(hist_data[0])):
        bincenters.append(hist_data[1][i]+(hist_data[1][i+1]-hist_data[1][i])/2)

    sf=np.sum(hist_data[0])/np.sum(hist_mc)

    #errorcheck before subtraction
    plt.plot(bincenters, hist_data[0],c="r")

    hist_data[0]=hist_data[0]-hists_bkg*sf

    if normalized:
        hist_data[0]=hist_data[0]/np.sum(hist_data[0]*binwidth)

    #errorcheck
    plt.plot(bincenters, hist_data[0],c="b")
    plt.plot(bincenters, hists_bkg*sf,c="k")
    plt.savefig("test.png")
    plt.clf()
    return hist_data[0], bincenters



def get_background_scale_factor(df_data, df_mc, dfs_bkg, bins, Range, normalized=True):
    #make weighted hists

    hist_data=list(np.histogram(df_data["mass_Z"], range=Range, bins=bins))


    hist_mc=np.histogram(df_mc["mass_Z"], weights=df_mc["weight"], range=Range, bins=bins)[0]
    hists_bkg=np.zeros(len(hist_data[0]))

    for df in dfs_bkg:

        hists_bkg+=np.histogram(df["mass_Z"], weights=df["weight"], range=Range, bins=bins)[0]

    hist_mc+=hists_bkg

    bincenters=[]
    for i in range(len(hist_data[0])):
        bincenters.append(hist_data[1][i]+(hist_data[1][i+1]-hist_data[1][i])/2)

    sf=np.sum(hist_data[0])/np.sum(hist_mc)

    return sf


def get_bkg_average_single(dfs_bkg, var, bins, Range,sf):

    #make weighted hists
    
    hists_bkg=np.zeros(bins)

    for df in dfs_bkg:
        hist_bkg = np.histogram(df[var], weights=df["weight"], range=Range, bins=bins)
        hists_bkg+=hist_bkg[0]

    bincenters=[]
    for i in range(len(hist_bkg[0])):
        bincenters.append(hist_bkg[1][i]+(hist_bkg[1][i+1]-hist_bkg[1][i])/2)

    
    mean = np.sum(np.array(bincenters)*hists_bkg)/np.sum(hists_bkg)
    
    histsum=np.sum(hists_bkg)
    n=histsum*sf
    
    return mean, n


def get_bkg_average_product(dfs_bkg, var1, var2, bins, Range):

    #make weighted hists

    hists_bkg=np.zeros(bins)

    for df in dfs_bkg:
        hist_bkg = np.histogram(df[var1]*df[var2], weights=df["weight"], range=Range, bins=bins)
        hists_bkg+=hist_bkg[0]

    bincenters=[]
    for i in range(len(hists_bkg)):
        bincenters.append(hist_bkg[1][i]+(hist_bkg[1][i+1]-hist_bkg[1][i])/2)

    mean = np.sum(np.array(bincenters)*hists_bkg)/np.sum(hists_bkg)
    
    return mean

def plot_histogram(ax, hist, color="k", alpha=1, label=None):
    #hist is output of datahist_minus_bkghist()
    binwidth=hist[1][1]-hist[1][0]
    ax.bar(hist[1], hist[0], width=binwidth, color=color, alpha=alpha, label=label)

def Poly2(x, a, b, c):
    return a*x*x + b*x +c

def CrystalBall(x, mu, sig, alpha, n):  
    return  crystalball.pdf(x=-x, beta=alpha, m=n, loc=-mu, scale=sig)

def Inv_mass(pt1, pt2, eta1, eta2, phi1, phi2):
    return np.sqrt( 2*pt1*pt2*(np.cosh(eta1-eta2)-np.cos(phi1-phi2)) )

def Minfinder(fun, args, bounds, n=10, N=10):
    i=0
    while i<N:
        X2=[]
        K=[]
        width=(bounds[1]-bounds[0])/n
        for k in np.linspace(bounds[0],bounds[1],n+1):
            X2.append(fun(k,*args))
            K.append(k)
        m=np.argmin(X2)

        if m==0:
            bounds=[bounds[0],K[m]+width]
        elif m==n:
            bounds=[K[m]-width,bounds[1]]
        else:
            bounds=[K[m]-width,K[m]+width]
        
        i+=1
        #print(m, K[m], X2[m])
    return K[m]

def load_from_root(fname):
    print(f"read {fname}")
    file = uproot.open(fname)
    tree = file["Events"]
    variables = tree.keys()
    df = tree.arrays(variables, library="pd")

    return df
