import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from python.utils import datahist_minus_bkghist


kappa_range=np.array([0.99, 1.01])
lambda_range=np.array([0.5,-0.5])*10**-3

def plot_lambda_fixed_eta(eta_bins, phi_bins, hdir, pdir):

    #make plotdir
    pdir = pdir+'closure/lambda/'
    os.makedirs(pdir, exist_ok=True)

    l_data = np.loadtxt(f"{hdir}lambda_data.csv", delimiter = ",")
    l_mc   = np.loadtxt(f"{hdir}lambda_mc.csv", delimiter = ",")

    for i in range(len(eta_bins)-1):
        eta_bin = [eta_bins[i],eta_bins[i+1]]

        for j in range(len(phi_bins)-1):
            if j==0:
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [l_mc[i][j],l_mc[i][j]],
                     c="b", label="mc")
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [l_data[i][j],l_data[i][j]],
                     c="r", label="data")
            else:
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [l_mc[i][j],l_mc[i][j]],
                     c="b")
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [l_data[i][j],l_data[i][j]],
                     c="r")
        plt.xticks(phi_bins)
        plt.xlim(left=min(phi_bins), right=max(phi_bins))
        plt.ylim(bottom=lambda_range[0], top=lambda_range[1])
        plt.xlabel("phi")
        plt.ylabel("lambda")
        plt.title(f"eta[{round(eta_bins[i],2)},{round(eta_bins[i+1],2)}]")
        plt.legend()
        plt.savefig(f"{pdir}lambda_eta_{i}.png")
        plt.clf()


def plot_lambda_fixed_phi(eta_bins, phi_bins, hdir, pdir):

    #make plotdir
    pdir = pdir+'closure/lambda/'
    os.makedirs(pdir, exist_ok=True)

    l_data = np.loadtxt(f"{hdir}lambda_data.csv", delimiter = ",")
    l_mc   = np.loadtxt(f"{hdir}lambda_mc.csv", delimiter = ",")

    for j in range(len(phi_bins)-1):
        phi_bin = [phi_bins[j], phi_bins[j+1]]

        for i in range(len(eta_bins)-1):
            if i==0:
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [l_mc[i][j],l_mc[i][j]],
                     c="b", label="mc")
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [l_data[i][j],l_data[i][j]],
                     c="r", label="data")
            else:
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [l_mc[i][j],l_mc[i][j]],
                     c="b")
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [l_data[i][j],l_data[i][j]],
                     c="r")
        plt.xticks(eta_bins)
        plt.xlim(left=min(eta_bins), right=max(eta_bins))
        plt.ylim(bottom=lambda_range[0], top=lambda_range[1])
        plt.xlabel("eta")
        plt.ylabel("lambda")
        plt.title(f"phi[{round(phi_bins[j],2)},{round(phi_bins[j+1],2)}]")
        plt.legend()
        plt.savefig(f"{pdir}lambda_phi_{j}.png")
        plt.clf()

def plot_kappa_fixed_eta(eta_bins, phi_bins, hdir, pdir):

    #make plotdir
    pdir = pdir+'closure/kappa/'
    os.makedirs(pdir, exist_ok=True)

    k_data = np.loadtxt(f"{hdir}kappa_data.csv", delimiter = ",")
    k_mc   = np.loadtxt(f"{hdir}kappa_mc.csv", delimiter = ",")

    for i in range(len(eta_bins)-1):
        for j in range(len(phi_bins)-1):
            if j==0:
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [k_mc[i][j],k_mc[i][j]],
                     c="b", label="mc")
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [k_data[i][j],k_data[i][j]],
                     c="r", label="data")
            else:
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [k_mc[i][j],k_mc[i][j]],
                     c="b")
                plt.plot([phi_bins[j], phi_bins[j+1]],
                     [k_data[i][j],k_data[i][j]],
                     c="r")
        plt.xticks(phi_bins)
        plt.xlim(left=min(phi_bins), right=max(phi_bins))
        plt.ylim(bottom=kappa_range[0], top=kappa_range[1])
        plt.xlabel("phi")
        plt.ylabel("kappa")
        plt.title(f"eta[{round(eta_bins[i],2)},{round(eta_bins[i+1],2)}]")
        plt.legend()
        plt.savefig(f"{pdir}kappa_eta_{i}.png")
        plt.clf()



def plot_kappa_fixed_phi(eta_bins, phi_bins, hdir, pdir):

    #make plotdir
    pdir = pdir+'closure/kappa/'
    os.makedirs(pdir, exist_ok=True)

    k_data = np.loadtxt(f"{hdir}kappa_data.csv", delimiter = ",")
    k_mc   = np.loadtxt(f"{hdir}kappa_mc.csv", delimiter = ",")

    for j in range(len(phi_bins)-1):
        phi_bin = [phi_bins[j], phi_bins[j+1]]

        for i in range(len(eta_bins)-1):
            if i==0:
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [k_mc[i][j],k_mc[i][j]],
                     c="b", label="mc")
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [k_data[i][j],k_data[i][j]],
                     c="r", label="data")
            else:
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [k_mc[i][j],k_mc[i][j]],
                     c="b")
                plt.plot([eta_bins[i], eta_bins[i+1]],
                     [k_data[i][j],k_data[i][j]],
                     c="r")
        plt.xticks(eta_bins)
        plt.xlim(left=min(eta_bins), right=max(eta_bins))
        plt.ylim(bottom=kappa_range[0], top=kappa_range[1])
        plt.xlabel("eta")
        plt.ylabel("kappa")
        plt.title(f"phi[{round(phi_bins[j],1)},{round(phi_bins[j+1],2)}]")
        plt.legend()
        plt.savefig(f"{pdir}kappa_phi_{j}.png")
        plt.clf()
