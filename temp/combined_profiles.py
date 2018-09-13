import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import frogress
import matplotlib.cm as cm

def make_plots(types):
    

    for ii, type_ in enumerate(types):
	dat=pd.read_csv("%s/%s/xi_2d.dat" % (input_dir, type_), header=None,sep=' ')

	dat = np.asarray(dat.values)
	rr = dat[:,0]
	yy = dat[:,1]
	yy_err = dat[:,2]

	idx = (dat[:,0] > 0.1) & (dat[:,0] < 10.0)
	rr=rr[idx]
	yy=yy[idx]
	yy_err=yy_err[idx]
	rrange=np.linspace(rr[0],rr[-1],rsteps)

	plt.errorbar(rr,yy,yy_err,fmt=".",label="Cut %s" % str(ii+1))

    plt.xlabel(r"$R$ ($h^{-1}$Mpc)")
    plt.ylabel(r"$\xi^{\rm 2d}$ ($h^{-1}$Mpc)")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.savefig("/work/dominik.zuercher/Output/temp/profiles_combined.pdf" )

    print("Plot done")




if __name__ == "__main__":

    input_dir = "/work/dominik.zuercher/Output/splashpipe"
    output_dir = "/work/dominik.zuercher/Output/Mest"

    modded=True
    mc=False

    ndim=8
    nwalkers=28 #98
    nwalkers_modded=28
    steps=100000
    rsteps=25
    types=['Planck_PS_21.5_blue_cut_1','Planck_PS_21.5_blue_cut_2','Planck_PS_21.5_blue_cut_3','Planck_PS_21.5_blue_cut_4','Planck_PS_21.5_blue_cut_5','Planck_PS_21.5_blue_cut_6','Planck_PS_21.5_blue_cut_7','Planck_PS_21.5_blue_cut_8','Planck_PS_21.5_blue_cut_9']
    make_plots(types)
