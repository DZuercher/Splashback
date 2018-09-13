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
    
    dat=pd.read_csv("%s/%s/xi_2d.dat" % (input_dir, types[0]), header=None,sep=' ')

    dat = np.asarray(dat.values)
    rr = dat[:,0]
    yy = dat[:,1]
    yy_err = dat[:,2]

    idx = (dat[:,0] > 0.1) & (dat[:,0] < 10.0)
    rr=rr[idx]
    yy=yy[idx]
    yy_err=yy_err[idx]
    rrange=np.linspace(rr[0],rr[-1],rsteps)

    f, axarr = plt.subplots(ncols=2,figsize=(30,15))

    for ii, type_ in enumerate(types):
	plotdat=pd.read_csv("%s/%s/plotsave.dat" % (output_dir, type_) ,sep=' ',header=None)
	plotdat=np.asarray(plotdat.values)
	sigmas=plotdat[0:rsteps,:]
	rhos=plotdat[rsteps:rsteps*2,:]
	devs=plotdat[rsteps*2:rsteps*3,:]
	rhodevs=plotdat[rsteps*3:rsteps*4,:]

	splashdat=pd.read_csv("%s/%s/values.dat" % (output_dir, type_),sep=' ',header=None)
	splashdat=np.asarray(splashdat.values)
	rsp2ds=splashdat[0,:]
	rsp3ds=splashdat[1,:]
	chisquares=splashdat[2,:]
	maxrhodevs=splashdat[3,:]
	maxinnerrhodevs=splashdat[4,:]


	axarr[0].plot(rrange,devs[:,0],'-',label="MCMC derivative")
	#axarr[1,0].plot(rrange,devs[:,0]+devs[:,1],'r--')
	#axarr[1,0].plot(rrange,devs[:,0]-devs[:,2],'r--')
	#axarr[0].axvline(rsp2ds[0],label="RSP 2D",lw=0.7)
	#axarr[1,0].axvline(rsp2ds[0]+rsp2ds[1],c='k',ls='dashed')
	#axarr[1,0].axvline(rsp2ds[0]-rsp2ds[2],c='k',ls='dashed')
	#axarr[0].axvline(rsp3ds[0],label="RSP 3D",alpha=0.5,lw=0.7)
	#axarr[1,0].axvline(rsp3ds[0]+rsp3ds[1],c='b',ls='dashed',alpha=0.5)
	#axarr[1,0].axvline(rsp3ds[0]-rsp3ds[2],c='b',ls='dashed',alpha=0.5)
	axarr[0].set_title("2D First Derivative")
	axarr[0].set_xscale("log")
	axarr[0].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
	axarr[0].set_ylabel(r"$\mathrm{d} log\xi^{\rm 2d} /\mathrm{d} logr$")
	#axarr[1,0].legend()


	axarr[1].plot(rrange,rhodevs[:,0],'-',label="MCMC derivative")
	#axarr[1,1].plot(rrange,rhodevs[:,0]+rhodevs[:,1],'r--')
	#axarr[1,1].plot(rrange,rhodevs[:,0]-rhodevs[:,2],'r--')
	#axarr[1].axvline(rsp2ds[0],label="RSP 2D",alpha=0.5,lw=0.7)
	#axarr[1,1].axvline(rsp2ds[0]+rsp2ds[1],c='k',ls='dashed',alpha=0.5)
	#axarr[1,1].axvline(rsp2ds[0]-rsp2ds[2],c='k',ls='dashed',alpha=0.5)
	#axarr[1].axvline(rsp3ds[0],label="RSP 3D",lw=0.7)
	#axarr[1,1].axvline(rsp3ds[0]+rsp3ds[1],c='b',ls='dashed')
	#axarr[1,1].axvline(rsp3ds[0]-rsp3ds[2],c='b',ls='dashed')
	axarr[1].set_title("3D First Derivative")
	axarr[1].set_xscale("log")
	axarr[1].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
	axarr[1].set_ylabel(r"$\mathrm{d} log\rho /\mathrm{d} logr$")
	#axarr[1,1].legend()

    #plt.tight_layout()
    plt.savefig("/work/dominik.zuercher/Output/temp/model_plots_combined.pdf" )

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
    types=['Planck_PS_21.5_blue_cut_1','Planck_PS_21.5_blue_cut_2','Planck_PS_21.5_blue_cut_3','Planck_PS_21.5_blue_cut_4','Planck_PS_21.5_blue_cut_5','Planck_PS_21.5_blue_cut_6']
    make_plots(types)
