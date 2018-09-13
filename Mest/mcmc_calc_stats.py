#! coding=utf8
import numpy as np
import pandas as pd
import emcee
from emcee.utils import MPIPool, sample_ellipsoid
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import functools
import subprocess
import argparse
import subprocess as sub
from subprocess import call


def MCMC_calc_stats(ndim,rank,type_,add,mc):
    call("mkdir -p %s/%s/data_parts" % (output_dir, type_), shell = 1)
    if mc==True:
        add=add+"_mc"
        ndim+=2
    #Reading of the data
    print("Doing "+str(type_)+" with prior "+str(add))
    data = pd.read_csv("%s/%s/xi_2d.dat" % (input_dir, type_), header=None, sep = ' ')
    cov_data = pd.read_csv("%s/%s/xi_2d_cov.dat" % (input_dir, type_), header=None, sep = ' ')

    data = np.asarray(data.values)
    rr = data[:,0]
    yy = data[:,1]
    yy_err = data[:,2]
    idx = (data[:,0] > 0.1) & (data[:,0] < radmax)
    cov = np.asarray(cov_data.values)
    cov=np.transpose(cov[idx])[idx]
    rr=rr[idx]
    yy=yy[idx]
    yy_err=yy_err[idx]
    rrange=np.linspace(rr[0],rr[-1],rsteps)

    inv_cov=np.linalg.inv(cov)
    if diag==True:
        inv_cov=np.diag(np.diag(inv_cov))

    chunky=int(steps/size)
    rest=steps-chunky*size
    mini=chunky*rank
    maxi=chunky*(rank+1)
    if rank>=(size-1)-rest:
        maxi+=1+rank-(size-1-rest)
        mini+=rank-(size-1-rest)
    if rank==size-1:
        maxi=steps
    mini=int(mini)
    maxi=int(maxi)

    sigma_array = np.zeros((1,rrange.size))
    rho_array = np.zeros((1,rrange.size))
    dev_array=np.zeros((1,rrange.size))
    rhodev_array=np.zeros((1,rrange.size))
    rsp2d_array=np.zeros((1,1))
    rsp3d_array=np.zeros((1,1))
    chisquare_array=np.zeros((1,1))
    maxrhodev_array=np.zeros((1,1))
    maxinnerrhodev_array=np.zeros((1,1))

    data = np.loadtxt("%s/%s/chainstate_full.txt" % (output_dir, type_))
    
    data=np.asarray(data)
    
    length=data.size/nwalkers/ndim
    samples = data.reshape((length,nwalkers,ndim))
    chainstate = samples[mini:maxi,:,:]
    print("Starting Calculations... "+str(chainstate.shape))

    for i in range(chainstate.shape[0]):
        print("Doing step %s out of %s" % (i, chainstate.shape[0]))
        for j in range(chainstate.shape[1]):
            pact=chainstate[i,j,0:8]
            chisquare=tools.get_chisquare_reload(chainstate[i,j,:],rr,yy,yy_err,inv_cov,prior="None",mc=mc,alt=False)

            #For first plot
            if mc == True:
		sigma = tools.total_surf_den(rrange,chainstate[i,j,:],False)
            else:
		sigma=tools.surf_den(rrange,chainstate[i,j,:])
            sigma_array = np.vstack((sigma_array,sigma))
            rho = tools.einasto_profile(rrange,pact)
            rho_array = np.vstack((rho_array,rho))
            #For second plot
            dev = tools.logdev_surf_den(rrange,pact,recursive=True,surf=sigma)
            dev_array = np.vstack((dev_array,dev))
            rhodev=tools.logdev_einasto_profile(rrange,pact)
            rhodev_array=np.vstack((rhodev_array,rhodev))
            #Additional values (splashbacks and chisquares)
            rsp2d = tools.mindev2_procedure(pact)
            rsp2d_array=np.vstack((rsp2d_array,rsp2d))
            rsp3d = tools.mindev3_procedure(pact)
            rsp3d_array=np.vstack((rsp3d_array,rsp3d))
            chisquare_array=np.vstack((chisquare_array,chisquare))
            maxrhodev=tools.logdev_einasto_profile(rsp3d,pact)
            maxinnerrhodev=tools.logdev_einasto_profile_inner(rsp3d,pact)
            maxrhodev_array=np.vstack((maxrhodev_array,maxrhodev))
            maxinnerrhodev_array=np.vstack((maxinnerrhodev_array,maxinnerrhodev))
    out=np.hstack((sigma_array, rho_array, dev_array, rhodev_array, rsp2d_array, rsp3d_array, chisquare_array, maxrhodev_array, maxinnerrhodev_array))
    out=out[1:,:]
    np.savetxt("%s/%s/data_parts/MCMC_data_%s.txt" % (output_dir, type_, rank),out)
    print("Done!") 


if __name__ == "__main__":

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    steps = 11000
    mc = True
    rsteps = 25
    ndim = 8 #Number of parameters
    radmax = 10.0
    nwalkers = 28 #Total Number of walkers (98)
    diag = False

    input_dir = "/work/dominik.zuercher/Output/splashpipe"
    output_dir = "/work/dominik.zuercher/Output/Mest"

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--type_", help="Type",default="Planck_PS_21")
    parser.add_argument("--add", help="Prior",default="")

    args=parser.parse_args()

    MCMC_calc_stats(ndim,rank,args.type_,args.add,mc)
