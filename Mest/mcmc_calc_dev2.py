#! coding=utf8
import numpy as np
import pandas as pd
import emcee
from emcee.utils import MPIPool, sample_ellipsoid
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/Splashback/toolbox')
import tools
import functools
import subprocess
import argparse
import subprocess as sub
from subprocess import call


def MCMC_calc_stats(ndim,rank,type_,add,mc):

    #Reading of the data
    print("Doing "+str(type_)+" with prior "+str(add))

    data = pd.read_csv("%s/%s/xi_2d.dat" % (input_dir, type_), header=None, sep = ' ')
    cov_data = pd.read_csv("%s/%s/xi_2d_cov.dat" % (input_dir, type_), header=None, sep = ' ')

    if mc == False:
	type_ += '_no_mc'

    call("mkdir -p %s/%s/dev2_parts" % (output_dir, type_), shell = 1)

    if mc==True:
        add=add+"_mc"
        ndim+=2

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

    rhodev2_array = np.zeros((1,rrange.size))

    if modded == False:
	data = np.loadtxt("%s/%s/chainstate_full.txt" % (output_dir, type_))
    else:
	data = np.loadtxt("%s/%s/chainstate_full_modded.txt" % (output_dir, type_))
	
    samples = np.asarray(data)
  
    """  
    length = data.size/nwalkers/ndim
    samples = data.reshape((length,ndim))
    """
    chainstate = samples[mini:maxi,:]
    print("Starting Calculations... "+str(chainstate.shape))

    for i in range(chainstate.shape[0]):
        print("Doing step %s out of %s" % (i, chainstate.shape[0]))
	pact = chainstate[i,0:8]
	rhodev2 = tools.logdev2_einasto_profile(rrange,pact)
	rhodev2_array = np.vstack((rhodev2_array,rhodev2))
    out = rhodev2_array[1:,:]
    np.savetxt("%s/%s/dev2_parts/MCMC_dev2_%s.txt" % (output_dir, type_, rank),out)
    print("Done!") 


if __name__ == "__main__":

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    mc = False
    rsteps = 25
    ndim = 8 #Number of parameters
    radmax = 10.0
    nwalkers = 28 #Total Number of walkers (98)
    steps = 2000000
    diag = False

    modded = True


    input_dir = "/work/dominik.zuercher/Output/splashpipe"
    output_dir = "/work/dominik.zuercher/Output/Mest"

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--type_", help="Type",default="Planck_PS_21")
    parser.add_argument("--add", help="Prior",default="")

    args=parser.parse_args()

    MCMC_calc_stats(ndim,rank,args.type_,args.add,mc)
