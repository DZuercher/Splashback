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
from subprocess import call



def MCMC_procedure(type_, add, prevrun, mc = False):
    if mc == True:
        add = add + "_mc"

    #Reading of the data
    print("Doing %s with prior %s" % (type_, add))
    data = pd.read_csv("%s/%s/xi_2d.dat" % (input_dir, type_), header = None, sep = ' ')
    cov_data = pd.read_csv("%s/%s/xi_2d_cov.dat" % (input_dir, type_), header = None, sep = ' ')
    if mc == False:
        mlep = pd.read_csv(output_dir + "/" + type_ +"/MLE_parameters.txt", header = None)
    else:
        mlep = pd.read_csv(output_dir + "/" + type_ + "/MLE_parameters_mc.txt", header = None)

    data = np.asarray(data.values)
    rr = data[:,0]
    yy = data[:,1]
    yy_err = data[:,2]
    idx = (data[:,0] > 0.1) & (data[:,0] < radmax)
    cov = np.asarray(cov_data.values)
    cov = np.transpose(cov[idx])[idx]
    rr = rr[idx]
    yy = yy[idx]
    yy_err = yy_err[idx]
    rrange = np.linspace(rr[0], rr[-1], rsteps)

    mlep = np.asarray(mlep.values).reshape(mlep.size)

    inv_cov = np.linalg.inv(cov)
    if diag == True:
        inv_cov = np.diag(np.diag(inv_cov))

    #Initialize parallization
    pool = MPIPool(loadbalance = True)

    #Check for previousruns and draw inital parameters around MLE 
    if pool.is_master():
        call("mkdir -p %s/%s/chainstates" % (output_dir, type_), shell = 1)
        call("mkdir -p %s/%s/pos" % (output_dir, type_), shell = 1)
        if prevrun != 0:
            p0 = np.loadtxt("%s/%s/pos/pos_%s%s_%s.dat" % (output_dir, type_, type_, add, prevrun))
        else:
            rhos_mle = mlep[0]
            alpha_mle = mlep[1]
            rs_mle = mlep[2]
            rho0_mle = mlep[3]
            se_mle = mlep[4]
            rt_mle = mlep[5]
            beta_mle = mlep[6]
            gamma_mle = mlep[7]
            if mc == True:
                fmin = mlep[8]
                sigma = mlep[9]
            
            snr_red = 263.0
            snr_sz = 60.0
            incr = snr_red/snr_sz
 
            rhos_sigma = 0.4/0.6*incr*abs(rhos_mle)
            alpha_sigma = 0.2/0.98*incr*abs(alpha_mle)
            rs_sigma = 0.27/0.32*incr*abs(rs_mle)
            rho0_sigma = 0.055/0.055*incr*abs(rho0_mle)
            se_sigma = 0.08/1.55*incr*abs(se_mle)
            rt_sigma = 0.04/0.087*incr*abs(rt_mle)
            beta_sigma = 0.1/1.0*incr*abs(beta_mle)
            gamma_sigma = 0.14/0.85*incr*abs(gamma_mle)

            from numpy.random import normal
            over = 0.4 #Inrease the standard deviation by 20% to achieve overdisperse sampling
            print("---------------------------------------------------------------------------")
            print("Initial standard deviation for draw:")
            print(rhos_sigma*over,alpha_sigma*over,rs_sigma*over,rho0_sigma*over,se_sigma*over,rt_sigma*over,beta_sigma*over,gamma_sigma*over)
            print("---------------------------------------------------------------------------")
            rhos_draws = normal(loc=rhos_mle,scale=rhos_sigma*over,size=nwalkers)
            alpha_draws = normal(loc=alpha_mle,scale=alpha_sigma*over,size=nwalkers)
            rs_draws = normal(loc=rs_mle,scale=rs_sigma*over,size=nwalkers)
            rho0_draws = normal(loc=rho0_mle,scale=rho0_sigma*over,size=nwalkers)
            #Tight initial draw of se
            se_draws = normal(loc=se_mle,scale=0.01,size=nwalkers)
            rt_draws = normal(loc=rt_mle,scale=rt_sigma*over,size=nwalkers)
            beta_draws = normal(loc=beta_mle,scale=beta_sigma*over,size=nwalkers)
            gamma_draws = normal(loc=gamma_mle,scale=gamma_sigma*over,size=nwalkers)

            p0 = [[rhos_draws[i],alpha_draws[i],rs_draws[i],rho0_draws[i],se_draws[i],rt_draws[i],beta_draws[i],gamma_draws[i]] for i in range(nwalkers)]

            if init_type == "tight":
                if mc == False:
                    p0 = [mlep + 1e-3*np.random.randn(ndim) for i in range(nwalkers)] #Draw from small ball around MLE
                else:
                    p0 = [mlep + 1e-3*np.random.randn(ndim+2) for i in range(nwalkers)] #Draw from small ball around MLE
    else:
        pool.wait()
        sys.exit(0)

    #Initialize sampler
    if mc == False:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, tools.lnlike_reload, args=[rr, yy, yy_err, inv_cov, prior,type_,mc], pool=pool)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim+2, tools.lnlike_reload, args=[rr, yy, yy_err, inv_cov, prior,type_,mc], pool=pool)

    print("Initialized")
    if prevrun==0:
        print("Performing Burn-in Phase")
        pos, prob, state = sampler.run_mcmc(p0, burn_in)
        sampler.reset()
        print("Burn-in Phase done")
    else:
        pos = p0

    ntot=0
    it = prevrun
    print("Starting MCMC...")
    while ntot<steps:
        it+=1
        chainstate=[]
        ii=0
        while ii<nsave:
            for result in sampler.sample(pos,iterations=1,storechain=False):
                pos=result[0]
                chainstate.append(pos)
                ii+=1
        chainstate=np.asarray(chainstate)
        newpos=chainstate[-1]
        ntot += nsave	

        np.savetxt("%s/%s/pos/pos_%s%s_%s.dat" % (output_dir, type_, type_, add, it),newpos)

        with file("%s/%s/chainstates/chainstate_%s%s_%s.txt" % (output_dir, type_, type_, add, it), 'w+') as outfile:
            outfile.write('# Array shape: {0}\n'.format(chainstate.shape))
            for data_slice in chainstate:
                np.savetxt(outfile, data_slice, fmt='%-7.9f')
                outfile.write('# New slice\n')


        print("Phase "+str(it)+" out of "+str(steps/nsave)+" phases done with a mean acceptance fraction of: "+str(np.mean(sampler.acceptance_fraction)))
        pos=newpos
    pool.close()



if __name__ == "__main__":

    mc = True

    init_type = "tight"
    diag = False #Reduces covariance matrix to diagonals
    rsteps = 25 #Iterations in each MCMC chain
    ndim = 8 #Number of parameters
    radmax = 10.0 #Maximal radius to consider
    nwalkers = 28 #Total Number of walkers (98)
    burn_in = 10000 #Number of MCMC burn in points used per walker (1000)
    prevrun = 0 #Number of the last runfile
    nsave = 1000 #steps in between savepoints (1000)
    prior = "best"
    input_dir = "/work/dominik.zuercher/Output/splashpipe"
    output_dir = "/work/dominik.zuercher/Output/Mest"

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--type_", help = "Type", default = "Planck_PS_21")
    parser.add_argument("--add", help = "Prior", default = "")
    parser.add_argument("--prevrun", help = "Previous runs", default = 0)
    parser.add_argument("--size", help = "Number of processors", default = 28)
    parser.add_argument("--steps", help = "Number of iterations", default = 100000)

    args = parser.parse_args()
    size = int(args.size)
    steps = int(args.steps)

    MCMC_procedure(args.type_, args.add, int(args.prevrun), mc)
