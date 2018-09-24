import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/Splashback/toolbox')
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


def vis_procedure(ndim,type_,add,modded=False,mc=False):
    if mc==True:
        ndim+=2
    #Reading data
    dat=pd.read_csv("%s/%s/xi_2d.dat" % (input_dir, type_), header=None,sep=' ')

    if mc == False:
	type_ += '_no_mc'


    dat = np.asarray(dat.values)
    rr = dat[:,0]
    yy = dat[:,1]
    yy_err = dat[:,2]

    idx = (dat[:,0] > 0.1) & (dat[:,0] < 10.0)
    rr=rr[idx]
    yy=yy[idx]
    yy_err=yy_err[idx]
    rrange=np.linspace(rr[0],rr[-1],rsteps)

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

    #Plotting graphs
    f, axarr = plt.subplots(ncols=2,nrows=2)
    axarr[0,0].errorbar(rr, yy, yy_err, fmt=".")
    axarr[0,0].plot(rrange,sigmas[:,0],'r-',label="MCMC fit")
    axarr[0,0].plot(rrange,sigmas[:,0]+sigmas[:,1],'r--')
    axarr[0,0].plot(rrange,sigmas[:,0]-sigmas[:,2],'r--')
    axarr[0,0].axvline(rsp2ds[0],c='k',label="RSP 2D")
    axarr[0,0].axvline(rsp2ds[0]+rsp2ds[1],c='k',ls='dashed')
    axarr[0,0].axvline(rsp2ds[0]-rsp2ds[2],c='k',ls='dashed')
    axarr[0,0].axvline(rsp3ds[0],c='b',label="RSP 3D",alpha=0.5)
    axarr[0,0].axvline(rsp3ds[0]+rsp3ds[1],c='b',ls='dashed',alpha=0.5)
    axarr[0,0].axvline(rsp3ds[0]-rsp3ds[2],c='b',ls='dashed',alpha=0.5)
    axarr[0,0].set_title("2D correlation signal")
    axarr[0,0].set_xscale("log")
    axarr[0,0].set_yscale("log")
    axarr[0,0].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    axarr[0,0].set_ylabel(r"$\xi^{\rm 2d}$ ($h^{-1}$Mpc)")
    axarr[0,0].legend(loc='lower left',fontsize="x-small")

    axarr[0,1].plot(rrange,rhos[:,0],'r-',label="MCMC fit")
    axarr[0,1].plot(rrange,rhos[:,0]+rhos[:,1],'r--')
    axarr[0,1].plot(rrange,rhos[:,0]-rhos[:,2],'r--')
    axarr[0,1].axvline(rsp2ds[0],c='k',label="RSP 2D",alpha=0.5)
    axarr[0,1].axvline(rsp2ds[0]+rsp2ds[1],c='k',ls='dashed',alpha=0.5)
    axarr[0,1].axvline(rsp2ds[0]-rsp2ds[2],c='k',ls='dashed',alpha=0.5)
    axarr[0,1].axvline(rsp3ds[0],c='b',label="RSP 3D")
    axarr[0,1].axvline(rsp3ds[0]+rsp3ds[1],c='b',ls='dashed')
    axarr[0,1].axvline(rsp3ds[0]-rsp3ds[2],c='b',ls='dashed')
    axarr[0,1].set_title("3D correlation signal")
    axarr[0,1].set_xscale("log")
    axarr[0,1].set_yscale("log")
    axarr[0,1].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    axarr[0,1].set_ylabel(r"$\xi^{\rm 2d}$ ($h^{-1}$Mpc)")
    #axarr[0,1].legend()

    axarr[1,0].plot(rrange,devs[:,0],'r-',label="MCMC derivative")
    axarr[1,0].plot(rrange,devs[:,0]+devs[:,1],'r--')
    axarr[1,0].plot(rrange,devs[:,0]-devs[:,2],'r--')
    axarr[1,0].axvline(rsp2ds[0],c='k',label="RSP 2D")
    axarr[1,0].axvline(rsp2ds[0]+rsp2ds[1],c='k',ls='dashed')
    axarr[1,0].axvline(rsp2ds[0]-rsp2ds[2],c='k',ls='dashed')
    axarr[1,0].axvline(rsp3ds[0],c='b',label="RSP 3D",alpha=0.5)
    axarr[1,0].axvline(rsp3ds[0]+rsp3ds[1],c='b',ls='dashed',alpha=0.5)
    axarr[1,0].axvline(rsp3ds[0]-rsp3ds[2],c='b',ls='dashed',alpha=0.5)
    axarr[1,0].set_title("2D First Derivative")
    axarr[1,0].set_xscale("log")
    axarr[1,0].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    axarr[1,0].set_ylabel(r"$\mathrm{d} log\xi^{\rm 2d} /\mathrm{d} logr$")
    #axarr[1,0].legend()


    axarr[1,1].plot(rrange,rhodevs[:,0],'r-',label="MCMC derivative")
    axarr[1,1].plot(rrange,rhodevs[:,0]+rhodevs[:,1],'r--')
    axarr[1,1].plot(rrange,rhodevs[:,0]-rhodevs[:,2],'r--')
    axarr[1,1].axvline(rsp2ds[0],c='k',label="RSP 2D",alpha=0.5)
    axarr[1,1].axvline(rsp2ds[0]+rsp2ds[1],c='k',ls='dashed',alpha=0.5)
    axarr[1,1].axvline(rsp2ds[0]-rsp2ds[2],c='k',ls='dashed',alpha=0.5)
    axarr[1,1].axvline(rsp3ds[0],c='b',label="RSP 3D")
    axarr[1,1].axvline(rsp3ds[0]+rsp3ds[1],c='b',ls='dashed')
    axarr[1,1].axvline(rsp3ds[0]-rsp3ds[2],c='b',ls='dashed')
    axarr[1,1].set_title("3D First Derivative")
    axarr[1,1].set_xscale("log")
    axarr[1,1].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    axarr[1,1].set_ylabel(r"$\mathrm{d} log\rho /\mathrm{d} logr$")
    #axarr[1,1].legend()
    plt.tight_layout()
    plt.savefig("%s/%s/model_plot_%s.pdf" % (output_dir, type_, type_))

    print("Plot done")
    
    #Drawing chains
    nochains=False
    if modded==False:
        data = np.loadtxt("%s/%s/chainstate_full.txt" % (output_dir, type_))
    else:
        data = pd.read_csv("%s/%s/chainstate_full_modded.txt" % (output_dir, type_),sep=' ',header=None,error_bad_lines=False)
        nwalkers=nwalkers_modded
    
    data=np.asarray(data)
    
    length=data.size/nwalkers/ndim
    try:
        samples = data.reshape((length,nwalkers,ndim))
    except:
        nochains=True
    if mc==False:
        param_names = [r'$\log_{10}(\rho_{\mathrm{s}})$',r'$\alpha$',r'$\log_{10}(r_{\mathrm{s}})$',r'$\log_{10}(\rho_{\mathrm{0}})$',r'$s_{\mathrm{e}}$',r'$\log_{10}(r_{\mathrm{t}})$',r'$\log_{10}(\beta)$',r'$\log_{10}(\gamma)$']
    else:
        param_names = [r'$\log_{10}(\rho_{\mathrm{s}})$',r'$\alpha$',r'$\log_{10}(r_{\mathrm{s}})$',r'$\log_{10}(\rho_{\mathrm{0}})$',r'$s_{\mathrm{e}}$',r'$\log_{10}(r_{\mathrm{t}})$',r'$\log_{10}(\beta)$',r'$\log_{10}(\gamma)$',r'f$_{\rm min}$',r'$\sigma$']

    if nochains==False:
	#Visualisation of chains
	fig=plt.figure(1)
	colors=np.linspace(0,1,int(nwalkers))
	f, axarr = plt.subplots(ndim,sharex=True)
	for i in range(ndim):
	    for j in range(nwalkers):
		x=np.linspace(1,samples.shape[0],samples.shape[0])
		axarr[i].plot(x,samples[:,j,i],c=cm.nipy_spectral(colors[j]),alpha=0.4)
		axarr[i].set_ylabel(param_names[i])

	plt.savefig("%s/%s/chains.pdf" % (output_dir, type_))
	plt.close(fig)
    print("chains drawn")
    

    #Printing value output
    flatsamples=data
    #Printing results to csv. Each column corresponds to one parameter. First value is 50% value, second is 84%-50% (upper std), third is 50%-16% (lower std)
    mcmc_par = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(flatsamples, [16, 50, 84], axis=0)))

    if mc==False:
        lrhos_mcmc, lalpha_mcmc, lrs_mcmc, lrho0_mcmc, lse_mcmc, lrt_mcmc, lbeta_mcmc, lgamma_mcmc = mcmc_par 
    else:
        lrhos_mcmc, lalpha_mcmc, lrs_mcmc, lrho0_mcmc, lse_mcmc, lrt_mcmc, lbeta_mcmc, lgamma_mcmc, fmin_mcmc, sigma_mcmc = mcmc_par 


    mcmc_par = np.asarray(mcmc_par)
    if mc==False:
        d = {'lrhos' : lrhos_mcmc, 'lalpha' : lalpha_mcmc, 'lrs' : lrs_mcmc, 'lrho0' : lrho0_mcmc, 'lse' : lse_mcmc, 'lrt' : lrt_mcmc, 'lbeta' : lbeta_mcmc, 'lgamma' : lgamma_mcmc, 'maxrhodev' : maxrhodevs,'maxinnerrhodev' : maxinnerrhodevs,'chisquare' :  chisquares, 'rsp2d' : rsp2ds, 'rsp3d' : rsp3ds}
    else:
        d = {'lrhos' : lrhos_mcmc, 'lalpha' : lalpha_mcmc, 'lrs' : lrs_mcmc, 'lrho0' : lrho0_mcmc, 'lse' : lse_mcmc, 'lrt' : lrt_mcmc, 'lbeta' : lbeta_mcmc, 'lgamma' : lgamma_mcmc, 'maxrhodev' : maxrhodevs,'maxinnerrhodev' : maxinnerrhodevs,'chisquare' :  chisquares, 'rsp2d' : rsp2ds, 'rsp3d' : rsp3ds, 'fmin': fmin_mcmc, 'sigma': sigma_mcmc}



    df = pd.DataFrame(data=d)
    df.to_csv("%s/%s/results.csv" % (output_dir, type_), sep=',', index=False)
    print("Results printed")


    #Corner Plot
    flatsamples=flatsamples.reshape([flatsamples.shape[0],flatsamples.shape[1]])
    fig = corner.corner(flatsamples,bins=50,labels=param_names,label_kwargs={'fontsize':16})
    fig.savefig("%s/%s/cornerplot.pdf" % (output_dir, type_))
    print("Cornerplots drawn")


if __name__ == "__main__":

    input_dir = "/work/dominik.zuercher/Output/splashpipe"
    output_dir = "/work/dominik.zuercher/Output/Mest"

    modded = True
    mc = False

    ndim = 8
    nwalkers = 28
    nwalkers_modded = 28
    steps = 100000
    rsteps = 25
    types = ['Planck_PS_21.5_blue_hard_spline','Planck_PS_21.5_red_hard_spline']
    adds = ["_best"]
    for add in adds:
        for type_ in types:
            vis_procedure(ndim,type_,add,modded,mc)
