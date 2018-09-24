#! coding=utf8
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/Splashback/toolbox')
import tools
import scipy.integrate as integrate
import scipy.special as special
from scipy.optimize import fmin
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import corner
import emcee
from scipy.linalg import inv


def MLE_procedure(type_, mc):
    input_dir = "/work/dominik.zuercher/Output/splashpipe/"+type_
    if mc == False:
        type_ += "_no_mc"
    output_dir = "/work/dominik.zuercher/Output/Mest/"+type_
    data = pd.read_csv("%s/xi_2d.dat" % (input_dir), header = None, sep = ' ')
    cov_data = pd.read_csv("%s/xi_2d_cov.dat" % (input_dir), header = None, sep = ' ')
    data = np.asarray(data.values)
    cov = np.asarray(cov_data.values)

    rr = data[:,0]
    yy = data[:,1]
    yy_err = data[:,2]
    idx = (data[:,0] > 0.0) & (data[:,0] < radmax) #radial cut
    rr = rr[idx]
    yy = yy[idx]
    yy_err = yy_err[idx]
    cov = np.transpose(cov[idx])[idx]
    invcov = inv(cov)
    if diag == True:
        invcov = np.diag(np.diag(invcov))
    """
    fig = plt.figure(1)
    plt.matshow(cov)
    plt.colorbar()
    plt.savefig("%s/%s_covariance.pdf" % (input_dir, type_))
    rij = np.zeros_like(cov)
    for i in range(cov.shape[0]):
        for j  in range(cov.shape[1]):
            rij[i,j] = cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])

    fig = plt.figure(2)
    plt.matshow(rij)
    plt.colorbar()
    plt.savefig("%s/%s_covmod.pdf" % (input_dir, type_))
    """
    if (type_=='Planck_PS_21.5_reloaded') & (mc==False):
        rhos=-1.
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.03
        se=1.2
        rt=0.4
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
    elif (type_=='Planck_PS_21_reloaded') & (mc==False):
        rhos=-0.8
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.02
        se=1.0
        rt=0.4
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
    elif (type_=='Planck_PS_19') & (mc==False):
        rhos=-0.8
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.02
        se=1.0
        rt=0.4
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
    elif (type_=='Planck_PS_22_reloaded') & (mc==False):
        rhos=-1.
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.02
        se=1.0
        rt=0.2
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
    elif (type_=='Planck_PS_21.5_red') & (mc==False):
        rhos=-0.4
        alpha=-1.1
        rs=0.03
        rho0=0.009
        se=0.84
        rt=-0.096
        beta=0.69
        gamma=0.13
    elif (type_=='Planck_PS_21.5_blue') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_low_2') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_low_3') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_low_4') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_low_5') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_low_6') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_cut_7') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_cut_8') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_blue_cut_9') & (mc==False):
        rhos=-1.0
        alpha=-0.67
        rs=0.35
        rho0=0.01
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_reloaded') & (mc==True):
        rhos=-0.8
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.04
        se=1.2
        rt=0.4
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
        fmin=0.1
        sigma=0.1
    elif (type_=='Planck_PS_21_reloaded') & (mc==True):
        rhos=-0.6
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.02
        se=1.0
        rt=0.4
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
        fmin=0.2
        sigma=0.4
    elif (type_=='Planck_PS_19') & (mc==True):
        rhos=-0.8
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.02
        se=1.0
        rt=0.4
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
        fmin=0.1
        sigma=0.1
    elif (type_=='Planck_PS_22_reloaded') & (mc==True):
        rhos=-0.8
        alpha=np.log10(0.2)
        rs=0.2
        rho0=0.025
        se=1.0
        rt=0.2
        beta=np.log10(4.0)
        gamma=np.log10(6.0)
        fmin=0.1
        sigma=0.1
    elif (type_=='Planck_PS_21.5_red_spline') & (mc==True):
        rhos=-0.4
        alpha=-1.1
        rs=0.03
        rho0=0.009
        se=0.84
        rt=-0.096
        beta=0.69
        gamma=0.13
        fmin=0.1
        sigma=0.1
    elif (type_=='Planck_PS_21.5_blue_spline') & (mc==True):
        rhos=-1.5
        alpha=-0.67
        rs=0.35
        rho0=0.001
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
        fmin=0.1
        sigma=0.1
    elif (type_=='Planck_PS_21.5_red_spline_no_mc') & (mc==False):
        rhos=-0.4
        alpha=-1.1
        rs=0.03
        rho0=0.009
        se=0.84
        rt=-0.096
        beta=0.69
        gamma=0.13
    elif (type_=='Planck_PS_21.5_blue_spline_no_mc') & (mc==False):
        rhos=-1.5
        alpha=-0.67
        rs=0.35
        rho0=0.001
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    elif (type_=='Planck_PS_21.5_red_hard_spline_no_mc') & (mc==False):
        rhos=-0.4
        alpha=-1.1
        rs=0.03
        rho0=0.009
        se=0.84
        rt=-0.096
        beta=0.69
        gamma=0.13
    elif (type_=='Planck_PS_21.5_blue_hard_spline_no_mc') & (mc==False):
        rhos=-1.5
        alpha=-0.67
        rs=0.35
        rho0=0.001
        se=0.58
        rt=0.144
        beta=0.57
        gamma=0.24
    else: 
        print("No prior defined!")

    if mc == False:
	p0 = np.array([rhos, alpha, rs, rho0, se, rt, beta, gamma]) #8 values
	yy_init = tools.surf_den(rr, p0)
	#yy_init_inner = tools.surf_den_inner(rr, p0)
	#yy_init_outer = tools.surf_den_outer(rr, p0,mc)
    else:
        p0 = np.array([rhos, alpha, rs, rho0, se, rt, beta, gamma,fmin,sigma]) #10 values
	yy_init = tools.total_surf_den(rr, p0)


    if mle == True:
        print "Performing MLE..."
        if prior_on == False:
	    if mc == False:
		p0 = np.loadtxt("%s/MLE_parameters.txt" % (output_dir))
	    else:
		p0 = np.loadtxt("%s/MLE_parameters_mc.txt" % (output_dir))
	print("The chi2 of the inital guess is: %s" % tools.get_chisquare_reload(p0, rr, yy, yy_err, invcov, prior, type_, mc))
	pMLE = tools.MLE_procedure_alt_reload(rr, yy, yy_err, invcov, p0, prior, type_, mc)

        if mc == False:
            np.savetxt("%s/MLE_parameters.txt" % (output_dir), pMLE)
        else:
            np.savetxt("%s/MLE_parameters_mc.txt" % (output_dir), pMLE)
	print("The chi2 of the MLE is: %s" % tools.get_chisquare_reload(pMLE, rr, yy, yy_err, invcov, prior,type_,mc))
	if mc == False:
	    yy_mle = tools.surf_den(rr, pMLE)
	    #yy_mle_inner = tools.surf_den_inner(rr, pMLE)
	    #yy_mle_outer = tools.surf_den_outer(rr, pMLE)
	else:
	    yy_mle = tools.total_surf_den(rr, pMLE)
	    #yy_mle_inner = tools.total_surf_den(rr, pMLE,part=-1)
	    #yy_mle_outer = tools.total_surf_den(rr, pMLE,part=1)


    fig = plt.figure(1)
    ax = plt.subplot(111)
    ax.errorbar(rr, yy, yy_err, fmt = ".")
    ax.plot(rr, yy_init, 'r.:', label = "INIT")
    #ax.plot(rr,yy_init_inner,'y.:',label="INIT INNER")
    #ax.plot(rr,yy_init_outer,'m.:',label="INIT OUTER")
    ax.plot(rr,yy_mle,'g.:',label="MLE")
    #ax.plot(rr,yy_mle_inner,color='yellow',linestyle='dashed',label="MLE INNER")
    #ax.plot(rr,yy_mle_outer,color='magenta',linestyle='dashed',label="MLE OUTER")

    ax.set_title("Estimation")
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    #ax.set_ylabel(r"$\xi^{\rm 2d}$ ($h^{-1}$Mpc)")    
    #ax.set_ylim([1e-4,1e1])
    plt.legend()
    plt.savefig("%s/Estimation.pdf" % (output_dir))
    plt.close()
    print("Plot done")




if __name__ == "__main__":
    prior_on = True
    mle = True #If true performs MLE search
    diag = False #If true reduces covariance matrix to diagonal entries
    radmax = 10.0 #MLE uses points out to radmax Mpc h^-1 
    prior = "best" #Which prior to use in get_chisquare
    mc = False #To use miscentering model or not


    
    types = ['Planck_PS_21.5_blue_hard_spline', 'Planck_PS_21.5_red_hard_spline']
    for type_ in types:
        MLE_procedure(type_,mc)
