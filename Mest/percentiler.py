from scipy.linalg import inv
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import pandas as pd
import numpy as np
from copy import deepcopy
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size


type='Planck_PS'
est='DR'


def procedure(rank):
	
	ndim=8
	steps=50000
	nwalkers=42
	runnum=50
	
	data = np.genfromtxt("/work/dominik.zuercher/MCMC_Store/chains/"+str(type)+"/chainstate_"+str(runnum)+".txt")
	samples = data.reshape((nwalkers,steps,ndim))

        chunky=int(steps/size)
        rest=steps-chunky*size
        mini=chunky*rank
        maxi=chunky*(rank+1)
        if rank>=(size-1)-rest:
                maxi+=1+rank-(size-1-rest)
                mini+=rank-(size-1-rest)
        if rank==size-1:
                maxi=steps-1
        mini=int(mini)
        maxi=int(maxi)
	
	data = pd.read_csv("/work/dominik.zuercher/DataStore/corr-pairs/"+str(type)+"/"+str(type)+"_plot("+str(est)+").dat", header=None, sep = ' ')
	cov_data = pd.read_csv("/work/dominik.zuercher/DataStore/corr-pairs/"+str(type)+"/"+str(type)+"_plot_cov("+str(est)+").dat", header=None, sep = ' ')
	data = np.asarray(data.values)
	cov = np.asarray(cov_data.values)
	cov=np.diag(np.diag(cov))
	rr = data[:,0]
	yy = data[:,1]
	yy_err = data[:,2]

	idx = (data[:,0] > 0.1) & (data[:,0] < 10.0)
	rr=rr[idx]
	yy=yy[idx]
	yy_err=yy_err[idx]
	cov=np.transpose(cov[idx])[idx]
	invcov=inv(cov)

	rrange=np.linspace(rr[0],rr[-1],200)
	sigma_array = np.zeros((1,rrange.size))
	rho_array = np.zeros((1,rrange.size))
	dev_array=np.zeros((1,rrange.size))
	rhodev_array=np.zeros((1,rrange.size))
	#dev2_array=np.zeros((1,rrange.size))
	#rhodev2_array=np.zeros((1,rrange.size))
	rsp2d_array=np.zeros((1,1))
	rsp3d_array=np.zeros((1,1))
	chisquare_array=np.zeros((1,1))
	chisquareconst_array=np.zeros((1,1))
	maxrhodev_array=np.zeros((1,1))
	maxinnerrhodev_array=np.zeros((1,1))

	#Calculate bands for modelplot
	it=mini
	print("Starting Kernel "+str(rank)+" with "+str(maxi-mini)+" objects")
	print("Calculating percentiles for plots...")
	while it<=maxi:
        	for j in range(nwalkers):        
			pact=samples[j,it,:]
			#For first plot
			sigma = tools.surf_den(rrange,pact)
			sigma_array = np.vstack((sigma_array,sigma))
			rho = tools.einasto_profile(rrange,pact)
			rho_array = np.vstack((rho_array,rho))
			#For second plot
			dev = tools.logdev_surf_den(rrange,pact,recursive=True,surf=sigma)
                	dev_array = np.vstack((dev_array,dev))
			rhodev=tools.logdev_einasto_profile(rrange,pact)
			rhodev_array=np.vstack((rhodev_array,rhodev))
                	"""
			#For third plot
			
			dev2 = tools.logdev2_surf_den(rrange,pact,recursive=True,surf=sigma,logdev=dev)
                	dev2_array = np.vstack((dev2_array,dev2))
			rhodev2 = tools.logdev2_einasto_profile(rrange,pact)
                	rhodev2_array = np.vstack((rhodev2_array,rhodev2))
                	"""
			#Additional values (splashbacks and chisquares)
			rsp2d = tools.mindev2_procedure(pact)
                	rsp2d_array=np.vstack((rsp2d_array,rsp2d))
			rsp3d = tools.mindev3_procedure(pact)
                	rsp3d_array=np.vstack((rsp3d_array,rsp3d))
			chisquare=tools.get_chisquare(pact,rr,yy,yy_err,invcov,prior="None")
			chisquare_array=np.vstack((chisquare_array,chisquare))
			chisquareconst=tools.get_chisquare(pact,rr,yy,yy_err,invcov,prior="None",constrained=True)
			chisquareconst_array=np.vstack((chisquareconst_array,chisquareconst))		
			maxrhodev=tools.logdev_surf_den(rsp3d,pact)
			maxinnerrhodev=tools.logdev_einasto_profile_inner(rsp3d,pact)
			maxrhodev_array=np.vstack((maxrhodev_array,maxrhodev))
			maxinnerrhodev_array=np.vstack((maxinnerrhodev_array,maxinnerrhodev))




		it+=1
	#dev2_array=dev2_array[1:,:]
	#rhodev2_array=rhodev2_array[1:,:]
	full_array=np.hstack((sigma_array,rho_array,dev_array,rhodev_array,rsp2d_array,rsp3d_array,chisquare_array,chisquareconst_array,maxrhodev_array,maxinnerrhodev_array))
	full_array=full_array[1:,:]	
	np.savetxt("/work/dominik.zuercher/MCMC_Store/MCMC_data/"+str(type)+"/parts/MCMC_data_"+str(type)+"_"+str(rank)+"_"+str(est)+".txt",full_array)
	
	print("Kernel "+str(rank)+" finished")


procedure(rank)
comm.Barrier()












