import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
from scipy.linalg import inv
import numpy as np
import pandas as pd

type_="Planck_PS_22"
est="DR"
#Reading the MLE values
mlep = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro/Mest/mles/MLE_parameters_"+str(type_)+"_"+str(est)+".txt", header=None)
mlep = np.asarray(mlep.values).reshape(mlep.size)
alpha_mle=mlep[1]
rs_mle=mlep[2]
rt_mle=mlep[5]
beta_mle=mlep[6]
gamma_mle=mlep[7]

snr_red=263.0
snr_sz=60.0 #SNR of Planck SZ
incr=snr_red/snr_sz

#get ranges for the new priors from the RedMaPPer errors
#The structure is new_error= RedMaPPer_error * SNR_correction

alpha_sigma=0.2*incr  # SZ_error= RedMaPPer_error/
rs_sigma=0.27*incr
rt_sigma=0.04*incr
beta_sigma=0.1*incr
gamma_sigma=0.14*incr

over=3.0

alpha_min=alpha_mle-over*alpha_sigma
alpha_max=alpha_mle+over*alpha_sigma
rs_min=rs_mle-over*rs_sigma
rs_max=rs_mle+over*rs_sigma

over_rt=over*6.0
rt_min=rt_mle-over_rt*rt_sigma
rt_max=rt_mle+over_rt*rt_sigma

beta_min=beta_mle-over*beta_sigma
beta_max=beta_mle+over*beta_sigma
gamma_min=gamma_mle-over*gamma_sigma
gamma_max=gamma_mle+over*gamma_sigma

print("Alpha flat prior from: "+str(alpha_min)+" to "+str(alpha_max))
print("rs flat prior from: "+str(rs_min)+" to "+str(rs_max))
print("rt flat prior from: "+str(rt_min)+" to "+str(rt_max))
print("Beta flat prior from: "+str(beta_min)+" to "+str(beta_max))
print("Gamma flat prior from: "+str(gamma_min)+" to "+str(gamma_max))


"""
prior="flat2"
p=0
rr=0
yy=0
yy_err=0
invcov=0
tools.get_chisquare(p,rr,yy,yy_err,invcov,prior,type_="Planck_PS_21",est="DR")


ptest=[-1.478296316147209577e+00,-6.477039135124851299e-01,3.338790721167303466e-01,1.856021404289147559e-05,-7.787123411120076000e-01,1.971757222801986620e+00,8.602966387981988983e-02,2.260197734649006307e+00]

rrange=np.linspace(0,10,100)

yy=tools.surf_den(rrange,ptest)
devs=tools.logdev_surf_den(rrange,ptest,recursive=True,surf=yy)
devs2=tools.logdev2_surf_den(rrange,ptest)
yyp=tools.einasto_profile(rrange,ptest)
devsp=tools.logdev_einasto_profile(rrange,ptest)
devs2p=tools.logdev2_einasto_profile(rrange,ptest)

rsp2d=tools.mindev2_procedure(ptest)
rsp3d=tools.mindev3_procedure(ptest)

f, axarr = plt.subplots(ncols=2,nrows=3)

axarr[0,0].plot(rrange,yy,label="Profile")
axarr[0,0].set_xscale("log")
axarr[0,0].set_yscale("log")
axarr[0,0].axvline(rsp2d)
axarr[1,0].plot(rrange,devs,label="Dev")
axarr[1,0].set_xscale("log")
axarr[1,0].axvline(rsp2d)
axarr[2,0].plot(rrange,devs2)
axarr[2,0].set_xscale("log")
axarr[2,0].axvline(rsp2d)
axarr[2,0].axhline(0)

axarr[0,1].plot(rrange,yyp,label="Profile")
axarr[0,1].set_xscale("log")
axarr[0,1].set_yscale("log")
axarr[0,1].axvline(rsp3d)
axarr[1,1].plot(rrange,devsp,label="Dev")
axarr[1,1].set_xscale("log")
axarr[1,1].axvline(rsp3d)
axarr[2,1].plot(rrange,devs2p)
axarr[2,1].set_xscale("log")
axarr[2,1].axvline(rsp3d)
axarr[2,1].axhline(0)

plt.savefig("Test.png")

"""
