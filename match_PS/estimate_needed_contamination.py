from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/Splashback/toolbox')
import tools

def NFW(r, rho0, rs):
    return rho0/(r/rs*(1 + r/rs)**2.0)

def einasto(r, rho0, A, alpha):
    return rho0*np.exp(-A*r**alpha)

contamination = 0.20
rsteps = 25
radmax = 10.0

red_params = np.genfromtxt("/work/dominik.zuercher/Output/Mest/Planck_PS_21.5_red_hard_spline_no_mc/results.csv", skip_header = 1, delimiter = ',')
red_params = red_params[0,:]
alpha = red_params[1]
beta = red_params[2]
gamma = red_params[3]
rho0 = red_params[4]
rhos = red_params[5]
rs = red_params[6]
rt = red_params[7]
se = red_params[8]

red_par_array = np.asarray([rhos, alpha, rs, rho0, se, rt, beta, gamma])


input_dir = "/work/dominik.zuercher/Output/splashpipe"
type_ = "Planck_PS_21.5_blue_hard_spline"
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
#weights = np.asarray([1.0 - 1.0/25.0*n for n in range(25)])


blue_array = np.genfromtxt("/work/dominik.zuercher/Output/Mest/Planck_PS_21.5_blue_hard_spline_no_mc/plotsave.dat")
blue_corr = blue_array[rsteps:rsteps*2,0]
blue_true_dev = blue_array[rsteps*3:rsteps*4,0]
#popt, pcov = curve_fit(einasto, rrange, blue_corr)
#blue_rho = einasto(rrange, 10.0, 0.5, 4.0)
#blue_dev = np.gradient(blue_rho, rrange)


#p = np.polyfit(rrange, np.log(blue_corr), deg=8)

#blue_rho = np.poly1d(p)
#blue_dev = blue_rho.deriv()

scale = 0.055

blue_rho = scale*rrange**(-1.5)
blue_dev = rrange*(-1.5*scale*rrange**(-2.5))/blue_rho


red_rho = tools.einasto_profile(rrange,red_par_array)
red_dev = tools.logdev_einasto_profile(rrange,red_par_array)
#blue_dev_proto = np.gradient(np.exp(blue_rho(rrange)),rrange)
#blue_dev = rrange/np.exp(blue_rho(rrange))*blue_dev_proto

full = contamination*red_rho + (1.0 - contamination)*blue_rho
full_dev = rrange/full*(contamination*tools.dev_einasto_profile(rrange,red_par_array) + (1.0 - contamination)*scale*-1.5*rrange**(-2.5))



f, ax = plt.subplots(ncols = 2, nrows = 1)
ax[0].plot(rrange, red_rho,'r-', label = r"Red 21.5 (6 $\sigma$)")
ax[0].plot(rrange, full,'k-', label = "%f contamination by Reds" % contamination)
ax[0].plot(rrange, blue_rho,'b-', label = r"$r^{-1.5}$ power law")
ax[0].plot(rrange, blue_corr,'b--', label = r"Blue 21.5 (6 $\sigma$)")
ax[1].plot(rrange, red_dev,'r-')
ax[1].plot(rrange, blue_dev,'b-')
ax[1].plot(rrange, full_dev,'k-')
ax[1].plot(rrange, blue_true_dev,'b--')
ax[0].set_xscale('log')
ax[1].set_xscale('log')
#ax[0].set_ylim([0,2])
ax[0].set_yscale('log')
#ax[1].set_xlim([1,6])
ax[0].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
ax[1].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
ax[0].set_ylabel(r"$\xi_{\rm 3D}(r)$")
f.legend()
f.savefig("/work/dominik.zuercher/Output/match_PS/add_corr.pdf")




