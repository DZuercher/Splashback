import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
import sys
import palettable
colors=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors


deproject=False

directory = sys.argv[1]

fname = directory + "/Clu=Norm_Gal=Norm/pairs.dat"
ranname = directory + "/Clu=Ran_Gal=Norm/pairs.dat"

rad, counts, avr, area = np.loadtxt(fname, usecols=(0, 1, 2, 3), unpack=1)
r_arr = np.unique(rad)
Nclu=np.loadtxt(directory+"/Clu=Norm_Gal=Norm/Ncluster.dat")
Nclu = np.tile(Nclu, r_arr.size)

if deproject==False:
    foo=np.loadtxt(ranname, unpack=1)
    rad, rancounts, avr, area, = np.loadtxt(ranname, usecols=(0, 1, 2, 3), unpack=1)
    Nclu_ran=np.loadtxt(directory+"/Clu=Ran_Gal=Norm/Ncluster.dat")
    Nclu_ran = np.tile(Nclu_ran, r_arr.size)

Nrad = r_arr.size
Njack = rad.size/Nrad
r_arr = rad[0::Njack]
xi = np.zeros(Nrad)
xi_err = np.zeros(Nrad)
xi_cov = np.zeros(Nrad*Nrad).reshape(Nrad, Nrad)
logrdiff=(r_arr[1]-r_arr[0])
r_arrbins = 10.0**np.append(r_arr, r_arr[-1]+logrdiff)
shell = 4./3.*np.pi*(r_arrbins[1:]**3-r_arrbins[:-1]**3)

if deproject:
    xijack = counts/np.repeat(shell,Njack)-1.
else:
    xijack = counts/(rancounts*Nclu/Nclu_ran) - 1.0

for i in range(r_arr.size):
    xi[i] = np.mean(xijack[Njack*i:Njack*(i+1)])

for i in range(r_arr.size):
    begi = Njack*i
    endi = Njack*(i+1)
    for j in range(r_arr.size):
        begj = Njack*j
        endj = Njack*(j+1)
        xi_cov[i][j] = np.mean((xijack[begi:endi]-xi[i])*(xijack[begj:endj]-xi[j]))*(Njack-1.)
    xi_err[i] = xi_cov[i][i]**0.5

ax = pl.subplot(221)
ax.errorbar(10.0**r_arr, xi, xi_err, fmt=".", color=colors[0])

from scipy.linalg import inv
print("SNR is: %s" % np.sqrt(np.dot(xi.T, np.dot(inv(xi_cov), xi))))

np.savetxt(directory + "/xi_2d.dat", np.transpose([10.0**r_arr, xi, xi_err]))
np.savetxt(directory + "/xi_2d_cov.dat", xi_cov)
ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_ylim(0.005, 1.0)
ax.set_xlim(0.1, 20.0)
pl.savefig(directory+"/plot.pdf")
