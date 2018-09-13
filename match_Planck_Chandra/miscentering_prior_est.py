# Uses Chandra and Planck clusters to find prior values for miscentering parameters fmin and sigma
import scipy.optimize as op
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cosmo = FlatLambdaCDM(H0 = 70 * u.km / u.s / u.Mpc, Om0 = 0.3)

def comoving_separation(z,theta):
    theta = theta/3600
    theta = theta*np.pi/180.0
    dc = cosmo.comoving_distance(z) #comoving distance
    da = dc/(1.0+z) #angular diameter distance
    rp = theta*da #physical separation
    rc = rp*(1.0+z) 
    rc = rc*0.7 #Output in comoving h**-1 Mpc
    return rc

def rayleigh_distribution(R,sigma):
    return R/sigma**2*np.exp(-R**2/(2*sigma**2))

def lnlike(sigma,R,yy):
    chisq = np.sum((yy - rayleigh_distribution(R,sigma))**2)
    return -0.5*chisq

def MLE_procedure(RR, yy, p0):
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, p0, args=(RR, yy))
    return result["x"]



if __name__=="__main__":

    data = pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/30/Planck_corrected_matched.csv")
    redshifts = data["REDSHIFT"].values
    distances = data["Chandra_distance"].values
    idx = distances > -10.0
    distances = distances[idx]
    redshifts = redshifts[idx]
   
    total=redshifts.size

    lowcut = 0.25

    com_dist = comoving_separation(redshifts,distances).value
    idx = com_dist > lowcut
    com_dist = com_dist[idx]

    reduced=com_dist.size
    fmin_est = float(reduced)/float(total)


    hist, bin_edges = np.histogram(com_dist,density=True,bins=20)
    bin_centers = (bin_edges[1:] - bin_edges[:-1])/2.0 + bin_edges[:-1]
    idx = hist > 0
    hist = hist[idx]
    bin_centers = bin_centers[idx]

    sigma_est = MLE_procedure(bin_centers, hist, 0.4)

    Rrange = np.linspace(np.min(bin_centers),np.max(bin_centers),100)
    modely = rayleigh_distribution(Rrange,sigma_est)

    plt.plot(Rrange,modely,'r--')
    plt.bar(bin_centers,hist/reduced,width=(bin_edges[1]-bin_edges[0]))
    plt.xlabel("R h**-1 Mpc (comoving)")
    plt.savefig("distance_fit.pdf") 

    print("fmin : "+str(fmin_est))
    print("sigma : "+str(sigma_est))




