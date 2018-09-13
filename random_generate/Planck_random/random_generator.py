#! coding=utf8
import numpy as np
import healpy
import pyfits
import matplotlib.pyplot as plt
from astropy import units as uwell
from astropy.coordinates import SkyCoord


if __name__=="__main__":

    seln_func = healpy.read_map("/work/dominik.zuercher/DataStore/Planck-SZ/HFI_PCCS_SZ-selfunc-union-cosmolog_R2.08.fits") #Planck survey selection function
    hdulist = pyfits.open("/work/dominik.zuercher/DataStore/Planck-SZ/HFI_PCCS_SZ-union_R2.08.fits") #Planck SZ catalog
    redshift_lim = 0.33 #Upper redshift limit for clusters
    Nsize = 10000 #Number of random clusters drawn
    sel_fig = "selectionfunction.pdf"
    bcg_fig = "onlyBCG.png"
    dist_fig = "dist.png"
    out_path = "/work/dominik.zuercher/DataStore/Planck-SZ/Planck_random_cosmo.csv"

    data = hdulist[1].data
    zred = data["redshift"]
    idx = (zred != -1.)
    zred = zred[idx]
    idy = (zred <= 0.33)
    zred = zred[idy]
    mu_r = -1.0 + 2.*np.random.random(size = Nsize)
    phi_r = 2.0*np.pi*np.random.random(size = Nsize)
    theta_r = np.arccos(mu_r)

    pix = healpy.pixelfunc.ang2pix(2048, theta_r, phi_r)
    mask = (seln_func[pix] == 1) #Select only points within the seln_func
    theta_r = theta_r[mask]
    phi_r = phi_r[mask]
    zred_r = zred[np.random.permutation(phi_r.size) % zred.size] #Produce 'random' redshifts

    fig = plt.figure(1)
    healpy.mollview(seln_func,norm = 'hist', cmap = 'PRGn')
    fig.savefig(sel_fig)

    fig = plt.figure(2)
    healpy.mollview(seln_func, norm = 'hist', cmap = 'PRGn')
    healpy.projscatter(planck_theta, planck_phi, s = 1, lw = 0)
    fig.savefig(bcg_fig)

    fig = plt.figure(3)
    healpy.mollview(seln_func, norm = 'hist', cmap = 'PRGn')
    healpy.projscatter(theta_r, phi_r, s = 1, lw = 0)
    fig.savefig(dist_fig)

    # Now convert using astropy coordinates
    c = SkyCoord(l = phi_r, b = np.pi/2 - theta_r, frame = 'galactic', unit = 'radian')
    output = np.transpose([c.icrs.ra.deg,c.icrs.dec.deg,zred_r])
    output = np.vstack((['RA','DEC','REDSHIFT'], output))
    np.savetxt(out_path, output, fmt = '%s', delimiter = ',')
