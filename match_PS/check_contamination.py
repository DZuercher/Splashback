import sys
sys.path.insert(0, "../../aum/install/lib64/python2.7/site-packages/")
import numpy as np
import cosmology as cc
import pandas
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,np.log10(8.0),1.0)

def get_app_magnitude(z):
    return -18.93629+25.0+5.0*np.log10(a.Dlofz(z))

def get_perturbed_magnitudes():
    # Read the matched file for SDSS z and PS magnitudes
    df = pandas.read_csv("/work/dominik.zuercher/Output/match_PS/matched_spec_new.dat", delim_whitespace=1, usecols=(0, 3, 4, 5, 6), header=None, names=(["zred", "rPS", "gPS", "iPS", "blue"]))
    
    # Read in Dominik's median errors as a function of g band magnitude and then perturb
    gmag, gmag_err = np.loadtxt("/work/dominik.zuercher/g_band_error.txt", unpack=1)
    rmag, rmag_err = np.loadtxt("/work/dominik.zuercher/r_band_error.txt", unpack=1)
    
    gspl = UnivariateSpline(gmag, np.log10(gmag_err))
    rspl = UnivariateSpline(rmag, np.log10(rmag_err))
    
    # Check the splines
    # xx = np.arange(15.0, 24.0, 0.01)
    # np.savetxt("test_gspl.dat", zip(xx, 10.0**gspl(xx)))
    # np.savetxt("test_rspl.dat", zip(xx, 10.0**rspl(xx)))
    
    
    df["gPS_pert"] = df.zred.values*0.0
    df["rPS_pert"] = df.zred.values*0.0
    import frogress
    for ii in frogress.bar(range(df.zred.size)):
        # First get the apparent magnitude at the redshift
        appmag_i = get_app_magnitude(df.zred.values[ii])
    
        appmag_g = df.gPS.values[ii] + (appmag_i-df.iPS.values[ii])
        appmag_r = df.rPS.values[ii] + (appmag_i-df.iPS.values[ii])
    
        # Now perturb according to the errors
        g_err = 10.**gspl(appmag_g)
        r_err = 10.**rspl(appmag_r)
    
        appmag_g += np.random.normal()*g_err
        appmag_r += np.random.normal()*r_err
    
        df.gPS_pert.values[ii] = appmag_g
        df.rPS_pert.values[ii] = appmag_r
    
    df.to_csv("Perturbed_match.dat", sep=" ", index=False)

def compute_contamination():
    df = pandas.read_csv("Perturbed_match.dat", delim_whitespace=1)

    zred, cut = np.loadtxt("/work/dominik.zuercher/Output/match_PS/spline_data.dat", unpack=1)
    spl = UnivariateSpline(zred, cut)

    bluepert = (df.gPS_pert.values-df.rPS_pert.values)<spl(df.zred.values)
    # for ii in range(df.gPS_pert.values.size):
    #     flag = (df.gPS_pert.values[ii]-df.rPS_pert.values[ii])<spl(df.zred.values[ii])
    #     print flag, df.gPS_pert.values[ii], df.rPS_pert.values[ii], spl(df.zred.values[ii])

    print np.sum( (df.blue.values<0.5) & (bluepert>0.5) )*1.0/np.sum(bluepert>0.5)

def get_red_fraction():
    df = pandas.read_csv("Perturbed_match.dat", delim_whitespace=1)
    
    zarr = np.arange(0.03, 0.33, 0.025)
    redfrac = np.zeros(zarr.size-1)
    for ii in range(redfrac.size):
        redfrac[ii] = np.sum( (df.zred.values>zarr[ii]) & (df.zred.values<zarr[ii+1]) & (df.blue.values<0.5) )*1.0/np.sum( (df.zred.values>zarr[ii]) & (df.zred.values<zarr[ii+1])  )
    np.savetxt("Redfraction_z.dat", np.transpose([zarr[:-1]/2+zarr[1:]/2, redfrac]))

get_perturbed_magnitudes()
compute_contamination()
get_red_fraction()
