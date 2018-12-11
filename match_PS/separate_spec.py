import frogress
import pandas as pd
import sys
sys.path.insert(0, "../../aum/install/lib64/python2.7/site-packages/")
import numpy as np
import cosmology as cc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from  scipy.optimize import curve_fit

def straight(x, a, m):
    return a + m*x

def double_gaussian(x, c1, mu1, sigma1, c2, mu2, sigma2):
    res =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) + c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )
    return res


def bin_redshifts(redshift, min_z, max_z, z_steps):
    z_edges = np.linspace(min_z, max_z, num = z_steps)
    bin_num = np.zeros(redshift.size)
    for ii,z in frogress.bar(enumerate(redshift)):
        it = 0
        upper_z = z_edges[it]
        while z > upper_z:
           it += 1
           upper_z = z_edges[it]
        bin_num[ii] = it
    return bin_num

#Not the same as in fitting.py !!!!
def get_cuts(color, bin_num, bins, z_steps, confidence, include_PS_errors):
    cuts = np.zeros(z_steps-1)
    for it in frogress.bar(range(1,z_steps)):
        colors = color[bin_num == it]
        hist, edges = np.histogram(colors, bins = bins)
        try:
            popt, pcov = curve_fit(double_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist)
        except:
            cuts[it-1] = -1000
            continue

        red_ind = np.argmax([popt[1], popt[4]])
        if red_ind == 0:
            red_mu = popt[1]
            blue_mu = popt[4]
            red_sigma = popt[2]
        else:
            red_mu = popt[4]
            blue_mu = popt[1]
            red_sigma = popt[5]

	xx = edges[:-1] + (edges[1] - edges[0])/2.

        res = double_gaussian(xx, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
        cut = xx[np.argmin(res[(xx<red_mu) & (xx>blue_mu)])]

        cuts[it - 1] = cut
        """
	plt.figure(it)
	plt.hist(colors, bins = bins)
	gr_range = edges[:-1] + (edges[1] - edges[0])/2.
	plt.plot(gr_range, double_gaussian(gr_range, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]), 'k-', lw = 0.3)
	plt.axvline(cut, color = 'r', lw = 0.3)
	plt.savefig("%s/test_hist_%s.pdf" % ("/work/dominik.zuercher/Output/match_PS_GAMA", it ))
        """
    return cuts




if __name__ == "__main__":

    a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,np.log10(8.0),1.0)

    dir_ = "/work/dominik.zuercher/Output/match_PS"

    spec_cat = "%s/matched_spec.dat" % dir_
    data = pd.read_csv(spec_cat, delim_whitespace=1,header=None,names=["zred","red","green","red_PS","green_PS","iband_PS"])
    cat = data.values

    """
    id_ = cat[:,0] <= 0.2
    cat = cat[id_,:]
    """

    Dl = np.zeros_like(cat[:,0])
    for ii,z in enumerate(cat[:,0]):
        Dl[ii] = np.log10(a.Dlofz(z))   

    #Convert to absoulute mag and k correct
    red = (cat[:,3] - cat[:,1]) - 25.0 - 5.0*Dl
    green = (cat[:,4] - cat[:,2]) - 25.0 - 5.0*Dl

    id_ = (green - red < 100) &( green - red > -100) & (red > -100)
    red = red[id_]
    green = green[id_]
    cat = cat[id_,:]

    #Custom color cut
    bin_num = bin_redshifts(red, -25., -10., 20)
    cuts = get_cuts(green-red, bin_num, 500, 20, 3, False)

    z_edges = np.linspace(-25., -10., num = 20)
    z_middles = z_edges[:-1] + (z_edges[1] - z_edges[0])/2.

    idx = (cuts > 0.0) & (cuts < 2.0)
    z_middles = z_middles[idx]
    cuts = cuts[idx]
    #spl = UnivariateSpline(z_middles, cuts)
    #popt = curve_fit(straight, z_middles, cuts)[0]
    #spl = lambda x: popt[0] + popt[1]*x


    idx = np.random.rand(green.size)
    boolar = idx<0.1
    pgreen = green[boolar]
    pred = red[boolar]

    xx = np.linspace(-30,0, 1000)
    #plt.plot(xx, 0.6 - 0.03*(xx + 20.0),'k-', lw=0.3)
    plt.plot(xx, 0.7 - 0.03*(xx + 20.0),'k--', lw=0.3)
    #plt.plot(xx, spl(xx), 'g--', lw=0.3)
    plt.scatter(pred, pgreen-pred, s = 0.01, alpha=0.5)
    plt.ylabel("g-r")
    plt.xlabel("r")
    plt.xlim([-22.5,-17])
    plt.ylim([0.0, 1.5])
    plt.savefig("%s/SDSS_test.pdf" % dir_)

    #using custom double goaussin fitted spline
    id_ = green - red <= 0.7 - 0.03*(red + 20.0)

    #Using SDSS directly
    #id_ = (green - red <= 0.8 - 0.03*(red + 20.0)) #True if identified as blue and False if red
    #Using a modified SDSS cut (shifted)
    #id_ = (green - red <= 0.6 - 0.03*(red + 20.0)) #True if identified as blue and False if red
    cat = np.hstack((cat,id_.reshape(id_.size,1))) 

    id_ = np.zeros(cat.shape[0])
    id_ = (cat[:,4] - cat[:,3] <= 0.8 - 0.03*(cat[:,3] + 20.0)) #True if identified as blue and False if red
    cat = np.hstack((cat,id_.reshape(id_.size,1))) 
    np.savetxt("%s/matched_spec_new.dat" % dir_,cat)

