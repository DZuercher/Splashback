from scipy.interpolate import UnivariateSpline
from  scipy.optimize import curve_fit
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import frogress

def cut(red, ax = 0.65):
    return ax + ((red - 0.1)*0.8/0.23)


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
            red_sigma = popt[2]
        else:
            red_mu = popt[4]
            red_sigma = popt[5]

        if include_PS_errors == False:
            cut = red_mu - confidence * np.abs(red_sigma)
        else:
	    cut = red_mu - confidence * np.sqrt( red_sigma**2.0  + red_err**2.0 + green_err**2.0 )

        cuts[it - 1] = cut 

        if (it % 10 == 0):
            plt.figure(it / 10)
            plt.hist(colors, bins = bins)
            gr_range = edges[:-1] + (edges[1] - edges[0])/2.
            plt.plot(gr_range, double_gaussian(gr_range, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]), 'k-', lw = 0.3)
            plt.axvline(cut, color = 'r', lw = 0.3)
  	    plt.savefig("%s/test_hist_%s.pdf" % (output_dir, it/ 10 ))
    return cuts


if __name__ == "__main__":

    min_z = 0.03
    max_z = 0.33
    z_step = 100
    bins = 500

    confidence = 3 # Confidence level to exclude Reds from Blues (in sigmas)

    catalog = "/work/dominik.zuercher/Output/match_PS_GAMA/matched_spec_new.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS_GAMA"


    data = np.loadtxt(catalog)
    redshift = data[:,0]
    red_PS = data[:,3]
    green_PS = data[:,4]
    iband_PS = data[:,5]
    id_ = data[:,6]

    #Remove Nans
    idx = (np.logical_not(np.isnan(data[:,1]))) & (np.logical_not(np.isnan(data[:,2])))
    redshift = redshift[idx]
    red_PS = red_PS[idx]
    green_PS = green_PS[idx]
    id_ = id_[idx]

    idx = (red_PS < 900) & (red_PS > -900) & (green_PS < 900) & (green_PS > -900) 
    redshift = redshift[idx]
    red_PS = red_PS[idx]
    green_PS = green_PS[idx]
    id_ = id_[idx]

    color = green_PS - red_PS
    print("Catalog read")


    bin_num = bin_redshifts(redshift, min_z, max_z, z_step)
    print("Redshift binning done")

    cuts = get_cuts(color, bin_num, bins, z_step, confidence, False)
    print("histograms made")

    """
    idx = np.random.rand(redshift.size)
    boolar = idx<0.05
    redshift = redshift[boolar]
    color = color[boolar]
    id_ = id_[boolar]
    """   

    ax = plt.subplot(221)
    ax.scatter(redshift[id_ == True], color[id_ == True], marker='.', s = 0.01, color = 'b', label = "Blue")
    ax.scatter(redshift[id_ == False], color[id_ == False], marker='.', s = 0.01, color = 'r', label = "Red")
    z_edges = np.linspace(min_z, max_z, num = z_step)
    z_middles = z_edges[:-1] + (z_edges[1] - z_edges[0])/2.

    idx = (cuts > 0.0) & (cuts < 1.5)
    z_middles = z_middles[idx]
    cuts = cuts[idx]

    spl = UnivariateSpline(z_middles, cuts)

    rrange = np.linspace(0, 0.35, 1000)

    output = np.hstack((z_middles.reshape(z_middles.size,1), cuts.reshape(cuts.size,1)))
    np.savetxt("%s/spline_data.dat" % output_dir, output)


    #Using the SDSS spline that was also used in the analysis
    foo = np.genfromtxt("/work/dominik.zuercher/Output/match_PS/spline_data.dat", unpack=1)
    spl = UnivariateSpline(foo[0,:], foo[1,:])
    ax.plot(rrange, spl(rrange), 'k-', lw=0.7, label = "spline")

    ax.set_ylim([0,2])
    ax.set_xlim([0.03,0.33])
    ax.set_xlabel("z")
    ax.set_ylabel(r"g$_{\mathrm{P1}}$ - r$_{\mathrm{P1}}$")
    plt.savefig("%s/color_vs_z.pdf" % output_dir)


