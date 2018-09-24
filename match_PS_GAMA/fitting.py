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


def get_cuts(color, bin_num, bins, z_steps, confidence):
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
        cut = red_mu - confidence * np.abs(red_sigma)
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


    catalog = "/work/dominik.zuercher/Output/match_PS_GAMA/matched_spec.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS_GAMA"

    rrange = np.linspace(0, 0.35, 1000)

    data = np.loadtxt(catalog)
    redshift = data[:,0]
    color_PS = data[:,1]

    idx = (color_PS < 900) & (color_PS > -900) 
    redshift = redshift[idx]
    color_PS = color_PS[idx]
    print("Catalog read")

    bin_num = bin_redshifts(redshift, min_z, max_z, z_step)
    print("Redshift binning done")
    cuts = get_cuts(color_PS, bin_num, bins, z_step, confidence)
    print("histograms made")

    plt.figure(-1)
    plt.scatter(redshift, color_PS, s = 0.01, color = 'b', label = "Blue")
    z_edges = np.linspace(min_z, max_z, num = z_step)
    z_middles = z_edges[:-1] + (z_edges[1] - z_edges[0])/2.

    idx = (cuts > 0.45) & (cuts < 1.5)
    z_middles = z_middles[idx]
    cuts = cuts[idx]

    output = np.hstack((z_middles.reshape(z_middles.size,1), cuts.reshape(cuts.size,1)))
    np.savetxt("%s/spline_data.dat" % output_dir, output)


    spl = UnivariateSpline(z_middles, cuts)

    SDSS_data = np.genfromtxt("/work/dominik.zuercher/Output/match_PS/spline_data.dat")
    z_middles_SDSS = SDSS_data[:,0]
    cuts_SDSS = SDSS_data[:,1]

    spl_SDSS = UnivariateSpline(z_middles_SDSS, cuts_SDSS)

    SDSS_with_err_data = np.genfromtxt("/work/dominik.zuercher/Output/match_PS/spline_data_with_errors.dat")
    z_middles_SDSS_with_err = SDSS_with_err_data[:,0]
    cuts_SDSS_with_err = SDSS_with_err_data[:,1]

    spl_SDSS_with_err = UnivariateSpline(z_middles_SDSS_with_err, cuts_SDSS_with_err)


    yy0 = cut(rrange)

    plt.plot(rrange, yy0, 'g-',lw=0.8, label="orig. cut")
    plt.plot(rrange, spl(rrange), 'k-', lw=0.8, label = "GAMA spline")
    plt.plot(rrange, spl_SDSS(rrange), 'r-', lw=0.8, label = "SDSS spline")
    plt.plot(rrange, spl_SDSS_with_err(rrange), 'm-', lw=0.8, label = "SDSS spline with err")

    plt.ylim([0,3])
    plt.xlim([0.03,0.33])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("g-r")
    plt.savefig("%s/color_vs_z.pdf" % output_dir)






