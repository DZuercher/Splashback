from scipy.interpolate import UnivariateSpline
from  scipy.optimize import curve_fit
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import frogress
import scipy.integrate as integrate

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

def get_cuts(color, bin_num, bins, z_steps, confidence, plotting):
    cuts = np.zeros(z_steps-1)
    blue_con = np.zeros(z_steps-1)
    red_con = np.zeros(z_steps-1)
    for it in frogress.bar(range(1,z_steps)):
        colors = color[bin_num == it]
        hist, edges = np.histogram(colors, bins = bins)
        try:
            popt, pcov = curve_fit(double_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist)
        except:
            cuts[it-1] = -1000
            blue_con[it-1] = -1000
            red_con[it-1] = -1000
            continue
        red_ind = np.argmax([popt[1], popt[4]])
        if red_ind == 0:
            red_c = popt[0]
            red_mu = popt[1]
            red_sigma = popt[2]
            blue_c = popt[3]
            blue_mu = popt[4]
            blue_sigma = popt[5]
        else:
            blue_c = popt[0]
            blue_mu = popt[1]
            blue_sigma = popt[2]
            red_c = popt[3]
            red_mu = popt[4]
            red_sigma = popt[5]

	#resacaling height to simulate cluster region
	red_integral = red_c*np.sqrt(2*np.pi)*np.abs(red_sigma)
	blue_integral = blue_c*np.sqrt(2*np.pi)*np.abs(blue_sigma)

	red_c_old = red_c
	red_c = red_c * 0.65*blue_integral/(red_integral*(1. - 0.65))

	red_integral = red_c*np.sqrt(2*np.pi)*np.abs(red_sigma)
	blue_integral = blue_c*np.sqrt(2*np.pi)*np.abs(blue_sigma)
        frac = red_integral/(red_integral + blue_integral)

	#Get cut
	cut = red_mu - confidence * np.abs(red_sigma)
        cuts[it - 1] = cut

	#broaden standard deviation by magnitude errors
	red_sigma_broad = np.sqrt( red_sigma**2.0  + red_err**2.0 + green_err**2.0 )
	blue_sigma_broad = np.sqrt( blue_sigma**2.0  + red_err**2.0 + green_err**2.0 )

	blues = integrate.quad( lambda x : blue_c * np.exp( - (x - blue_mu)**2.0 / (2.0 * blue_sigma_broad**2.0) ), -np.inf, cut)[0]
	reds = integrate.quad( lambda x : red_c * np.exp( - (x - red_mu)**2.0 / (2.0 * red_sigma_broad**2.0) ), -np.inf, cut)[0]
	blue_con[it - 1] = blues
	red_con[it - 1] = reds

        if (it % 10 == 0) & (plotting == True):
            plt.figure(it / 10)
            plt.hist(colors, bins = bins)
            gr_range = edges[:-1] + (edges[1] - edges[0])/2.
            plt.plot(gr_range, blue_c* np.exp( - (gr_range - blue_mu)**2.0 / (2.0 * blue_sigma**2.0) ), 'b-', lw = 0.3)
            plt.plot(gr_range, red_c_old* np.exp( - (gr_range - red_mu)**2.0 / (2.0 * red_sigma**2.0) ), 'r-', lw = 0.3)
            plt.plot(gr_range, blue_c* np.exp( - (gr_range - blue_mu)**2.0 / (2.0 * blue_sigma_broad**2.0) ), 'b--', lw = 0.3)
            plt.plot(gr_range, red_c* np.exp( - (gr_range - red_mu)**2.0 / (2.0 * red_sigma_broad**2.0) ), 'r--', lw = 0.3)
            plt.axvline(cut, color = 'k', lw = 0.3)
            plt.savefig("%s/test_hist_%s.pdf" % (output_dir, it/ 10 ))

    return cuts, red_con, blue_con


if __name__ == '__main__':

    #max values (3sigma in 20.0 to 21.5 mag)
    red_err = 0.12
    green_err = 0.09
    iband_err = 0.13

    #Error on magnitudes at 20.0 to 21.0 mag
    red_err = 0.132
    green_err = 0.261
    iband_err = 0.099

    min_z = 0.03
    max_z = 0.33
    z_step = 100
    bins = 500

    confidence = 6 # Confidence level to exclude Reds from Blues (in sigmas)

    catalog = "/work/dominik.zuercher/Output/match_PS/matched_spec_new.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS"


    data = np.loadtxt(catalog)
    redshift = data[:,0]
    red_PS = data[:,3]
    green_PS = data[:,4]
    iband_PS = data[:,5]
    id_ = data[:,6]
    idx = (red_PS < 900) & (red_PS > -900) & (green_PS < 900) & (green_PS > -900) & (iband_PS < 900) & (iband_PS > -900)

    redshift = redshift[idx]
    red_PS = red_PS[idx]
    green_PS = green_PS[idx]
    iband_PS = iband_PS[idx]
    id_ = id_[idx]
    print("Catalog read")


    colors = [green_PS - iband_PS, red_PS - iband_PS, green_PS - red_PS]
    labels = ["g-i", "r-i", "g-r"]
    for i in range(3):

	bin_num = bin_redshifts(redshift, min_z, max_z, z_step)
	print("Redshift binning done")
	cuts, red_con, blue_con = get_cuts(colors[i], bin_num, bins, z_step, confidence, plotting = (i == 2))
	print("histograms made")

	plt.figure(-1)
	plt.scatter(redshift[id_ == True], colors[i][id_ == True], s = 0.01, color = 'b', label = "Blue")
	plt.scatter(redshift[id_ == False], colors[i][id_ == False], s = 0.01, color = 'r', label = "Red")
	z_edges = np.linspace(min_z, max_z, num = z_step)
	z_middles = z_edges[:-1] + (z_edges[1] - z_edges[0])/2.

	idx = (cuts > 0.0) & (cuts < 1.5)
	z_middles = z_middles[idx]
	cuts = cuts[idx]
	red_con = red_con[idx]
	blue_con = blue_con[idx]
	if(i == 2):	
	    spl = UnivariateSpline(z_middles, cuts)
	    rrange = np.linspace(0, 0.35, 1000)

	    plt.plot(rrange, spl(rrange), 'k-', lw=0.3, label = "spline")

	    plt.ylim([0,3])
	    plt.xlim([0.03,0.33])
	    plt.legend()
	    plt.xlabel("z")
	    plt.ylabel(labels[i])
	    plt.savefig("%s/%s_vs_z.pdf" % (output_dir, labels[i]))

	idx = (red_con > 0.0) & (blue_con > 0.0)
	red_con = red_con[idx]
	blue_con = blue_con[idx]
	contamination = np.sum(red_con)/np.sum(red_con + blue_con)
	print("Contamination calculated with color %s amounts to %s" % (labels[i], 100*contamination))
