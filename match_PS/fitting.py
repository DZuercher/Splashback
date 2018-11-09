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

    catalog = "/work/dominik.zuercher/Output/match_PS/matched_spec_new.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS"


    data = np.loadtxt(catalog)
    redshift = data[:,0]
    red_PS = data[:,3]
    green_PS = data[:,4]
    iband_PS = data[:,5]
    id_ = data[:,6]

    idx = (red_PS < 900) & (red_PS > -900) & (green_PS < 900) & (green_PS > -900) 
    redshift = redshift[idx]
    red_PS = red_PS[idx]
    green_PS = green_PS[idx]
    id_ = id_[idx]
    color = red_PS - green_PS
    print("Catalog read")


    bin_num = bin_redshifts(redshift, min_z, max_z, z_step)
    print("Redshift binning done")
    cuts = get_cuts(color, bin_num, bins, z_step, confidence, False)
    #if calc_contamination == True:
#	cuts_err = get_cuts(color, bin_num, bins, z_step, confidence, True)
    print("histograms made")




    plt.figure(-1)
    plt.scatter(redshift[id_ == True], color[id_ == True], s = 0.01, color = 'b', label = "Blue")
    plt.scatter(redshift[id_ == False], color[id_ == False], s = 0.01, color = 'r', label = "Red")
    z_edges = np.linspace(min_z, max_z, num = z_step)
    z_middles = z_edges[:-1] + (z_edges[1] - z_edges[0])/2.

    idx = (cuts > 0.0) & (cuts < 1.5)
    z_middles_1 = z_middles[idx]
    cuts = cuts[idx]
    #cuts_err = cuts_err[idx]

    spl = UnivariateSpline(z_middles_1, cuts)

    #spl_err = UnivariateSpline(z_middles_1, cuts_err)


    rrange = np.linspace(0, 0.35, 1000)
    yy0 = cut(rrange)

    
    output = np.hstack((z_middles.reshape(z_middles.size,1), cuts.reshape(cuts.size,1)))
 #   if include_PS_errors == False:
    np.savetxt("%s/spline_data.dat" % output_dir, output)
 #   else:
  #      np.savetxt("%s/spline_data_with_errors.dat" % output_dir, output)

    #np.savetxt("%s/redshifts.dat" % output_dir, z_middles.reshape(z_middles.size,1))
    #np.savetxt("%s/colors.dat" % output_dir, cuts.reshape(z_middles.size,1))

    #plt.plot(z_middles, cuts, 'k-', lw = 0.3, label = "cut")
    #plt.plot(rrange, yy0, 'g-',lw=0.3, label="orig. cut")
    plt.plot(rrange, spl(rrange), 'k-', lw=0.3, label = "spline")
    #plt.plot(rrange, spl_err(rrange), 'k-', lw=0.3, label = "spline (inc. PS errors)")

    plt.ylim([0,3])
    plt.xlim([0.03,0.33])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("g-r")
    plt.savefig("%s/color_vs_z.pdf" % output_dir)


    """
    #Calculate contaminations for mixed characterization (best estimate?)
    #Falses
    reds_below_cut = (color < spl(redshift) ) & (color > spl_err(redshift))
    #blues_above_cut = (id_ == True) & (color >= spl(redshift) )
    #Rights
    #reds_above_cut = (id_ == False) & (color >= spl(redshift) )
    blues_below_cut = (color < spl_err(redshift) )

    print("Reds below cut: %s" % np.sum(reds_below_cut))
    print("Blues below cut: %s" % np.sum(blues_below_cut))
    print("Contamination: %f" % (np.sum(reds_below_cut)/(np.sum(reds_below_cut) + np.sum(blues_below_cut)))*100.0 )
    #print("Reds above cut: %s" % np.sum(reds_above_cut))
    #print("Blues above cut: %s" % np.sum(blues_above_cut))
    """





