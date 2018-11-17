import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np



if __name__ == "__main__":

    min_z = 0.03
    max_z = 0.33
    z_step = 100
    bins = 500

    confidence = 3 # Confidence level to exclude Reds from Blues (in sigmas)

    catalog = "/work/dominik.zuercher/Output/match_PS_GAMA/matched_spec_new.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS_GAMA"

    data     = np.loadtxt(catalog)
    redshift = data[:,0]
    red_PS   = data[:,1]
    green_PS = data[:,2]
    iband_PS = data[:,3]
    id_      = data[:,4]

    idx = (red_PS < 900) & (red_PS > -900) & (green_PS < 900) & (green_PS > -900) & (iband_PS < 900) & (iband_PS > -900)
    redshift = redshift[idx]
    red_PS = red_PS[idx]
    green_PS = green_PS[idx]
    iband_PS = iband_PS[idx]
    id_ = id_[idx]
    print("Catalog read")

    idx = np.random.rand(redshift.size)
    boolar = idx<0.05
    redshift = redshift[boolar]
    red_PS = red_PS[boolar]
    green_PS = green_PS[boolar]
    iband_PS = iband_PS[boolar]
    id_ = id_[boolar]

    color = green_PS - red_PS
    fig = plt.figure(1)
    plt.scatter(red_PS, color, s = 0.01)
    fig.savefig("%s/g-r_vs_r.pdf" % output_dir)

    color = green_PS - iband_PS
    fig = plt.figure(2)
    plt.scatter(iband_PS, color, s = 0.01)
    fig.savefig("%s/g-i_vs_i.pdf" % output_dir)

    color = red_PS - iband_PS
    fig = plt.figure(3)
    plt.scatter(iband_PS, color, s = 0.01)
    fig.savefig("%s/r-i_vs_i.pdf" % output_dir)
