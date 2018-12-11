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

    data = np.loadtxt(catalog)

    #Remove Nans
    idx = (np.logical_not(np.isnan(data[:,1]))) & (np.logical_not(np.isnan(data[:,2])))
    data = data[idx]

    redshift = data[:,0]
    red   = data[:,1]
    green = data[:,2]

    print("Catalog read")

    color = green - red
    fig = plt.figure(1)
    plt.scatter(red, color, s = 0.01)
    fig.savefig("%s/g-r_vs_r.pdf" % output_dir)

