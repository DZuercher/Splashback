#Combines extinction corrections from different filters into one file
import pandas as pd
import numpy as np
import glob

if __name__ == "__main__":

    size = 980 #Number of processors used in get_dust.py
    dust_dir = "./" #Directory containing the dust catalogs produced by get_dust.py
    outfile = "test.dat"

    output = np.zeros((1, 5))

    for ii in range(size):
        redfile = np.loadtxt(outdir + "redparts/dust_catalog_red_%s.csv" % ii)
        greenfile = np.loadtxt(outdir + "greenparts/dust_catalog_green_%s.csv" % ii)
        ifile = np.loadtxt(outdir + "ibandparts/dust_catalog_iband_%s.csv" % ii)
        green = greenfile[:,2]
        iband = ifile[:,2]
        green = green.reshape(green.size,1)
        iband = iband.reshape(iband.size,1)
        foo = np.hstack((redfile, green, iband))
        output = np.vstack((output, foo))

    np.savetxt(outfile, output[1:,:])





