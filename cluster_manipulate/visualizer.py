# Plots a map of the (random) cluster positions

import healpy as hp
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy import units as uwell
from astropy.coordinates import SkyCoord


cat_path = "catalog.csv"
rand_cat_path = "random_catalog.csv"
sel_func_path =" HFI_PCCS_SZ-selfunc-union-survey_R2.08.fits"
title = "Planck Clusters"
out_path = "Planck_scatter.pdf"


if __name__=="__main__":
    catalog = pd.read_csv(cat_path)
    random_catalog = pd.read_csv(rand_cat_path)

    ra = np.asarray(catalog["RA"].values)
    dec = np.asarray(catalog["DEC"].values)
    rand_ra = np.asarray(random_catalog["RA"].values)
    rand_dec = np.asarray(random_catalog["DEC"].values)

    seln_func = hp.read_map(sel_func_path)
    hp.visufunc.mollview(seln_func, title = title, coord=["G","C"])
    hp.projscatter(ra, dec, lonlat=True, s = 3, c = 'r')
    #hp.projscatter(rand_ra, rand_dec, lonlat=True, s = 3, c = 'b')
    hp.graticule()
    plt.savefig(out_path)
