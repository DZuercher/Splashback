import healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
cat_path = "/work/dominik.zuercher/DataStore/Planck-SZ/Catalogs/21_bins/Planck_21_bins.csv"

planck = pd.read_csv(cat_path)
planck = planck[(planck['REDSHIFT'] <= 0.33) & (planck['REDSHIFT'] >= 0.03)]
ra = np.asarray(planck["RA"].values)
dec = np.asarray(planck["DEC"].values)

f, axarr = plt.subplots(ncols=3, nrows=1, figsize=(18,6))

sky1 = hp.read_map("/work/dominik.zuercher/Output/depth_analysis/bads_21.fits")
plt.axes(axarr[0])
hp.mollview(sky1, hold=True, cbar=False, title="PS 21")
hp.projscatter(ra, dec, lonlat=True, s = 1, c = 'r')
#hp.projscatter(rand_ra, rand_dec, lonlat=True, s = 3, c = 'b')
hp.graticule()

sky1 = hp.read_map("/work/dominik.zuercher/Output/depth_analysis/bads_21.5.fits")
plt.axes(axarr[1])
hp.mollview(sky1, hold=True, cbar=False, title="PS 21.5")
hp.projscatter(ra, dec, lonlat=True, s = 1, c = 'r')
#hp.projscatter(rand_ra, rand_dec, lonlat=True, s = 3, c = 'b')
hp.graticule()

sky1 = hp.read_map("/work/dominik.zuercher/Output/depth_analysis/bads_22.fits")
plt.axes(axarr[2])
hp.mollview(sky1, hold=True, cbar=False, title="PS 22")
hp.projscatter(ra, dec, lonlat=True, s = 1, c = 'r')
#hp.projscatter(rand_ra, rand_dec, lonlat=True, s = 3, c = 'b')
hp.graticule()

plt.subplots_adjust(wspace=0.5)
f.savefig("/work/dominik.zuercher/Output/skymaps.pdf")
