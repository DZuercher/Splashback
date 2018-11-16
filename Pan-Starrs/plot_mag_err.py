from scipy import stats
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def half_gaussian(x, c, mu, sigma):
    return c * np.exp(-(x-mu)**2.0/ (2.0*sigma**2.0))

catalog = pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Originals/short_catalog.csv", usecols = [7,8,11,12], names = ['red','rederr','green','greenerr'])

plt.figure(1)


red = np.asarray(catalog['red'].values)
rederr = np.asarray(catalog['rederr'].values)
idx = (red > -900) & (red < 21.9) & (rederr > -900) & (red > 15.4) #+ 0.4
red = red[idx]
rederr = rederr[idx]
red_med = np.median(rederr)


hist, edges = np.histogram(rederr,bins=1000)
n, plotbins, patches = plt.hist(rederr, bins=50, color='r', alpha=0.5, normed=True)
popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist, bounds = ([0.0, 0.0, 0.0], [ np.inf, 0.01, 0.1]) )
foo = popt[2]
plt.axvline(3*foo,color='r')

print("Red error: %s" % red_med)

green = np.asarray(catalog['green'].values)
greenerr = np.asarray(catalog['greenerr'].values)
idx = (green > -900) & (green < 22.8) & (greenerr > -900) & (green > 16.3) #+1.3
green = green[idx]
greenerr = greenerr[idx]
green_med = np.median(greenerr)

hist, edges = np.histogram(greenerr,bins=1000)
plt.hist(greenerr, bins=plotbins, color='g', alpha=0.5, normed=True)
popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist, bounds = ([0.0, 0.0, 0.0], [ np.inf, 0.01, 0.1]) )
foo = popt[2]
print("Green error: %s" % green_med)
plt.axvline(3*foo,color='g')

catalog = pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Originals/short_catalog_onlz_imag.csv", usecols =[11,12], names=['i','ierr'])

iband = np.asarray(catalog['i'].values)
ibanderr = np.asarray(catalog['ierr'].values)
idx = (iband > -900) & (iband < 21.5) & (ibanderr > -900) & (iband > 15.0) 
iband = iband[idx]
ibanderr = ibanderr[idx]
i_med = np.median(ibanderr)

hist, edges = np.histogram(ibanderr,bins=1000)
plt.hist(ibanderr, bins=plotbins,color='k', alpha=0.5, normed=True)
popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist, bounds = ([0.0, 0.0, 0.0], [ np.inf, 0.01, 0.1]) )
foo = popt[2]
print("iband error: %s" % i_med)
plt.axvline(3*foo,color='k')

text = "$i_{PS}$ error: 0.063 mag \n $r_{PS}$ error: 0.064 mag \n $g_{PS}$ error: 0.081 mag"
plt.text(0.18, 38, s=text)
plt.xlabel("mag. err.")
plt.xlim([0.0, 0.25])
plt.savefig("/work/dominik.zuercher/Output/Pan-Starrs/errors.pdf")


rvals, rbins, foo = stats.binned_statistic(red, rederr, 'median', bins=100)
rbins = rbins[:-1] + 0.5*(rbins[1:] - rbins[:-1])
gvals, gbins, foo = stats.binned_statistic(green, greenerr, 'median', bins=100)
gbins = gbins[:-1] + 0.5*(gbins[1:] - gbins[:-1])
ivals, ibins, foo = stats.binned_statistic(iband, ibanderr, 'median', bins=100)
ibins = ibins[:-1] + 0.5*(ibins[1:] - ibins[:-1])


np.savetxt("/work/dominik.zuercher/r_band_error.txt",np.hstack((rbins.reshape(-1, 1), rvals.reshape(-1,1))))
np.savetxt("/work/dominik.zuercher/g_band_error.txt",np.hstack((gbins.reshape(-1, 1), gvals.reshape(-1,1))))
np.savetxt("/work/dominik.zuercher/i_band_error.txt",np.hstack((ibins.reshape(-1, 1), ivals.reshape(-1,1))))
red_med = rvals[-1]
green_med = gvals[-1]
i_med = ivals[-1]


f, axarr = plt.subplots(ncols=3,nrows=1,figsize=(11,3), tight_layout=True)

axarr[0].scatter(red,rederr,c='r',s=0.01,label='r band',marker='.')
axarr[1].scatter(green,greenerr,c='g',s=0.01,label='g band',marker='.')
axarr[2].scatter(iband,ibanderr,c='k',s=0.01,label='i band',marker='.')

axarr[0].plot(rbins, rvals, 'k--')
axarr[1].plot(gbins, gvals, 'k--')
axarr[2].plot(ibins, ivals, 'w--')

"""
axarr[0].axhline(red_med, c='k', ls='--', lw=0.5)
axarr[1].axhline(green_med, c='k', ls='--', lw=0.5)
axarr[2].axhline(i_med, c='k', ls='--', lw=0.5)
"""

#axarr[0].text(15.5,0.28,"Mag Error r-band : %d mmag" % (red_med*1000),fontsize=7)
#axarr[1].text(16.5,0.28,"Mag Error g-band: %d mmag" % (green_med*1000),fontsize=7)
#axarr[2].text(15.2,0.28,"Mag Error i-band: %d mmag" % (i_med*1000),fontsize=7)

axarr[2].set_xlim([15.0, 21.5])
axarr[1].set_xlim([16.3, 22.8])
axarr[0].set_xlim([15.4, 21.9])

axarr[0].set_ylim([0.0, 0.3])
axarr[1].set_ylim([0.0, 0.3])
axarr[2].set_ylim([0.0, 0.3])

axarr[0].set_xlabel(r'r$_{\mathrm{P1}}$')
axarr[0].set_ylabel(r'r$_{\mathrm{P1}}$ error')
axarr[1].set_xlabel(r'g$_{\mathrm{P1}}$')
axarr[1].set_ylabel(r'g$_{\mathrm{P1}}$ error')
axarr[2].set_xlabel(r'i$_{\mathrm{P1}}$')
axarr[2].set_ylabel(r'i$_{\mathrm{P1}}$ error')
f.savefig('mag_errors.png')
