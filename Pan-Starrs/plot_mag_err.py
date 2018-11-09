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
idx = (red > -900) & (red < 20.4) & (rederr > -900) & (red > 15.4) #+ 0.4
red = red[idx]
rederr = rederr[idx]

hist, edges = np.histogram(rederr,bins=1000)
n, plotbins, patches = plt.hist(rederr, bins=50, color='r', alpha=0.5, normed=True)
popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist, bounds = ([0.0, 0.0, 0.0], [ np.inf, 0.01, 0.1]) )
foo = popt[2]
plt.axvline(3*foo,color='r')

print("Red error: %s" % (3*foo))



green = np.asarray(catalog['green'].values)
greenerr = np.asarray(catalog['greenerr'].values)
idx = (green > -900) & (green < 21.3) & (greenerr > -900) & (green > 16.3) #+1.3
green = green[idx]
greenerr = greenerr[idx]

hist, edges = np.histogram(greenerr,bins=1000)
plt.hist(greenerr, bins=plotbins, color='g', alpha=0.5, normed=True)
popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist, bounds = ([0.0, 0.0, 0.0], [ np.inf, 0.01, 0.1]) )
foo = popt[2]
print("Green error: %s" % (3*foo))
plt.axvline(3*foo,color='g')

catalog = pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Originals/short_catalog_onlz_imag.csv", usecols =[11,12], names=['i','ierr'])

iband = np.asarray(catalog['i'].values)
ibanderr = np.asarray(catalog['ierr'].values)
idx = (iband > -900) & (iband < 20.0) & (ibanderr > -900) & (iband > 15.0) 
iband = iband[idx]
ibanderr = ibanderr[idx]

hist, edges = np.histogram(ibanderr,bins=1000)
plt.hist(ibanderr, bins=plotbins,color='k', alpha=0.5, normed=True)
popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist, bounds = ([0.0, 0.0, 0.0], [ np.inf, 0.01, 0.1]) )
foo = popt[2]
print("iband error: %s" % (3*foo))
plt.axvline(3*foo,color='k')

text = "$i_{PS}$ error: 0.063 mag \n $r_{PS}$ error: 0.064 mag \n $g_{PS}$ error: 0.081 mag"
plt.text(0.18, 38, s=text)
plt.xlabel("mag. err.")
plt.xlim([0.0, 0.25])
plt.savefig("/work/dominik.zuercher/Output/Pan-Starrs/errors.pdf")




"""
plt.scatter(red,rederr,c='r',s=0.2,label='r band')
plt.scatter(green,greenerr,c='g',s=0.2,label='g band')
plt.scatter(iband,ibanderr,c='k',s=0.2,label='i band')
plt.legend()
plt.xlabel('mag')
plt.ylabel('mag err')
plt.savefig('/work/dominik.zuercher/Output/Pan-Starrs/err_vs_col.pdf')
"""
