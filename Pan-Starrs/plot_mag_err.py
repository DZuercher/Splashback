import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def half_gaussian(x, c, sigma):
    return c * np.exp(-(x)**2.0/ (2.0*sigma**2.0))

catalog = pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Originals/short_catalog.csv", usecols = [7,8,11,12], names = ['red','rederr','green','greenerr'])

red = np.asarray(catalog['red'].values)
rederr = np.asarray(catalog['rederr'].values)
idx = (red > -900) & (red < 21.5) & (rederr > -900) & (rederr < 5)
red = red[idx]
rederr = rederr[idx]

hist, edges = np.histogram(rederr,bins=1000)
idx = hist < 50
middles = edges[:-1] + 0.5 * (edges[1]-edges[0])
foo = middles[idx]
foo = foo[1]
#popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist)
print("Red error: %s" % foo)



green = np.asarray(catalog['green'].values)
greenerr = np.asarray(catalog['greenerr'].values)
idx = (green > -900) & (green < 21.5) & (greenerr > -900) & (greenerr < 5)
green = green[idx]
greenerr = greenerr[idx]

hist, edges = np.histogram(greenerr,bins=1000)
idx = hist < 50
middles = edges[:-1] + 0.5 * (edges[1]-edges[0])
foo = middles[idx]
foo = foo[1]
#popt, pcov = curve_fit(half_gaussian, edges[:-1] + (edges[1] - edges[0])/2., hist)
print("Green error: %s" % foo)

catalog = pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Originals/short_catalog_onlz_imag.csv", usecols =[11,12], names=['i','ierr'])

iband = np.asarray(catalog['i'].values)
ibanderr = np.asarray(catalog['ierr'].values)
idx = (iband > -900) & (iband < 21.5) & (ibanderr > -900) & (ibanderr < 5)
iband = iband[idx]
ibanderr = ibanderr[idx]
print(np.sum(ibanderr > 0.1))

plt.scatter(red,rederr,c='r',s=0.2,label='r band')
plt.scatter(green,greenerr,c='g',s=0.2,label='g band')
plt.scatter(iband,ibanderr,c='k',s=0.2,label='i band')
plt.legend()
plt.xlabel('mag')
plt.ylabel('mag err')
plt.savefig('/work/dominik.zuercher/Output/Pan-Starrs/err_vs_col.pdf')
