import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

h=0.7

cat = pd.read_csv('/work/dominik.zuercher/DataStore/Planck-SZ/Catalogs/Selected_Planck_clusters.csv', sep=',')

MSZ = cat['MSZ'].values*h
z = cat['REDSHIFT'].values

MSZ = MSZ[((z <= 0.3) & (z > 0.03))]
z = z[((z <= 0.3) & (z > 0.03))]


fig = plt.figure(tight_layout = True)
gs = gridspec.GridSpec(5,5)

ax = fig.add_subplot(gs[:4,:4])
ax.plot(z, MSZ, 'k.')
ax.set_ylabel(r'M$_{500\mathrm{c}}$ [h$^{-1} 10^{14}$ M$_{\odot}$]', fontsize=20)
ax.xaxis.set_ticklabels([])
ax.tick_params(labelsize=15)

ax = fig.add_subplot(gs[-1,:-1])
ax.hist(z, bins=50,histtype='step',color='k')
ax.set_xlabel('z', fontsize=20)
ax.tick_params(labelsize=15)

ax = fig.add_subplot(gs[:-1,-1])
ax.hist(MSZ, bins=50, orientation='horizontal',histtype='step',color='k')
ax.yaxis.set_ticklabels([])
ax.tick_params(labelsize=15)

fig.savefig('Planck_summary.pdf')

"""
fig = plt.figure(1)
plt.hist(MSZ, bins=50)
plt.xlabel(r'MSZ [h$^{-1} 10^{14}$ Msol]')
plt.savefig('MSZ_histogram.png')

ig = plt.figure(2)
plt.hist(z, bins=50)
plt.xlabel('z')
plt.savefig('z_histogram.png')

fig = plt.figure(3)
plt.plot(z, MSZ, '.')
plt.xlabel('z')
plt.ylabel('MSZ [h$^{-1} 10^{14}$ Msol]')
plt.savefig('MSZ_vs_z.png')
"""


