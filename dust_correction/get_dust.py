#Uses mwdust module to get the extinction corrections from MW dust. Run which each color band separately.
#! coding=utf8
#NOTE: In order for mwdust module to work require Python 2.7.5 and Numpy version 1.13 EXACTLY!!!

import mwdust
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import frogress
import csv
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
from glob import glob
from mpi4py import MPI
import os.path



def procedure(rank, size):
 
	steps = 2643.0 - 635.0
	chunky = int(steps/size)
        rest = steps - chunky*size
        mini = chunky*rank
        maxi = chunky*(rank + 1)
        if rank >= (size - 1 - rest):
                maxi += 2 + rank - size + rest
		mini += rank - size + 1 + rest
	mini += 635
	maxi += 635
	mini = int(mini)
	maxi = int(maxi)

        if color == "iband":
            filter_ = "SDSS i"
        elif color == "red":
            filter_ = "SDSS r"
        elif color == "green":
            filter_ = "SDSS g"
        else:
            print("Defined filter not available. Choose red, green or iband")
            return

	imap = mwdust.Combined15(filter = filter_, sf10 = True)


	print("Creating Grid...")
	
	foos = np.arange(mini, maxi, 1)
	full = tools.create_grid(foos)
	print("Kernel %s doing %s up to %s" % (rank, mini, maxi))
	print("Calculating borders of the cells")
	minras, maxras, mindecs, maxdecs = tools.get_centers(full[:,0], full[:,1], borders = True, maskdetect = True)

	print("Using %s objects in Kernel %s" % (full.shape[0], rank))


	print("Converting coordinates to galactic system")
	idy = np.where(np.abs(minras + 99.0) < 0.0001))[0]
	idy2 = np.where(np.abs(minras + 999.0) < 0.0001))[0]
	idy = np.append(idy, idy2)
	idy = np.unique(idy)
	for it in idy:
		minras[it] = 0.0
		mindecs[it] = 0.0
		maxras[it] = 0.0
		maxdecs[it] = 0.0
	c_maxs = SkyCoord(ra = maxras*u.degree, dec = maxdecs*u.degree, frame = 'icrs')
	c_mins = SkyCoord(ra = minras*u.degree, dec = mindecs*u.degree, frame = 'icrs')
	c_maxs = c_maxs.galactic
	c_mins = c_mins.galactic
	ramins = c_mins.l.degree
	ramaxs = c_maxs.l.degree
	decmaxs = c_mins.b.degree
	decmins = c_maxs.b.degree
	print("Calculating extinction...")
	dustyness = np.zeros(0)
	for idx in frogress.bar(range(full.shape[0])):
		if idx in idy:
			dustyness = np.append(dustyness, -99.0)
			continue
		dustynesscolor = np.zeros(0)
		for b in np.linspace(decmins[idx], decmaxs[idx], 5):
                	for l in np.linspace(ramins[idx], ramaxs[idx], 5):
				if color == "red":
                                    foo = rmap(l, b, depth)
				elif color == "iband":
                                    foo = imap(l, b, depth)
                                elif color == "green":
				    foo = gmap(l, b, depth)
                        	dustynesscolor = np.append(dustynesscolor, foo)
        	dustynesscolor = np.average(dustynesscolor)
		dustyness = np.append(dustyness, dustynesscolor)
        np.savetxt(str(outdir)+str(color)+"parts/dust_catalog_"+str(color)+"_"+str(rank)+".csv", np.hstack((full, dustyness.reshape(dustyness.size, 1))))


if __name__ == "__main__":

    color = "iband" #Can be iband, red or green
    depth = 1000 #distance from galactic center in kpc used for dust absorbtion evalution
    resolution = 0.01 #Number of steps per skycellside
    outdir = "./"

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    
    procedure(rank,size)
    
    comm.Barrier()
