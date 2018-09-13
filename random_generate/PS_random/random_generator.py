#! coding=utf8
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as uwell
import astropy.wcs as wcs
import pandas as pd
import matplotlib.pyplot as plt
import math
import healpy as hp
import os
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import argparse
from mpi4py import MPI
from glob import glob


def init_procedure(rank):	
    total_randoms = int(total_num*4.0/3.0/size) #Randoms generated per thread

    #Generate the random objects
    ran_dec  = -1.0 + 2.*np.random.random(size = total_randoms)
    ran_ra = 2.0*np.pi*np.random.random(size = total_randoms)
    ran_dec = np.arccos(ran_dec)
    print("generated")

    #Throw away those in bad skycells using bads.fits mask
    bads = hp.read_map("%s/bads_%s.fits" % (mask_path, mag_cut))

    print("converting")
    pix = hp.pixelfunc.ang2pix(NSIDE, ran_dec, ran_ra)
    pixmask = (bads[pix] == 1)
    ran_ra = ran_ra[pixmask]
    ran_dec = ran_dec[pixmask]

    print("selected")
    ran_dec = np.degrees(np.pi/2.0 - ran_dec)
    ran_ra = np.degrees(ran_ra)
    output = np.hstack((ran_ra.reshape((ran_ra.size, 1)), ran_dec.reshape((ran_dec.size, 1))))
    output = np.hstack((output, np.zeros((output.shape[0], 2))))
    projectioncell = tools.decide_pro(output[:,0], output[:,1])

    print("Decided on Projectioncell")
    output[:,2] = projectioncell.reshape((output.shape[0],))
    outdict = {key:np.zeros((1, 4)) for key in range(635, 2644)}

    print("Writing Directory...")	
    for j in range(output.shape[0]):
        if j == int(output.shape[0]/4):
            print("%s at 25%" % rank)
        if j == int(output.shape[0]/2):
            print("%s at 50%" % rank)
        cell = output[j,2]
        outdict[int(cell)] = np.vstack((outdict[int(cell)], output[j,:]))

    print("Saving")
    for i in range(635, 2644): 
        if outdict[int(i)][1:,:].size != 0:
            out = pd.DataFrame(outdict[int(i)][1:,:])
            try:
                np.savetxt("%s/random_cell_%s_rank_%s.dat" % (outdir, int(i), rank), outdict[int(i)][1:,:])
            except:
                print("cannot save")
	
def main_procedure(rank):

    steps = 2643.0 - 635.0
    
    chunky = int(steps/size)
    rest = steps - chunky*size
    mini = chunky*rank
    maxi = chunky*(rank + 1)
    if rank >= (size - 1 -rest):
        maxi += 2 + rank - size + rest
        mini += rank - size + 1 + rest
    mini += 635
    maxi += 635
    mini = int(mini)
    maxi = int(maxi)

    for cellnumber in np.arange(mini, maxi, 1):
        try:
            df = pd.read_csv("%s/random_cell_%s.dat" % (indir, cellnumber), sep=' ', header=None, names=(['ra','dec','pro','sky']))
            output = df.values
        except:
            continue
        if output.size == 0:
            continue
        print("Deciding on Skycell...")

        skycells = tools.decide_sky(output[:,0], output[:,1], cellnumber)
        if skycells.size == 0:
            continue
        output[:,3] = skycells.reshape((skycells.shape[0],))
        np.savetxt("%s/random_cell_%s.dat" % (middir, cellnumber), output)

        try:
            df = pd.read_csv("%s/random_cell_%s.dat" % (middir, cellnumber), sep = ' ', header = None, names=(['ra','dec','pro','sky']))
            output = df.values
        except:
            continue
        if output.size == 0:
            continue
        
        selectedobjects = tools.mask_select(output, cellnumber)
        
        np.savetxt("%s/random_cell_%s.csv" % (outdir, cellnumber), selectedobjects, delimiter = ',')
    

if __name__ == "__main__":

    
    NSIDE = 2048
    total_num = 1155286477.0  #Total number of objects in final PS catalog (approx ok)
    hdulist = fits.open('/home/dominik.zuercher/Documents/RSP_Pro/BCG_correction/Pan-STARRS_Pics/ps1grid.fits') #Pan-STARRS grid file
    mask_path =  "/home/dominik.zuercher/Documents/RSP_Pro/depth_analysis"

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--init", help = "If true new randoms are generated")
    parser.add_argument("--mag_cut", help = "Mag cut (21, 21.5 or 22)")
    parser.add_argument("--indir", help="Input Directory (only for init == False")
    parser.add_argument("--outdir", help="Output Directory")
    parser.add_argument("--middir", help="Intermediate Directory (only for init == False)")
    args = parser.parse_args()
    init = args.init
    mag_cut = args.mag_cut
    outdir = args.outdir
    indir = args.indir
    middir = args.middir

    if init == True:	
            init_procedure(rank)
    else:
            main_procedure(rank)
    print("%s finished" % rank)
    comm.Barrier()

