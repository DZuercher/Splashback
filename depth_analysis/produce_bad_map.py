#NOTE: The adding cuts part requires Numpy==1.13 since there is a bug in Numpy 1..14 causing Skycoord to fail if compiled with Python 2.7.5
import healpy as hp
import pandas as pd
from astropy.table import Table
import numpy as np 
import frogress
from astropy import units as uwell
from astropy.coordinates import SkyCoord
from itertools import product
from astropy.io import fits
import astropy.wcs as wcs
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse



def part1():
    print("Reading dustmap...")
    dust = pd.read_csv(dust_path, sep = ' ', usecols = [0,1,-1], names = ['procell','skycell','SDSS_i_avg'], header = None)
    dustprocells = np.asarray(dust["procell"].values).astype(int)
    dustskycells = np.asarray(dust["skycell"].values).astype(int)
    dust_avg = np.asarray(dust["SDSS_i_avg"].values)

    print("Reading Mag lim catalog")
    data = Table.read(cut_path, format = 'fits')
    data = data.to_pandas()
    skycells = np.asarray(data["skyCellID"].values).astype(int)
    procells = np.asarray(data["projectionid"].values).astype(int)
    mag = np.asarray(data["max_iKronMag"].values)

    for ii in frogress.bar(range(skycells.shape[0])):
        skycell = skycells[ii]
        procell = procells[ii]
        idx = np.where((dustprocell s== procell) & (dustskycells == skycell))[0]
        if idx.size != 1:
            print(idx, procell, skycell)
            print("Failed")
            break
        if dust_avg[idx] > -20.0:
            mag[ii] -= dust_avg[idx]	
        else:
            print("No Mask")
            mag[ii] = 0.0
    idx = (mag < mag_cut)
    outpro = procells[idx]
    outpro = outpro.reshape((outpro.size, 1))
    outsky = skycells[idx]
    outsky = outsky.reshape((outsky.size, 1))
    inpu = np.hstack((outpro, outsky))

    print("Calculating ra, dec borders...")
    outpro = inpu[:,0]
    outsky = inpu[:,1]
    minras, maxras, mindecs, maxdecs = tools.get_centers(outpro, outsky, borders = True)
    print("Check fully missing masks...")
    idx = np.where(np.abs(minras - (-99)) < 0.0001)[0]
    idbroken = np.where(np.abs(minras - (-999) ) <0.0001)[0]
    print("There are %s masks that cannot be processed" % idbroken.size)
    brokenpro = np.take(outpro, idx)
    startcell, M, zone, xsub, ysub, xsize, ysize, mindec, maxdec, prodec = tools.get_cell_info(brokenpro)
    step = 360.0/M
    cen = step*(brokenpro - startcell)
    highra = cen + step/2.0
    lowra = cen - step/2.0
    bx = lowra < 0.0
    lowra[bx] = lowra[bx] + 360.0
    minras[idx] = lowra
    mindecs[idx] = mindec
    maxras[idx] = highra
    maxdecs[idx] = maxdec

    output = pd.DataFrame({'pro':outpro,'sky':outsky,'minras':minras,'maxras':maxras,'mindecs':mindecs,'maxdecs':maxdecs})

    output.to_csv(outfile_1, index = False)


def part2(rank):

    print("Calculating Healpy map...")
    output = pd.read_csv(outfile_1)
    outpro = np.asarray(output["pro"].values)
    outsky = np.asarray(output["sky"].values)
    minras = np.asarray(output["minras"].values)
    maxras = np.asarray(output["maxras"].values)
    mindecs = np.asarray(output["mindecs"].values)
    maxdecs = np.asarray(output["maxdecs"].values)

    steps = minras.size        
    chunky = int(steps/size)
    rest = steps - chunky*size
    mini = chunky*rank
    maxi = chunky*(rank + 1)
    if rank >= (size - 1) - rest:
        maxi += 2 + rank - size +rest)
        mini += rank - size + 1 + rest
    if rank == size - 1:
        maxi = steps
    mini = int(mini)
    maxi = int(maxi)

    bad_pixels = np.zeros(0)

    for ii in range(mini, maxi):
        #completely removing the stupid 2643 cell
        if int(outpro[ii]) == 2643:
            minras[ii] = 0.0
            maxras[ii] = 360.0
            mindecs[ii] = 87.969
            maxdecs[ii] = 90.0
        if minras[ii] > maxras[ii]:
            print("ra reverted (at 0-360 border case)")
            print(minras[ii], maxras[ii])
            minras1 = minras[ii]
            maxras1 = 360.0
            minras2 = 0.0
            maxras2 = maxras[ii]
            rarange = np.arange(minras1, maxras1, resolution) #adjusting for healpy resolution of 0.02 deg
            decrange = np.arange(mindecs[ii], maxdecs[ii], resolution)
            c = np.asarray(list(product(rarange, decrange)), dtype = float)
            c[:,1] = np.subtract(90.0, c[:,1]) #conversion to coangle
            for n in range(c.shape[0]):
                if (c[n,1] < 0): #Accounts for convertion problem at poles
                    c[n,1] = 0
            c[:,1] = np.radians(c[:,1])
            c[:,0] = np.radians(c[:,0])
            pix = hp.pixelfunc.ang2pix(NSIDE, c[:,1], c[:,0])
            bad_pixels = np.append(bad_pixels, pix)
            
            rarange = np.arange(minras2, maxras2, resolution) #adjusting for healpy resolution of 0.02 deg
            decrange = np.arange(mindecs[ii], maxdecs[ii], resolution)
            c = np.asarray(list(product(rarange, decrange)), dtype = float)
            c[:,1] = np.subtract(90.0, c[:,1]) #conversion to coangle
            for n in range(c.shape[0]):
                if (c[n,1] < 0): #Accounts for convertion problem at poles
                    c[n,1] = 0
            c[:,1] = np.radians(c[:,1])
            c[:,0] = np.radians(c[:,0])
            pix = hp.pixelfunc.ang2pix(NSIDE, c[:,1], c[:,0])
            bad_pixels = np.append(bad_pixels, pix)
        else:
            rarange = np.arange(minras[ii], maxras[ii], resolution) #adjusting for healpy resolution of 0.02 deg
            decrange = np.arange(mindecs[ii], maxdecs[ii], resolution)
            c = np.asarray(list(product(rarange, decrange)), dtype = float)
            try:
                c[:,1] = np.subtract(90.0, c[:,1]) #conversion to coangle
            except:
                print("Why?")
                print(minras[ii], maxras[ii], mindecs[ii], maxdecs[ii], outpro[ii], outsky[ii])
            for n in range(c.shape[0]):
                if (c[n,1] < 0): #Accounts for convertion problem at poles
                    c[n,1] = 0
            try:
                c[:,1] = np.radians(c[:,1])
                c[:,0] = np.radians(c[:,0])
            except:
                print("Why?")
                print(minras[ii], maxras[ii], mindecs[ii], maxdecs[ii], outpro[ii], outsky[ii])
            pix = hp.pixelfunc.ang2pix(NSIDE, c[:,1], c[:,0])
            bad_pixels = np.append(bad_pixels, pix)

    bad_pixels = np.unique(bad_pixels)
    np.savetxt("%s/bad_pixels_%s.dat" % (outdir_2, rank) , bad_pixels)


def part3():
    bad_pixels = np.genfromtxt(infile_3)
    bad_pixels = np.unique(bad_pixels)
    print("Adding cuts...")
    data = np.zeros(0)

    print("Adding +/- 20 deg band cut for MW")
    #Apply a bandcut of +/-20 in galactic coords to get rid of MW completely
    lspace = np.arange(0.0, 360.0, resolution)
    bspace = np.arange(-20.0, 20.0, resolution) 
    listspace = np.asarray(list(product(lspace, bspace)), dtype = float)
    c = SkyCoord(l = listspace[:,0], b = listspace[:,1], frame = 'galactic', unit = 'degree')
    raspace = c.icrs.ra.deg
    decspace = c.icrs.dec.deg
    raspace = np.radians(raspace)
    decspace = np.radians(np.subtract(90.0, decspace))
    pixcut = hp.pixelfunc.ang2pix(NSIDE, decspace, raspace)

    print("Cut below -31 deg latitiude...")
    raspace = np.arange(0.0, 360.0, resolution)
    raspace = np.radians(raspace)
    decspace = np.arange(-90.0, -31.0, resolution)
    decspace = np.subtract(90.0, decspace)
    decspace = np.radians(decspace)
    listspace = np.asarray(list(product(raspace, decspace)), dtype = float)
    lowcut = hp.pixelfunc.ang2pix(NSIDE, listspace[:,1], listspace[:,0])

    data = np.append(bad_pixels, pixcut)
    data = np.append(data, lowcut)
    data = data.astype(int)
    data = np.unique(data)

    mapi = np.arange(hp.nside2npix(NSIDE))
    mapi.fill(1)
    mapi[data] = 0
    hp.write_map("%s/bads_%s.fits" % (outdir_3, mag_cut), mapi, fits_IDL = False, overwrite = True)
    print("Visualization and producing Healpy map...")

    fig = plt.figure(1)
    hp.mollview(mapi, title = "badpixels")
    plt.savefig("%s/vis_bads_21.pdf" % outdir_3)


if __name__ == "__main__":

    resolution = 0.01 #adjusted for NSIDE = 2048
    NSIDE = 2048
    dust_path = "../dust_correction/dust_catalog.csv"
    cut_path = "sot_maglim_surhud.fit"

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--part", help = "Either 1,2 or 3")
    parser.add_argument("--outfile_1", help = "Output path for part 1")
    parser.add_argument("--outdir_2", help = "Output directory for part 2")
    parser.add_argument("--infile_3", help = "Input path for part 3")
    parser.add_argument("--outdir_3", help = "Output directory for part 3")
    parser.add_argument("--mag_cut", help = "Magnitude cut")

    args = parser.parse_args()
    part = args.part
    outfile_1 = args.outfile_1
    outdir_2 = args.outdir_2
    infile_3 = args.infile_3
    outdir_3 = args.outdir_3
    mag_cut = args.mag_cut


    if part == 2:
        from mpi4py import MPI
        from glob import glob
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size

    if part == 1:
        part1()
    elif part == 2:
        part2(rank)
        comm.Barrier()
    elif part == 3:
        part3()
    else:
        print("Part parameter wrong!")
