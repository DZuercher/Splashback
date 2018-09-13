#Requires a splitted version of the PS catalog
import numpy as np
import pandas as pd
import healpy as hp
from operator import itemgetter, attrgetter, methodcaller
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import os
from glob import glob
import argparse


def preparation1(rank):
    print("Reading original PS catalog...")
    catalog=pd.read_csv("%s/split_deep_part_%s.csv" % (catalog_dir, rank), delimiter = ',', names=['id','pro','sky','ira','idec','iPSFMag','iKronMag','SDSS_i_avg'], usecols = [0,1,2,5,7,9,11,16], header = None, dtype = {'id':str})
    print("Throw away stars...")

    catalog = catalog.loc[catalog.pro != -999.0]
    catalog = catalog.loc[catalog.sky != -999.0]
    catalog = catalog.loc[catalog.ira != -999.0]
    catalog = catalog.loc[catalog.idec != -999.0]
    catalog = catalog.loc[catalog.iPSFMag != -999.0]
    catalog = catalog.loc[catalog.iKronMag != -999.0]
    catalog = catalog.loc[catalog.SDSS_i_avg != -999.0]
    
    try:
        catalog['mag'] = catalog["iKronMag"].values - catalog["SDSS_i_avg"].values
    except:
        print("Broken " + str(rank))
        print(catalog["SDSS_i_avg"])
    catalog['iPSFMag'] = np.asarray(catalog["iPSFMag"].values) - np.asarray(catalog["SDSS_i_avg"].values)

    newcat = catalog[((catalog["iPSFMag"] - catalog['mag']) > 0.05) & (catalog['iPSFMag'] > 15.0) & (catalog['mag'] <= mag_lim)]
    #For star catalog production
    #newcat = catalog[((catalog["iPSFMag"] - catalog['mag']) < 0.05) & (catalog['mag'] <= 21.5)]

    #Get rid of unnecessary columns
    del newcat['iKronMag']
    del newcat['iPSFMag']
    del newcat['SDSS_i_avg']

    print("Throw away objects in bad skycells (low mag limit or no mask available)...")
    bads = hp.read_map("%s/bads_%s.fits" % (mask_path, mag_lim))

    ra = np.asarray(newcat["ira"].values)
    dec = np.asarray(newcat["idec"].values)
            
    ra = np.radians(ra)
    dec = np.radians(np.subtract(90.0, dec))

    print("Performing Healpy selection...")
    #Perform selection using Healpy
    pix = hp.pixelfunc.ang2pix(NSIDE, dec, ra)
    pixmask = (bads[pix] == 1)
    newcat2 = newcat[pixmask]

    ra = newcat2["ira"].values
    dec = newcat2["idec"].values
    pros = tools.decide_pro(ra, dec)
    newcat2.loc[:,"selfpro"] = pd.Series(pros, index = newcat2.index)

    for i in range(635, 2644):
        out = newcat2.loc[newcat2.selfpro.values.astype(int) == int(i)]
        out.to_csv("%s/PS_cell_%s_rank_%s.dat" % (outdir_1, i, rank), index = False, header = False, sep = ' ')


def preparation2(rank):
	
    steps = 2643.0 - 635.0

    #Good splitter
    chunky = int(steps/size)
    rest = steps - chunky*size
    mini = chunky*rank
    maxi = chunky*(rank + 1)
    if rank >= (size - 1) - rest:
        maxi += 2 + rank - size + rest)
        mini += rank - size + 1 + rest
    mini += 635
    maxi += 635
    mini = int(mini)
    maxi = int(maxi)

    for cellnumber in np.arange(mini, maxi, 1):
        try:
            catalog = pd.read_csv("%s/PS_cell_%s.dat" % (indir_2, cellnumber), sep = ' ', header = None, names = ['id','pro','sky','ira','idec','mag','selfpro'], dtype = {'id':str})
        except:
            continue
        procell = np.asarray(catalog["pro"].values)
        if procell.size == 0:
            continue
        procell = procell.reshape((procell.size, 1)).astype(int)
        skycell = np.asarray(catalog['sky'].values).reshape((procell.size, 1)).astype(int)
        ra = np.asarray(catalog["ira"].values).reshape((procell.size, 1))
        dec = np.asarray(catalog["idec"].values).reshape((procell.size, 1))
        all_ids = np.asarray(catalog["id"].values).reshape((procell.size, 1))

        print("Deciding on Skycells...")		
        selfprocell = np.asarray(catalog['selfpro'].values).reshape((procell.size, 1)).astype(int)
        selfskycell = tools.decide_sky(ra, dec, cellnumber)
        selfskycell = selfskycell.reshape((procell.size, 1)).astype(int)
        catalog["selfsky"] = selfskycell
        
        del catalog['selfpro']
        del catalog['selfsky']

        procell = np.asarray(catalog["pro"].values)
        procell = procell.reshape((procell.size, 1)).astype(int)
        skycell = np.asarray(catalog['sky'].values).reshape((procell.size, 1)).astype(int)
        ra = np.asarray(catalog["ira"].values).reshape((procell.size, 1))
        dec = np.asarray(catalog["idec"].values).reshape((procell.size, 1))
        mag = np.asarray(catalog["mag"].values).reshape((procell.size, 1))
        id_ = np.asarray(catalog["id"].values).reshape((procell.size, 1))

        alginput = np.hstack((ra, dec, procell, skycell, mag))
        print("Shape %s and %s" % (alginput.shape, id_.shape))

        print("Performing masking...")
        newcat, newids = tools.mask_select(alginput, cellnumber, id_)
       
        print("Shape %s and %s" % (newcat.shape, newids.shape))
        newcat = np.hstack((newcat, newids))
        newcat = pd.DataFrame(newcat)
        newcat.to_csv("%s/PS_part_%s.csv" % (outdir_2, cellnumber), index = False, header = ['ra','dec','pro','sky','mag','id'])
		


if __name__ == "__main__":

    mask_path = "/home/dominik.zuercher/Documents/RSP_Pro/depth_analysis"
    NSIDE=2048

    from mpi4py import MPI
    comm=MPI.COMM_WORLD
    rank=comm.rank
    size=comm.size

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--step", help="Step of the process (1 or 2)")
    parser.add_argument("--indir", help="Input directory")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--mag_lim", help="Either 21, 21.5 or 22")
    args = parser.parse_args()

    mag_lim = args.mag_lim
    step = args.step
    indir = args.indir
    outdir = args.outdir

    if step == 1:
        preparation1(rank)
    if step == 2:
        preparation2(rank)
    comm.Barrier()

