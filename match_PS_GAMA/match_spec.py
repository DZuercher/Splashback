#Matches RedMaPPer spectroscopic galaxies with PS galaxies. Contains redshift (from RedMaPPer) and color (from PS)
from astropy.io import fits
import scipy.spatial as kdtree
import sys
import numpy as np
from astropy.io import fits as pyfits
from astroML.crossmatch import crossmatch_angular
import pandas
from mpi4py import MPI


def get_assigned_targets(cat, ps_cat):
    ra_e = cat[:,0]
    dec_e = cat[:,1]
    X_e = np.empty((ra_e.size, 2), dtype = np.float64)
    X_e[:, 0] = ra_e
    X_e[:, 1] = dec_e

    data_m = ps_cat

    ra_m = data_m[:,0]
    dec_m = data_m[:,1]
    X_m = np.empty((ra_m.size, 2), dtype = np.float64)
    X_m[:, 0] = ra_m
    X_m[:, 1] = dec_m
    print "Out of:", ra_m.size

    max_radius = 0.3 / 3600  # 0.01 arcsec
    dist, ind = crossmatch_angular(X_m, X_e, max_radius)
    idx = (~np.isinf(dist))

    return data_m[idx], dist, ind


def with_kdtree(cat, ps_cat):
    #Prepare cat
    ra_e = cat[:,0]*np.pi/180.
    dec_e = cat[:,1]*np.pi/180.
    X_e = np.empty((ra_e.size, 3), dtype = np.float64)
    X_e[:, 0] = np.cos(ra_e)
    X_e[:, 1] = np.sin(dec_e)
    X_e[:, 2] = np.cos(dec_e)

    #Make kdtree
    tree = kdtree.KDTree(X_e)

    #Prepare ps_cat
    data_m = ps_cat
    ra_m = data_m[:,0]*np.pi/180.
    dec_m = data_m[:,1]*np.pi/180.
    X_m = np.empty((ra_m.size, 3), dtype = np.float64)
    X_m[:, 0] = np.cos(ra_m)
    X_m[:, 1] = np.sin(dec_m)
    X_m[:, 2] = np.cos(dec_m)
 
    max_radius = 0.3/3600*np.pi/180.0 #In arcsecs
    dd, ii = tree.query(X_m, distance_upper_bound = max_radius) #dd = Distance between matched objects (inf if no match found), ii == Index of matched pair in cat, botth have length of ps_cat 
    id_red = (~np.isinf(dd)) #Is true if distance dd is not infinite ( == matched object)
    #idloc=np.where(np.isnan(dd)==False)
    return dd, ii, id_red



def procedure(rank):
    print("started")
    #read spectroscopic GAMA galaxies
    hdul = fits.open(catalog_1)
    data = hdul[1].data
    ra = data['RA']
    dec = data['DEC']
    redshift = data['Z']

    idx = (redshift > 0.03) & (redshift < 0.33)
    ra = ra[idx]
    dec = dec[idx]
    redshift = redshift[idx]
    ra = ra.reshape(ra.size,1)  
    dec = dec.reshape(dec.size,1)  
    redshift = redshift.reshape(redshift.size,1)  
    cat = np.hstack((ra, dec, redshift))

    print(cat)

    #read PS 21.5 galaxies
    ps_dat = pandas.read_csv("%s/PS_catalog.csv%03d" % (catalog_2, rank), sep = ',', header = None, usecols = (1, 2, 5, 6, 7), names = (["ra", "dec", 'rband', 'gband', "mag_auto"]) )
    ra = ps_dat.ra.values
    dec = ps_dat.dec.values
    rband = ps_dat.rband.values
    gband = ps_dat.gband.values
    iband = ps_dat.mag_auto.values
    color = gband - rband
    ps_cat = np.hstack((ra.reshape(ra.size,1), dec.reshape(ra.size,1), color.reshape(color.size,1)))
    print("PS read")

    dd, ii, id_ = with_kdtree(cat, ps_cat)
    print("Matching done!")
    z_matched = np.zeros(0)
    iband_PS_matched = np.zeros(0)
    red_PS_matched = np.zeros(0)
    green_PS_matched = np.zeros(0)
    for el in range(dd.size):
	if id_[el] == True:
	    z_now = redshift[ii[el]]
            red_PS_now = rband[el]
            green_PS_now = gband[el]
            iband_PS_now = iband[el]
	    z_matched = np.append(z_matched, z_now)
            red_PS_matched = np.append(red_PS_matched, red_PS_now)
            green_PS_matched = np.append(green_PS_matched, green_PS_now)
            iband_PS_matched = np.append(iband_PS_matched, iband_PS_now)


    if z_matched.size > 0:
	np.savetxt("%s/matched_%03d" % (output_dir, rank), np.hstack((z_matched.reshape(z_matched.size,1), red_PS_matched.reshape(red_PS_matched.size,1), green_PS_matched.reshape(green_PS_matched.size,1), iband_PS_matched.reshape(iband_PS_matched.size,1) )))
	print("saved!")
    return



if __name__ == "__main__":

    output_dir = "/work/dominik.zuercher/Output/match_PS_GAMA/spec_parts"
    catalog_1= "/work/dominik.zuercher/DataStore/GAMA/SpecObj.fits" #RedMaPPer
    catalog_2 = "/work/dominik.zuercher/DataStore/Pan-Starrs/Chunked_galaxies/PS_catalog_21.5" # Chunked Pan-Starrs

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    procedure(rank)
    comm.Barrier()
    sys.stderr.write("We finished!!!")
