#Matches RedMaPPer spectroscopic galaxies with PS galaxies. Contains redshift (from RedMaPPer) and color (from PS)
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
    #read redmapper galaxies
    dat = pandas.read_csv(catalog_1, delim_whitespace = 1)
    idx = (dat.REDSHIFT.values > 0.03) & (dat.REDSHIFT.values < 0.33)
    ra = dat.RA.values[idx]
    dec = dat.DEC.values[idx]
    zred = dat.REDSHIFT.values[idx]
    ra = ra.reshape(ra.size,1)
    dec = dec.reshape(dec.size,1)
    zred = zred.reshape(zred.size,1)
    cat = np.hstack((ra, dec, zred))
    print("redmapper read")

    #read PS 21.5 galaxies
    ps_dat = pandas.read_csv("%s/PS_catalog.csv%03d" % (catalog_2, rank), sep = ',', header = None, usecols = (1, 2, 5, 6, 7), names = (["ra", "dec", 'rband', 'gband', "mag_auto"]) )
    ra = ps_dat.ra.values
    dec = ps_dat.dec.values
    rband = ps_dat.rband.values
    gband = ps_dat.gband.values
    color = gband - rband
    ps_cat = np.hstack((ra.reshape(ra.size,1), dec.reshape(ra.size,1), color.reshape(color.size,1)))
    print("PS read")

    #ps_cat_matched,dist,ind=get_assigned_targets(cat,ps_cat)
    dd, ii, id_red = with_kdtree(cat, ps_cat)
    print("Matching done!")
    z_red = np.zeros(0)
    color_red = np.zeros(0)
    z_blue = np.zeros(0)
    color_blue = np.zeros(0)
    for el in range(dd.size):
        #Problem: if dd = inf then in ii is just last entry in redmapper catalog...
        print(dd)
	if id_red[el] == True:
	    zred_now = zred[ii[el]]
	    color_now = color[el]
	    z_red = np.append(zred_out, zred_now)
	    color_red = np.append(color_out, color_now)
	else:
	    zred_now = zred[ii[el]]
	    color_now = color[el]
	    z_blue = np.append(zred_out, zred_now)
	    color_blue = np.append(color_out, color_now)

    if z_red.size > 0:
	np.savetxt("%s/reds_%03d" % (output_dir, rank), np.hstack((z_red.reshape(z_red.size,1), color_red.reshape(color_red.size,1))))
	print("saved!")
    if z_blue.size > 0:
	np.savetxt("%s/blues_%03d" % (output_dir, rank), np.hstack((z_blue.reshape(z_blue.size,1), color_blue.reshape(color_blue.size,1))))
	print("saved!")
    return



if __name__ == "__main__":

    output_dir = "/work/dominik.zuercher/Output/match_PS/parts"
    catalog_1= "/work/dominik.zuercher/DataStore/RedMaPPer/Catalogs/redmapper_galaxies.dat" #RedMaPPer
    catalog_2 = "/work/dominik.zuercher/DataStore/Pan-Starrs/Chunked_galaxies/PS_catalog_21.5" # Chunked Pan-Starrs

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    procedure(rank)
    comm.Barrier()
    sys.stderr.write("We finished!!!")
