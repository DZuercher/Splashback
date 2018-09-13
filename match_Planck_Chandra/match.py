#Matches Planck clusters to Chandra analogs in order to get miscentering assuming Chandra clusters at correct posisitons
import numpy as np
from astropy.io import fits as pyfits
import pandas as pd

def dec2deg(sign, hd, min, sec):
    """
    Converts declination to degrees

    Parameters
    ----------
    hd : int
    degrees
    m : int
    arcminutes
    s : float
    arcseconds

    Returns
    -------
    hd : float
    A decimal number

    """
    return sign*(np.absolute(hd) + min/60.0 + sec/3600.0)

def ra2deg(hd, min, sec):
    """
    Converts ra to degrees

    Parameters
    ----------
    hd : int
    hours
    m : int
    minutes
    s : float
    seconds

    Returns
    -------
    hd : float
    A decimal number

        """
    return (np.float(hd) + min/60.0 + sec/3600.0)*360./24.

def convert_ra_from_list(chandrara):
    outframe=np.zeros(0)
    for obj in chandrara:
        splitra=obj.split(":")
	hd = float(splitra[0])
	mts = float(splitra[1])
	sec = float(splitra[2])
        outra = ra2deg(hd, mts, sec)
        outframe=np.append(outframe,outra)
    return outframe

def convert_dec_from_list(chandradec):
    outframe=np.zeros(0)
    for obj in chandradec:
        splitdec = obj.split(":")
	sign = np.sign(float(splitdec[0]+"1"))
	hd = float(splitdec[0])
	mts = float(splitdec[1])
	sec = float(splitdec[2])
        outdec = dec2deg(sign, hd, mts, sec)
        outframe=np.append(outframe,outdec)
    return outframe

def match(rra, ddec, xxra, xxdec,cut):
	rra = np.pi/180.*rra
    	ddec = np.pi/180.*ddec
    	xx = np.cos(ddec)*np.cos(rra)
    	yy = np.cos(ddec)*np.sin(rra)
    	zz = np.sin(ddec)

    	xxra = np.pi/180.*xxra
    	xxdec = np.pi/180.*xxdec
    	xxp = np.cos(xxdec)*np.cos(xxra)
    	yyp = np.cos(xxdec)*np.sin(xxra)
    	zzp = np.sin(xxdec)

    	from scipy.spatial import cKDTree

    	tree = cKDTree(zip(xxp, yyp, zzp))

    	dist, idx = tree.query(zip(xx, yy, zz), k=1)

    	# Convert distance to arcseconds from radians
    	dist = dist * 180./np.pi * 3600.0
	
	passidx = dist <= cut
	notpassedidx=np.invert(passidx)
	matched=idx[passidx]
	notmatched=idx[notpassedidx]
    	return matched, notmatched, passidx,dist


if __name__=="__main__":

    #Load catalogs
    print("Loading Planck...")
    planck=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/30/Planck_corrected.dat",sep=",")
    planckra=np.asarray(planck["RA"].values)
    planckdec=np.asarray(planck["DEC"].values)

    print("Loading Chandra...")
    chandra=pd.read_csv("/work/dominik.zuercher/DataStore/Chandra/Chandra_original.tab",delim_whitespace=True,header=None,skiprows=2,names=["Name","RA","Dec","z","K0","K100","alpha","Tcl","Lbol","UL1","LHa","UL2","Lrad"])
    chandrara=np.asarray(chandra["RA"].values)
    chandrara=convert_ra_from_list(chandrara)
    chandradec=np.asarray(chandra["Dec"].values)
    chandradec=convert_dec_from_list(chandradec)

    planck_beam = 7.5 #arcmin
    cut = planck_beam*60.0

    print("matching...")
    matched, notmatched, passidx, dist = match(planckra,planckdec,chandrara,chandradec,cut)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f=plt.figure(1)
    plt.hist(dist,100,range=(0.0,cut))
    plt.savefig("distances.pdf")

    print("Found " + str(matched.size) + " matching objects")

    matchlist = ["NO"]*planckra.size
    matchdist = [-999]*planckra.size
    for idx,entry in enumerate(passidx):
        if entry==True:
           matchlist[idx]="YES"
           matchdist[idx]=dist[idx]
    planck["Chandra_match"] = matchlist
    planck["Chandra_distance"] = matchdist
    planck.to_csv("/work/dominik.zuercher/DataStore/Planck-SZ/30/Planck_corrected_matched.csv")

