#Reads original Planck cluster catalog and gets Pan-STARRS picture-names
# coding=utf-8
import pandas as pd
import sys
sys.path.insert(0,'../toolbox')
import tools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits


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

def convert_ra_from_list(listobj):
    hd = float(listobj[0])
    mts = float(listobj[1])
    sec = float(listobj[2])
    return ra2deg(hd, mts, sec)


def convert_dec_from_list(listobj):
    #sign = float(spl[0][0:1]+"1")
    sign = np.sign(float(listobj[0]+"1"))
    hd = float(listobj[0][1:])
    mts = float(listobj[1])
    sec = float(listobj[2])
    return dec2deg(sign, hd, mts, sec)


def read_Planck(inputfile, output_full):
    #Reading FITS data into npy table
    print("Reading data from original catalog")
    hdulist = fits.open(inputfile)
    head = hdulist[0].header
    body = hdulist[1].data
    full = pd.DataFrame(body)
    full.to_csv(output_full, index=False)
    cols = hdulist[1].columns
    data = np.zeros((1,8))
    length = body['INDEX'][-1]

    ra = np.asarray(body['RA'])
    dec = np.asarray(body['DEC'])
    redshift = np.asarray(body['REDSHIFT'])
    mass = np.asarray(body['MSZ'])
    return ra, dec, redshift, mass

def read_MCXC(inputfile, output_full):
    print("Reading data from MCXC catalog")

    #Convert ra/dec from hourangles to degrees
    with open(inputfile) as f:
	for i,line in enumerate(f):
	    if i<=4:
		continue
	    linelist = line.split("  ")
	    linelist = filter(None, linelist)
	    radec = linelist[1]
	    radeclist = radec.split(" ")
	    radeclist = filter(None, radeclist)
	    ra = radeclist[0:3]
	    dec = radeclist[3:6]

	    newra = convert_ra_from_list(ra)
	    newdec = convert_dec_from_list(dec)

	    radius = linelist[3]
	    name = linelist[0]
	    datlist = linelist[2].split(" ")
	    datlist = filter(None, datlist)
	    redshift = datlist[0]
	    #lx = datlist[1]
	    mass = datlist[2]
	    #rastring = "%:%:%" % (ra[0], ra[1], ra[2]) )
	    #decstring = "%:%:%" % (dec[0], dec[1], dec[2]) )

	    adder = pd.DataFrame([[name,newra,newdec,redshift,lx,mass,radius[:-2]]], columns = ["NAME","RA","DEC","REDSHIFT","lx_500","MSZ","radius_500"])
	    outframe = outframe.append(adder)
    outframe.to_csv(output_full, index=False)

    return np.asarray(outframe['RA'].values), np.asarray(outframe['DEC'].values), np.asarray(outframe['REDSHIFT'].values), np.asarray(outframe['MSZ'].values)

def generate_outputfile(ra, dec, redshift, mass):
    print("Collecting filenames...")
    data, goods = tools.get_filenames(ra=ra, dec=dec, filter="any")

    ridx = (data[:,4] == 'r')
    gidx = (data[:,4] == 'g')
    iidx = (data[:,4] == 'i')
    yidx = (data[:,4] == 'y')
    zidx = (data[:,4] == 'z')

    rdata = data[ridx,:]
    base = rdata[:,:-2]
    rdata = rdata[:,-2:].astype(str)
    gdata = data[gidx,-2:].astype(str)
    idata = data[iidx,-2:].astype(str)
    ydata = data[yidx,-2:].astype(str)
    zdata = data[zidx,-2:].astype(str)

    rdata[:,1] = [x[:-1] for x in rdata[:,1]]
    gdata[:,1] = [x[:-1] for x in gdata[:,1]]
    idata[:,1] = [x[:-1] for x in idata[:,1]]
    zdata[:,1] = [x[:-1] for x in zdata[:,1]]
    ydata[:,1] = [x[:-1] for x in ydata[:,1]]

    ra = ra[goods]
    dec = dec[goods]
    redshift = redshift[goods]
    mass = mass[goods]

    ra = ra.reshape((ra.size, 1))
    dec = dec.reshape((dec.size, 1))
    redshift = redshift.reshape((redshift.size, 1))
    mass = mass.reshape((mass.size, 1))
    base = np.delete(base, [2,3,4,5], axis=1)

    outdata = np.hstack((ra, dec, redshift, mass, base, rdata, gdata, idata, ydata, zdata))
    outdata = outdata.astype(str)

    print("Writing output-file")
    out = pd.DataFrame(outdata)
    out.to_csv(outfile, header=['RA', 'DEC', 'REDSHIFT', 'MASS', 'PRO', 'SKY', 'TYPE', 'RNAME', 'RSHORTNAME', 'GNAME', 'GSHORTNAME', 'INAME', 'ISHORTNAME', 'YNAME', 'YSHORTNAME', 'ZNAME', 'ZSHORTNAME'], index=False)



def print_mass_redshift(outfile):
    print("Generating Mass vs. redshift diagram...")
    inputdata = pd.read_csv(outfile, dtype=str)

    x = inputdata['REDSHIFT'].values
    y = inputdata['MASS'].values

    fig = plt.figure(1)
    plt.plot(x, y, 'o')
    plt.title("Halo Mass vs. Redshift")
    plt.xlabel('z')
    plt.xlim([0,0.4])
    plt.ylabel("MSZ")
    fig.savefig()

def write_files(outfile):
    print("Writing name-files")
    inputdata = pd.read_csv(outfile, dtype=str)
    #Files which are used by download.sh to get the pictures of the clusters
    ra = inputdata['RA'].values
    dec = inputdata['DEC'].values
    redshift = x

    color1 = inputdata['INAME'].values
    color2 = inputdata['RNAME'].values
    color3 = inputdata['GNAME'].values

    names = tools.generate_names(ra=ra, dec=dec, redshift=redshift)
    https = tools.generate_https(color1, color2, color3, ra=ra, dec=dec, redshift=redshift)

    np.savetxt(http_file, https, fmt='%s')
    np.savetxt(name_file, names, fmt='%s')


if __name__ == "__main__":

    catalog_type = "Planck" #Either 'Planck' or 'MCXC'


    inputfile = 'HFI_PCCS_SZ-union_R2.08.fits' #Path to input catalog in fits format
    output_full = 'Full_Planck.csv' #Generated csv file corresponding to original catalog (only used if Planck catalog is read)
    outfile = 'configfile.csv' #Path to csv file that gets generated and which
    #will contain the selected clusters with the corresponding picture names
    diagram = 'Mvsz.pdf' #Path to generated plot (pdf) for Mass vs. Redshift relation
    http_file = "PS_https.dat" #Path to generated file containing the addresses
    #of the pictures
    name_file = "PS_names.dat" #Path to generated file containing the more 
    #conveniant names which will be used for the downloaded pictures


    if catalog_type=="Planck":
	ra, dec, redshift, mass = read_Planck(inputfile, output_full):
    elif catalog_type=="MCXC":
	ra, dec, redshift, mass = read_MCXC(inputfile, output_full):
    else:
	print("catalog_type not available")

    generate_outputfile(ra, dec, redshift, mass)
    print_mass_redshift(outfile)
    write_files(outfile)
    print("Done")
