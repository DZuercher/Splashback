import numpy as np
import pyfits
import pandas
import sys
import fitsio
import glob

def lens_select(lensargs):
    np.random.seed(10)
    Njack = lensargs['Njack']
    rank = lensargs["rank"]
    size = lensargs["size"]

    if lensargs['type'] == "Planck":
        df = pandas.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/Catalogs/%02d_bins/Planck_%02d_bins.csv" % (Njack, Njack))
        ra = df["RA"].values
        dec = df["DEC"].values
        zred = df["REDSHIFT"].values
        wt = ra/ra
        jackreg = df["jackreg"].values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        jackreg = jackreg[idx]
        jackreg=jackreg.astype(int)
	sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra, dec, zred, wt, jackreg

    if (lensargs['type'] == "Planck_random"):
        df = pandas.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/Randoms/%02d_bins/Planck_randoms_%02d_bins_%s.csv" % (Njack, Njack, str(lensargs["mag_limit"])), sep = ",")
        ra = df["RA"].values
        dec = df["DEC"].values
        zred = df["REDSHIFT"].values
        wt = ra/ra
        jackreg = df["jackreg"].values
        #rankarr = np.random.randint(size, size=ra.size)
        #idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (rankarr==rank)
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        jackreg = jackreg[idx]
        jackreg=jackreg.astype(int)
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra, dec, zred, wt, jackreg

    if lensargs['type'] == "redmapper":
        hdulist = pyfits.open("/work/dominik.zuercher/DataStore/RedMaPPer/Originals/redmapper_dr8_public_v5.10_clusters.fits")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda"]
        wt = ra/ra

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        jackreg = np.random.randint(Njack, size=ra.size)
        return ra, dec, zred, jackreg

    #Further options are not enabled currently
    """
    if lensargs['type'] == "redmapper-random":
        ra, dec, zred, lamda, wt = np.loadtxt("/home/surhud/DataStore/redmapper/redmapper_public_v5.10_randoms_%05d.dat" % (lensargs["rannum"]), unpack=1)
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d randoms \n" % (np.sum(idx)))
        jackreg = np.random.randint(Njack, size=ra.size)
        return ra[idx], dec[idx], zred[idx], jackreg[idx]

    if lensargs['type'] == "camira":
        ra, dec, zred, lam = np.loadtxt("./Atsushi_data/camira_%s_wide_v%d.dat" % (lensargs['hsc-release'], lensargs['version']), unpack=1, usecols=(0, 1, 2, 3))
        wt = lam/lam

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lam>lensargs["lammin"]) & (lam<=lensargs["lammax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))

        jackreg = np.random.randint(Njack, size=ra.size)
        return ra[idx], dec[idx], zred[idx], jackreg[idx]

    if lensargs['type'] == "dr12mgII":
        hdulist = pyfits.open("DataStore/Trimmed_BOSS_DR12_107.fits")
        data = hdulist[1].data
        ra = data["RA"]
        dec = data["Dec"]
        zabs = data["ZABS"].astype('float64')
        err_zabs = data["ERR_ZABS"].astype('float64')
        zqso = data["ZQSO"].astype('float64')
        rew_mg27 = data["REW_MGII_2796"].astype('float64')
        rew_mg28 = data["REW_MGII_2803"].astype('float64')
        crit = data["CRITERION_MGII_FEII"]

        idx = (zabs>lensargs["zmin"]) & (zabs<=lensargs["zmax"]) & (zqso>=lensargs["zmax"]) & (rew_mg27>lensargs["ew_mgii_2796_min"]) & (crit>0.5)
        ra = ra[idx]
        dec = dec[idx]
        zabs = zabs[idx]
        err_zabs = err_zabs[idx]
        zqso = zqso[idx]
        rew_mg27 = rew_mg27[idx]
        rew_mg28 = rew_mg28[idx]
        sys.stdout.write("Selecting %d absorbers \n" % (np.sum(idx)))
        jackreg = np.random.randint(Njack, size=ra.size)
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zabs, jackreg

    if lensargs['type'] == "planck-sz-panstarrs":
        df = pandas.read_csv("./DataStore/Planck-SZ/%02d/Planck-SZ.csv" % Njack)
        ra = df["RA"].values
        dec = df["DEC"].values
        zred = df["REDSHIFT"].values
        wt = ra/ra
        jackreg = df["jackreg"].values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        jackreg = jackreg[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra, dec, zred, jackreg

    if lensargs['type'] == "planck-sz-panstarrs-random":
        df = pandas.read_csv("./DataStore/Planck-SZ/%02d/Planck_randoms.csv" % Njack)
        ra = df["ra"].values
        dec = df["dec"].values
        zred = df["redshift"].values
        wt = ra/ra
        jackreg = df["jackreg_ran"].values
        rankarr = np.random.randint(size, size=ra.size)

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (rankarr==rank)
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        jackreg = jackreg[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra, dec, zred, jackreg

    if lensargs['type'] == "planck-sz":
        df = pandas.read_csv("/home/surhud/PromisePegasus/DataStore/Planck/planck_sz_2015.filtered.dat", delim_whitespace=1)
        ra = df["ra"].values
        dec = df["dec"].values
        zred = df["zred"].values
        wt = ra/ra

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        jackreg = np.random.randint(Njack, size=ra.size)
        return ra, dec, zred, jackreg

    if lensargs['type'] == "redmapper-des-sva1":
        hdulist = pyfits.open("../redmapper_sva1_public_v6.3_catalog.fits")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda"]
        wt = ra/ra

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        jackreg = np.random.randint(Njack, size=ra.size)
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, jackreg
    """


def source_select(sourceargs, chunksize):


    if sourceargs['type'] == "Pan-Starrs" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Catalogs/PS_catalog_%s.csv" % (str(sourceargs["mag_limit"])), sep = ',', chunksize = chunksize, header = None, usecols=(1, 2, 7), names=(["ra", "dec", "mag_auto"]))
	try:
	    datagal = sourceargs['dfchunks'].next()
	    Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    if sourceargs['type'] == "Pan-Starrs_random" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/PS_randoms_%s.csv" % (str(sourceargs["mag_limit"])), sep = ',', usecols=[0, 1], chunksize = chunksize, header = None, names = (["ra", "dec"]))
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    if sourceargs['type'] == "Pan-Starrs_chunks" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Chunked_galaxies/PS_catalog_%s/PS_catalog.csv%03d" % (str(sourceargs["mag_limit"]), sourceargs["rank"]), chunksize = chunksize, header = None, usecols=(1, 2, 5, 6, 7), names=(["ra", "dec", "rband", "gband", "mag_auto"]), sep= ",")
	try:
	    datagal = sourceargs['dfchunks'].next()
	    Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1


    if sourceargs['type'] == "Pan-Starrs_chunks_random" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Chunked_randoms/PS_randoms_with_mag_%s/PS_randoms_%s.dat%03d.dat" % (str(sourceargs["mag_limit"]), str(sourceargs["mag_limit"]), sourceargs["rank"]), usecols = [0,1], chunksize = chunksize, header = None, names = (["ra", "dec"]), sep=",")
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    #Further options are currently not enabled
    """
    if sourceargs['type'] == "atsushi_hsc_s15b" and sourceargs['filetype'] == "fits":
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits'] = fitsio.FITS('./Atsushi_data/wide_s15B.fits')
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        if status:
            return datagal, sourceargs, Ngal, status

        wheremask = sourceargs['fits'][1].where("cmodel_i_corr < %.5f && i_countinputs >= %d && #row>%d && #row<=%d" % (sourceargs['mag_limit'], sourceargs['countinputs_limit'], itern*chunksize, (itern+1)*chunksize))
        try:
            datagal = sourceargs['fits'][1]['ra2000','decl2000','cmodel_i_corr'][wheremask]
        except:
            status = 1
            sourceargs['iter'] = sourceargs['iter'] + 1
            return datagal, sourceargs, Ngal, status

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]

    if sourceargs['type'] == "atsushi_hsc_s15b_random" and sourceargs['filetype'] == "fits":
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits'] = fitsio.FITS('./Atsushi_data/random.fits')
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        if status:
            return datagal, sourceargs, Ngal, status

        wheremask = sourceargs['fits'][1].where("i_countinputs >= %d && #row>%d && #row<=%d" % (sourceargs['countinputs_limit'], itern*chunksize, (itern+1)*chunksize))
        try:
            datagal = sourceargs['fits'][1]['ra2000','decl2000'][wheremask]
        except:
            status = 0
            sourceargs['iter'] = sourceargs['iter'] + 1
            return datagal, sourceargs, Ngal, status

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]

    if sourceargs['type'] == "des_sva1_gold_masked_z22.1" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("../sva1_gold_r1.0_catalog.masked.z22.1_rest.dat", chunksize=chunksize, delim_whitespace=1)
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1


    if sourceargs['type'] == "des_sva1_gold_masked_z22.1_random" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("../sva1_gold_r1.0_catalog.masked.z22.1_random_rest.dat", chunksize=chunksize, delim_whitespace=1, names=(["ra", "dec"]), header=None)
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    if sourceargs['type'] == "sdss_chunks" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("./DataStore/SDSS/Chunked_galaxies/cross-correlation.21i_wcol.filtered.dat%03d.dat" % (sourceargs["rank"]), chunksize=chunksize, delim_whitespace=1, usecols=(0, 1, 2), skiprows=1, header=None, names=(["ra", "dec", "mag_auto"]))
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    if sourceargs['type'] == "sdss_chunks_random" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("./DataStore/SDSS/Chunked_randoms/Randoms_sdss.dat%03d.dat" %(sourceargs["rank"]), chunksize=chunksize, delim_whitespace=1, skiprows=1, header=None, names=(["ra", "dec"]), usecols=(2, 3))
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    if sourceargs['type'] == "sdss" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("./DataStore/SDSS/cross-correlation.21i_wcol.filtered.dat", chunksize=chunksize, delim_whitespace=1, usecols=(0, 1, 2), skiprows=1, header=None, names=(["ra", "dec", "mag_auto"]))
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1

    if sourceargs['type'] == "sdss_random" and sourceargs['filetype'] == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("./DataStore/SDSS/Randoms_sdss.dat", chunksize=chunksize, delim_whitespace=1, skiprows=1, header=None, names=(["ra", "dec"]), usecols=(2, 3))
        try:
            datagal = sourceargs['dfchunks'].next()
            Ngal = datagal.ra.size
            status = 0
        except:
            datagal = 0
            Ngal = 0
            status = 1
        sourceargs['iter'] = sourceargs['iter'] + 1
    """
    return datagal, sourceargs, Ngal, status
