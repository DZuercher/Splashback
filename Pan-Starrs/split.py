import sys

Catalog = "/work/dominik.zuercher/DataStore/Pan-Starrs/Catalogs/PS_catalog_22.csv"

Outdir = "/work/dominik.zuercher/DataStore/Pan-Starrs/Chunked_galaxies/PS_catalog_22/PS_catalog.csv"
lines_per_file = 80972 #Number of lines in Catalog/980 and roundup

with open(Catalog, "r") as f:
    itern = 0
    fnum = 0
    for line in f:
        if itern==0:
            if fnum!=0:
                fw.close()
            print "Writing fnum", fnum
            fw = open("%s" % Outdir + "%03d" % fnum,"w")
            fnum = fnum + 1
        fw.write(line)
        itern = itern + 1
        itern = itern%lines_per_file
