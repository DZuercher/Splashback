# Completes pictures located on edges of skycells
import os
import numpy as np
from glob import glob
import numpy as np
import sys
sys.path.insert(0,'../toolbox')
import tools
import frogress


if __name__=="__main__":

    indir = "./Pan-STARRS_Pics/original_pictures/" #Input directory containing downloaded pictures
    outdir = "./Pan-STARRS_Pics/edge_corrected/" #Output directory for corrected pictures
    conf_path = "configfile.csv" #Output file of the readin process
    it = 0 #Pictures to skip


    picnames = glob(indir+"*") 
    names, ras, decs, widths, heights, pros, skys, outs = tools.edge_check(picnames, conf_path)
    print("There are %s objects that need to be corrected" % names.size)

    np.savetxt(outdir + "names.txt", names, fmt = '%.s')
    np.savetxt(outdir + "ras.txt", ras)
    np.savetxt(outdir + "decs.txt", decs)
    np.savetxt(outdir + "widths.txt", widths)
    np.savetxt(outdir + "heights.txt", heights)
    np.savetxt(outdir + "pros.txt", pros)
    np.savetxt(outdir + "skys.txt", skys)
    np.savetxt(outdir + "outs.txt", outs)

    names = np.genfromtxt(outdir + "names.txt", dtype=str)
    ras = np.genfromtxt(outdir + "ras.txt")
    decs = np.genfromtxt(outdir + "decs.txt")
    widths = np.genfromtxt(outdir + "widths.txt")
    heights = np.genfromtxt(outdir + "heights.txt")
    pros = np.genfromtxt(outdir + "pros.txt")
    skys = np.genfromtxt(outdir + "skys.txt")
    outs = np.genfromtxt(outdir + "outs.txt")

    names = names[it:]
    ras = ras[it:]
    decs = decs[it:]
    widths = widths[it:]
    heights = heights[it:]
    pros = pros[it:]
    skys = skys[it:]
    outs = outs[it:]


    for j in range(len(names)):
            names[j] = str(outdir) + str(os.path.basename(names[j]))

    for i in frogress.bar(range(len(names))):
            tools.edge_correct(names[i],pros[i],skys[i],ras[i],decs[i],outs[i],widths[i],heights[i],indir,outdir)


