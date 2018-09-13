import sys
print(sys.path)
sys.path.insert(0, "install/lib64/python2.7/site-packages/")
print(sys.path)
import splashback
print("Look here:", splashback.__file__)
import frogress
import pandas
import numpy as np
import sys
import pyfits
import fitsio
import argparse
import yaml
from splashback_select import lens_select, source_select
from subprocess import call

def run_pipeline(config):

    rmin = config["rmin"]
    rmax = config["rmax"]
    rbin = config["rbin"]
    dirout = config["dirout"]
    randomsources = config["randomsources"]
    colored = config["colored"]
    deproject=config["deproject"]

    if randomsources:
        pairout = config["dirout"]+"/random_pairs.dat"
    else:
        pairout = config["dirout"]+"/pairs.dat"

    # No race condition here as the file is written at the end
    if not randomsources:
        call("mkdir -p %s" % (config["dirout"]), shell=1)
    else:
        # Wait until directory is written
        import os
        import time
        while not os.path.exists("./%s/" % (config["dirout"])):
            time.sleep(1)
    
    mag_limit = float(config["source"]["mag_limit"])
    z_limit = config["source"]["z_limit"]
    
    splashpipe = splashback.splashback(rmin, rmax, rbin, pairout, mag_limit, z_limit, config["lens"]["Njack"], deproject=deproject, colored=colored, verbose=False)
    
    # First get all the lens galaxy data
    print(config["lens"])
    ra, dec, zred, wt, jackreg = lens_select(config["lens"])

    #counting jackknife regions
    if "rank" in config["lens"]:
        fp=open(config["dirout"]+"/Ncluster_"+str(config["lens"]["rank"])+".dat","w+")
    else:
        fp = open(config["dirout"]+"/Ncluster.dat", "w")
    for ii in range(config["lens"]["Njack"]):
        idx = (ii != jackreg)
        fp.write("%d\n" % (np.sum(idx)))
    fp.close()
    splashpipe.allocate_lens_memory(ra.size)
    for i in range(ra.size):
        splashpipe.process_lens(ra[i], dec[i], zred[i], jackreg[i])
    tree = splashpipe.finalize_lenses()
    
    # Ready to pounce on the source data, a million galaxies at a time
    itern=0
    done = False
    sourceargs = config["source"]
    sourceargs["iter"] = 0
    if randomsources:
        sourceargs['type'] = sourceargs['type']+"_random"
    
    while not done:
        itern = itern + 1
        datagal, sourceargs, Ngal, status = source_select(sourceargs, 10000000)
        print(status)
        if status!=0 :
            break
    
        #if itern==10:
        #    break
    
        # For every source, query the clusters around it
        print("\n", sourceargs["iter"])
        #for i in frogress.bar(range(Ngal)):
        for i in range(Ngal):
        #for i in range(10):
    
            if sourceargs['filetype'] == "fits":
                if randomsources:
                    ragal, decgal = datagal[i]
                    maggal = 20.0
                else:
                    ragal, decgal, maggal = datagal[i]
            elif sourceargs['filetype'] == "ascii":
                if randomsources:
                    ragal, decgal = datagal.ra.values[i], datagal.dec.values[i]
                    maggal = 20.0
                else:
                    if colored!=0:
                        ragal, decgal, rband, gband, maggal = datagal.ra.values[i], datagal.dec.values[i], datagal.rband.values[i], datagal.gband.values[i], datagal.mag_auto.values[i]
                    else:
                        ragal, decgal, maggal = datagal.ra.values[i], datagal.dec.values[i], datagal.mag_auto.values[i]
            if colored!=0:
                if (gband<-100.0) | (rband<-100.0):
                    continue
                color=gband-rband
                splashpipe.process_source(ragal, decgal, maggal, randomsources, color)
            else:
                splashpipe.process_source(ragal, decgal, maggal, randomsources)
        splashpipe.finalize_results(writeok=False)
    
    splashpipe.finalize_results(writeok=True)
    return splashpipe

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", help="Configuration file")
    parser.add_argument("--rmax", help="Max radius", type=float, default=10.0)
    parser.add_argument("--rmin", help="Min radius", type=float, default=0.1)
    parser.add_argument("--rbin", help="Radial bin", type=int, default=8)
    parser.add_argument("--dirout", help="Output filename with pairs information", default="debug")
    parser.add_argument("--random", help="To use random sources or not", type=int, default=0)
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as ymlfile:
        config = yaml.load(ymlfile)
    print(config)
    print(args)

    # Convert arguments to variables
    rmax = (args.rmax)
    rmin = (args.rmin)
    rbin = (args.rbin)

    config["rmin"] = rmin
    config["rmax"] = rmax
    config["rbin"] = rbin
    config["dirout"] = args.dirout
    config["randomsources"] = args.random==1
    
    run_pipeline(config)
