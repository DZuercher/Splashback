import sys
print(sys.path)
sys.path.insert(0, "install/lib64/python2.7/site-packages/")
print(sys.path)
import numpy as np
import socket
from mpi4py import MPI
import mpi4py
from splashback_aroundsource import run_pipeline
import argparse
import sys

def getsource_dict(rank, size, magcut, colored, galaxydir):
    d = {}
    d["type"] = galaxydir
    d["filetype"] = "ascii"
    d["mag_limit"] = magcut
    d["z_limit"] = 0.33
    d["rank"] = rank
    d["size"] = size
    d["colored"] = colored
    return d

def getlens_dict(rank, size, magcut, Njack, clusterdir):
    d = {}
    d["type"] = clusterdir
    d["zmin"] = 0.03
    d["zmax"] = 0.33
    d["Njack"] = Njack
    d["mag_limit"] = magcut
    d["rank"] = rank
    d["size"] = size
    return d

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--rmax", help="Max radius", type=float, default=10.0)
    parser.add_argument("--rmin", help="Min radius", type=float, default=0.1)
    parser.add_argument("--rbin", help="Radial bin", type=int, default=8)
    parser.add_argument("--dirout", help="Directory", default="debug")
    parser.add_argument("--random", help="To use random sources or not", type=int, default=0)
    parser.add_argument("--colored", help="To use color code or not (0,1,-1)", type=int, default=0)
    parser.add_argument("--magcut", help="Magnitude cut", default="21.5")
    parser.add_argument("--clusterdir", help="Cluster directory", default="Planck")
    parser.add_argument("--galaxydir", help="Galaxy directory", default="Pan-Starrs_chunks")
    parser.add_argument("--deproject", help="Deproject or not (True or False)", type=bool, default=False)
    parser.add_argument("--Njack", help="Number of Jackknife bins to use", default=30)

    args = parser.parse_args()

    ssample = getsource_dict(rank, size, args.magcut, args.colored, args.galaxydir)
    lsample = getlens_dict(rank, size, args.magcut, args.Njack, args.clusterdir)

    config = {}
    config["rmin"] = args.rmin
    config["rmax"] = args.rmax
    config["rbin"] = args.rbin
    config["dirout"] = args.dirout+"/%03d" % rank
    config["randomsources"] = (args.random == 1)
    if args.random == 1:
        galaxydir += "_random"
    config["colored"] = args.colored
    config["deproject"] = args.deproject

    config["lens"] = lsample
    config["source"] = ssample

    print(rank, config)
    sys.stderr.write("Rank: %d \n" % rank+socket.gethostname())

    a=run_pipeline(config)
    comm.Barrier()
    sys.stderr.write("We are at barrier : 1")
    comm.Barrier()
    sys.stderr.write("We are at barrier : 2")
