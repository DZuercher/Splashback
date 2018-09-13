#Adds a random magnitude to the random catalogs (needed if chunked version of randoms is used)
import sys
import pandas as pd
import numpy as np
from mpi4py import MPI
comm=MPI.COMM_WORLD
rank=comm.rank
size=comm.size
mag_limit = sys.argv[1]

def add_mag(rank,mag_limit):
    print("Reading random chunk")
    dfrand = pd.read_csv("Chunked_randoms/PS_randoms_%s.dat%03d.dat"%(mag_limit,rank),sep=",",names=["ra","dec","pro","sky"])

    randsize=len(dfrand.ra.values)

    np.random.seed()
    chunkrandom=np.random.randint(0,980)

    print("Reading galaxy chunk")
    df = pd.read_csv("../Chunked_galaxies/PS_catalog_%s.csv%03d.dat"%(mag_limit,chunkrandom),names=["mag"],usecols=[4])

    np.random.seed()
    linerandoms=np.random.randint(0,len(df.mag.values),size=randsize)
    
    print("Drawing random magnitudes")
    randmags=df.mag[linerandoms]
    dfrand['mag']=randmags.values

    print("Saving")
    dfrand.to_csv("Chunked_randoms_magnitudes/PS_randoms_%s.dat%03d.dat"%(mag_limit,rank),sep=',',index=False,header=False)




########################
#MAIN
########################

add_mag(rank,mag_limit)

