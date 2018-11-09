#!/bin/bash

#PBS -lnodes=35:ppn=28,walltime=10:00:00
#PBS -N Planck_PS
#PBS -q large

cd /home/dominik.zuercher/Documents/Splashback/splashpipe
pwd
module load mpi/mpich-3.2-x86_64


mpirun python mpi_splashback.py --rmax 50.0 --rmin 0.1 --rbin 15 --colored 1 --clusterdir "Planck" --galaxydir "Pan-Starrs_chunks" --magcut "21.5" --dirout "/work/dominik.zuercher/Output/splashpipe/Planck_PS_21.5_red_3sigma/Clu=Norm_Gal=Norm" > /work/dominik.zuercher/Output/splashpipe/logs/output_1.dat 2>&1
