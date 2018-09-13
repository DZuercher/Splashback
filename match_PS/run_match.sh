#!/bin/bash
#PBS -lnodes=35:ppn=28,walltime=10:00:00
#PBS -N match_PS
#PBS -q large
cd /home/dominik.zuercher/Documents/RSP_Pro/match_PS
pwd
module load mpi/mpich-3.2-x86_64

mpirun python match_spec.py > /work/dominik.zuercher/Output/match_PS/logs/output_spec.log 2>&1

