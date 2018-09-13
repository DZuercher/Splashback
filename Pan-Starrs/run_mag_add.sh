#PBS -lnodes=35:ppn=28,walltime=10:00:00
#PBS -N 1_Planck_PS_redux
#PBS -q large

cd $PBS_O_WORKDIR
pwd
module load mpi/mpich-x86_64
mpirun python add_random_magnitude.py 22 > output.log 2>&1 
