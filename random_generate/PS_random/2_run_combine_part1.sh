#PBS -lnodes=1:ppn=28,walltime=48:00:00
#PBS -N Random combiner
#PBS -q small
cd /home/dominik.zuercher/Documents/RSP_Pro/random_generate/PS_random/
module load mpi
mpirun ./combine_part1.sh >output.log 2>&1

