#PBS -lnodes=35:ppn=28,walltime=10:00:00
#PBS -N random_generator
#PBS -q large
cd /home/dominik.zuercher/Documents/RSP_Pro/random_generate/PS_random
module load mpi/mpich-x86_64
echo "Started!!!" 
mpirun python random_generator.py --init 1 --mag_cut 21 --outdir "/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms/parts" >output.txt 2>&1
