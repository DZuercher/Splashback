#PBS -lnodes=35:ppn=28,walltime=10:00:00
#PBS -N random_generator
#PBS -q large
cd /home/dominik.zuercher/Documents/RSP_Pro/random_generate/PS_random
module load mpi/mpich-x86_64
echo "Started!!!" 
mpirun python random_generator.py --init 0 --mag_cut 21 --indir "/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms/finished" --middir "/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms2" --outdir "/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/randoms3" >output.txt 2>&1
