#PBS -lnodes=30:ppn=28,walltime=10:00:00
#PBS -N preparator
#PBS -q large
cd /home/dominik.zuercher/Documents/RSP_Pro/Pan-Starrs
module load mpi/mpich-x86_64
echo "Started!!!"
mpirun python PS_preparator.py --step 2 --mag_lim 21 --indir "/work/dominik.zuercher/DataStore/Pan-Starrs/PS_parts1/finished" --outdir "/work/dominik.zuercher/DataStore/Pan-Starrs/PS_parts2/parts" >output2.txt 2>&1
