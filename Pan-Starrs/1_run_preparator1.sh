#PBS -lnodes=25:ppn=28,walltime=10:00:00
#PBS -N preparator
#PBS -q large
cd /home/dominik.zuercher/Documents/RSP_Pro/Pan-Starrs
module load mpi/mpich-x86_64
echo "Started!!!"
mpirun python PS_preparator.py --step 1 --mag_lim 22 --indir "/work/dominik.zuercher/DataStore/Pan-Starrs/original/split_deep" --outdir "/work/dominik.zuercher/DataStore/Pan-Starrs/PS_parts1/parts" >output1.txt 2>&1
