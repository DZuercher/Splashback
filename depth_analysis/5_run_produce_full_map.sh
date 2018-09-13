#PBS -lnodes=1:ppn=28,walltime=48:00:00
#PBS -N Separator
#PBS -q tiny
cd /home/dominik.zuercher/Documents/RSP_Pro/depth_analysis
module load mpi/mpich-x86_64
mpirun python produce_full_map.py --mid_directory "mask_pix_parts" >output.txt 2>&1
