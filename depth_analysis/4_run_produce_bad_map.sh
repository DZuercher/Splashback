#PBS -lnodes=6:ppn=28,walltime=48:00:00
#PBS -N Separator
#PBS -q small
cd /home/dominik.zuercher/Documents/RSP_Pro/depth_analysis
module load mpi/mpich-x86_64
mpirun python produce_bad_map.py --part 3 --mag_cut 21 --infile_3 "bad_pixels.dat" --outdir_3 "." > output.txt 2>&1
