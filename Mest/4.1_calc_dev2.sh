#PBS -lnodes=35:ppn=28,walltime=48:00:00
#PBS -N calc_stats_1
#PBS -q large
cd /home/dominik.zuercher/Documents/Splashback/Mest
module load mpi/mpich-x86_64
mpirun python mcmc_calc_dev2.py --type_ "Planck_PS_21.5_red_hard_spline" --add "_best"  > /work/dominik.zuercher/Output/Mest/logs/output_calc_red.log 2>&1
