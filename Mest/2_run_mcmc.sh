#PBS -lnodes=1:ppn=28,walltime=48:00:00
#PBS -N MCMC_1
#PBS -q mini
cd /home/dominik.zuercher/Documents/Splashback/Mest
module load mpi/mpich-x86_64
mpirun python mcmc_sampler.py --steps 100000 --prevrun 0 --size 28 --type_ "Planck_PS_21.5_blue_hard_spline" --add "_best"  >/work/dominik.zuercher/Output/Mest/logs/output_blue.log 2>&1

