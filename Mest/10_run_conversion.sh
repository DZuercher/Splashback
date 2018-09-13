#PBS -lnodes=1:ppn=1,walltime=48:00:00
#PBS -N MCMC_1
#PBS -q tiny
cd /home/dominik.zuercher/Documents/RSP_Pro/Mest
python conversion_stats.py --modded 1 --type_ "Planck_PS_21_reloaded" --add "_best"  >output1.log 2>&1

