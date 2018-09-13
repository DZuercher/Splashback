#!/bin/bash

#PBS -lnodes=35:ppn=28,walltime=10:00:00
#PBS -N Planck_PS
#PBS -q large

cd /home/dominik.zuercher/Documents/RSP_Pro/splashpipe
pwd
module load mpi/mpich-3.2-x86_64


mpirun python mpi_splashback.py --rmax 50.0 --rmin 0.1 --rbin 15 --colored -1 --clusterdir "Planck" --galaxydir "Pan-Starrs_chunks" --magcut "21.5" --dirout "/work/dominik.zuercher/Output/splashpipe/Planck_PS_21.5_blue_spline/Clu=Norm_Gal=Norm" > /work/dominik.zuercher/Output/splashpipe/logs/output_1.dat 2>&1
#mpirun python mpi_splashback.py --rmax 50.0 --rmin 0.1 --rbin 15 --colored -1 --clusterdir "Planck" --galaxydir "Pan-Starrs_chunks" --magcut "21.5" --dirout "/work/dominik.zuercher/Output/splashpipe/Planck_PS_21.5_blue_cut_9/Clu=Norm_Gal=Norm" > /work/dominik.zuercher/Output/splashpipe/logs/output_1.dat 2>&1


#Old ones left for doc purposes
#Deprojection
#mpirun python mpi_splashback_galaxies.py --rmax 10.0 --rmin 0.1 --rbin 8 --random 0 --colored 0 --deproject 1 --clusterdir "Planck" --magcut "21.5" --dirout "reloaded/Planck_PS_21.5_deprojected" > output.log 2>&1 
#Colored
#mpirun python mpi_splashback_galaxies.py --rmax 50.0 --rmin 0.1 --rbin 15 --random 0 --colored -1 --clusterdir "Planck_random" --magcut "21.5" --dirout "reloaded/Planck_PS_21.5_blue_random" > output.log 2>&1 
#mpirun python mpi_splashback_galaxies.py --rmax 50.0 --rmin 0.1 --rbin 15 --random 0 --colored -1 --clusterdir "Planck" --magcut "21.5" --dirout "debug_Planck_PS_21.5_blue" > output.log 2>&1 
#mpirun python mpi_splashback_galaxies.py --rmax 10.0 --rmin 0.1 --rbin 8 --clusterdir "Planck_random" --magcut "22" --dirout "reloaded/Planck_PS_22_random_reloaded" > output.log 2>&1
#mpirun python mpi_splashback_galaxies.py --rmax 10.0 --rmin 0.1 --rbin 8 --dirout "debug_zmin_0.33_magcut_21.5_sz"
#mpirun python mpi_splashback_galaxies.py --rmax 10.0 --rmin 0.1 --rbin 8 --dirout "debug_zmin_0.33_magcut_21.5_szrandoms"
#mpirun python mpi_splashback_galaxies_sdssxrmapper.py --rmax 10.0 --rmin 0.1 --rbin 8 --dirout "debug_zmin_0.33_magcut_21.0_sdssxrmapper_buggy"
#mpirun python mpi_splashback_galaxies_sdssxrmapper.py --rmax 10.0 --rmin 0.1 --rbin 8 --dirout "debug_zmin_0.33_magcut_21.0_sdssxrmapper_buggy" --random 1
