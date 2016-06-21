#!/bin/bash
#SBATCH -J honeycomb_3c                    # job name
#SBATCH -o honeycomb_3c.o%j                # output and error file name (%j expands to jobID)
#SBATCH -N 1                          # number of nodes requested
#SBATCH -n 16                         # total number of mpi tasks requested
#SBATCH -p development                    # queue (partition) -- normal, development, etc.
#SBATCH -t 2:00:00                   # run time (hh:mm:ss) - 1.5 hours

#load some modules needed by Gromacs
module load intel/13.0.2.146  mvapich2/1.9a2
module load gromacs/4.6.5

#start the runscript
./runscript.sh > output