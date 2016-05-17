#!/bin/bash
#SBATCH -J test_re                    # job name
#SBATCH -o test_re.o%j                # output and error file name (%j expands to jobID)
#SBATCH -N 1                          # number of nodes requested
#SBATCH -n 24                         # total number of mpi tasks requested
#SBATCH -p development                # queue (partition) -- normal, development, etc.
#SBATCH -t 0:30:00                    # run time (hh:mm:ss) - 1.5 hours

#load some modules needed by Gromacs
module load boost
module load cxx11
module load gromacs/5.0.6

#start the runscript
./runscript.sh > output

