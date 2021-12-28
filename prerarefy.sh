#!/bin/bash

#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=36        # the number of CPU cores per node
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=48:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --account=uni

#SBATCH --job-name=RarefyDataA         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=kleineba@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED

module load palma/2020a
module load foss
module load R
...
# START THE APPLICATION
Rscript /home/k/kleineba/prerarefy.R /scratch/tmp/kleineba/PreRarefaction/Data/B _1 _2 21092


...