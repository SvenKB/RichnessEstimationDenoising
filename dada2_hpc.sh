#!/bin/bash

#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=36        # the number of CPU cores per node
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=2-00:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --account=uni
#SBATCH --mem=192000

#SBATCH --job-name=Data_A_pooled_prior_rf         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=kleineba@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED

module load palma/2020a
module load foss
module load R
...
# START THE APPLICATION
Rscript /home/k/kleineba/pipeline_hpc.R /scratch/tmp/kleineba/PreRarefaction/Data/A/Rarefied/rf11150 11150 TRUE _1 _2 32 FALSE /home/k/kleineba/pooled_priors_A.RDS
...