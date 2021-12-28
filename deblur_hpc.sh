#!/bin/bash

#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=36        # the number of CPU cores per node
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=2-00:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --account=uni
#SBATCH --mem=64000

#SBATCH --job-name=Data_B_deblur_rf       # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=kleineba@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED
conda activate qiime2-2021-11

qiime deblur denoise-16S \
  --i-demultiplexed-seqs /scratch/tmp/kleineba/PreRarefaction/Data/B/Rarefied/rf21092/filtered/Manifest_B_rarefied.qza \
  --o-table table-deblur_B_rf.qza \
  --p-trim-length 229 \
  --p-sample-stats \
  --o-stats deblur-stats_B_rf.qza \
  --p-min-size 2 \
  --p-jobs-to-start 32 \
  --o-representative-sequences rep-seqs-deblur_B_rf.qza