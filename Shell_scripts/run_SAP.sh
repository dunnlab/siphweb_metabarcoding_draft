#!/bin/bash

#SBATCH --mem-per-cpu=4g
#SBATCH -t 10:10:00
#SBATCH --array=0-100
#SBATCH --job-name=SAP_jobfile.txt
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --ntasks=1

/ysm-gpfs/apps/software/dSQ/0.92/dSQBatch.py SAP_jobfile.txt
