#!/bin/bash
#SBATCH --job-name=deadSAP  
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-999
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=0-01:00:00

module load dSQ
dSQ --jobfile SAP_jobfile.txt --mem-per-cpu=4g -t 12:12:00 > runSAP.sh
