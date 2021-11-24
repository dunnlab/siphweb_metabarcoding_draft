#!/bin/bash
#SBATCH --job-name=dada2

#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=0-13:00:00

Rscript rscriptdada.R
