#!/bin/bash
#SBATCH --job-name=miseq_fastqc

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-06:00:00

fastqc *.fastq
