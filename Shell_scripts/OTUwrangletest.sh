#!/bin/bash
#SBATCH --job-name=Q2_pipe  
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=0-01:00:00

module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

runid="test_3seqs"
barcode="ALL"
pipeline_dir="fullpipelinetrial"

Rscript --vanilla "$pipeline_dir"/OTUwrangle.R $runid $barcode noprimers_allseqs/assignments.csv featuretable3lines.tsv "."
