#!/bin/bash
#SBATCH --job-name=load_qiime2_qza
#SBATCH --partition=general
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1-23:00:00

qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path manifest \
   --output-path paired-end-demux.qza \
   --input-format PairedEndFastqManifestPhred33