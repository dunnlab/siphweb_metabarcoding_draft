#!/bin/bash
#SBATCH --job-name=SAP-to-Taxa 
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=0-05:00:00

module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1
runid="test_3seqs"
pipeline_dir="fullpipelinetrial"
barcodes=("152" "166" "272" "179" "261" "134")

echo "Starting SAP loop!"

for barcode in ${barcodes[@]}; do
  
  cat "$pipeline_dir"/Q2_"$barcode"/splitSAP"$runid""$barcode"/*/assignments.csv > "$pipeline_dir"/Q2_"$barcode"/"$runid"_"$barcode"_pooled_assignments.csv

  pathtoAssignments="$pipeline_dir"/Q2_"$barcode"/"$runid"_"$barcode"_pooled_assignments.csv
  pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-table.tsv

  #Create table with OTU IDs and frequency in each sample
  Rscript --vanilla "$pipeline_dir"/OTUwrangle.R $runid $barcode $pathtoAssignments $pathtoFeatures "$pipeline_dir"/Q2_"$barcode"
  tail -n +2 "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""IDtable.tsv" > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""_IDtable.tsv"
  rm "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""IDtable.tsv"

done
