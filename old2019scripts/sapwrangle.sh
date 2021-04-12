#!/bin/bash
#SBATCH --job-name=SAP-to-Taxa 
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-00:10:00

module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1
runid="singlend"
pipeline_dir="singlendtrial"
barcodes=("134" "152" "166" "179" "261" "272")

echo "Starting SAP loop!"

for barcode in ${barcodes[@]}; do
  
  cat "$pipeline_dir"/Q2_"$barcode"/splitSAP"$runid""$barcode"/*/assignments.csv > "$pipeline_dir"/Q2_"$barcode"/"$runid"_"$barcode"_pooled_assignments.csv
  tail -n +2 "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-095-table.tsv | sed  's/^#OTU //'  > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-095-table-1.tsv

  pathtoAssignments="$pipeline_dir"/Q2_"$barcode"/"$runid"_"$barcode"_pooled_assignments.csv
  pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-095-table-1.tsv

  #Create table with OTU IDs and frequency in each sample
  #								  Arg1     Arg2          Arg3             Arg4                     Arg5
  Rscript --vanilla OTUwrangle.R $runid $barcode $pathtoAssignments $pathtoFeatures "$pipeline_dir"/Q2_"$barcode"/

done
