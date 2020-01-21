#!/bin/bash
#SBATCH --job-name=assignmentR
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=1-00:00:00

runid="RUN1"
pipeline_dir="fullpipelineRUN1"
pathtoTARGZs="/SAY/standard2/cwd7-CC0522-FASEEB/data/sequences/illumina/SIPHWEB_MiSeq_Lane_1"
pathtoFASTQs=$pipeline_dir"/fastq_files"
pathtoHTMLS=$pipeline_dir"/fastqc_HTMLs"
barcodes=("152" "166" "272" "179" "261" "134")
queries=("^TGACGGAAGGGCACCACCAG|^TCCACCAACTAAGAACGGCC|^CTGGTGGTGCCCTTCCGTCA|^GGCCGTTCTTAGTTGGTGGA" "^AACGGCTACCACATCCAAGG|^CACCAGACTTGCCCTCCAAT|^CCTTGGATGTGGTAGCCGTT|^ATTGGAGGGCAAGTCTGGTG" "^AAACGATGCCGACTAGCGAT|^TCCACCAACTAAGAACGGCC|^ATCGCTAGTCGGCATCGTTT|^GGCCGTTCTTAGTTGGTGGA" "^GGCCGTTCTTAGTTGGTGGA|^TGCGGCCCAGAACATCTAAG|^TCCACCAACTAAGAACGGCC|^CTTAGATGTTCTGGGCCGCA" "^AACAGGTCTGTGATGCCCTT|^TGTGTACAAAGGGCAGGGAC|^AAGGGCATCACAGACCTGTT|^GTCCCTGCCCTTTGTACACA" "^CTTTGTACACACCGCCCGTC|^CCTTGTTACGACTTTTACTTCCTCT|^GACGGGCGGTGTGTACAAAG|^AGAGGAAGTAAAAGTCGTAACAAGG")
grep_output=$pipeline_dir/"grep_output"
bb_output=$pipeline_dir"/bb_output"
bb_orphans=$pipeline_dir"/bb_orphans"
cutadapt_o=$pipeline_dir"/cutadapt_o"
manifest="$pipeline_dir"/manifest
metadata="$pipeline_dir"/metadata.tsv

module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

for barcode in ${barcodes[@]}; do

  cat "$pipeline_dir"/Q2_"$barcode"/splitSAP"$runid""$barcode"/*/assignments.csv > "$pipeline_dir"/Q2_"$barcode"/"$runid"_"$barcode"_pooled_assignments.csv

  pathtoAssignments="$pipeline_dir"/Q2_"$barcode"/"$runid"_"$barcode"_pooled_assignments.csv
  pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-table.tsv

  #Create table with OTU IDs and frequency in each sample
  Rscript --vanilla "$pipeline_dir"/OTUwrangle.R $runid $barcode $pathtoAssignments $pathtoFeatures "$pipeline_dir"/Q2_"$barcode"
  tail -n +2 "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""IDtable.tsv" > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""_IDtable.tsv"
  rm "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""IDtable.tsv"

done
