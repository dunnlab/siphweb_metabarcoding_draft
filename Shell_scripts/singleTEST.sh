#!/bin/bash
#SBATCH --job-name=SingleTEST  
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=7-00:00:00

source activate qiime2-2018.11

#Declare temporary variables
runid="singlend"
pipeline_dir="singlendtrial"
pathtoFASTQs=fullpipelinetrial"/fastq_files"
barcodes=("152" "166" "272" "179" "261" "134")
queries=("^TGACGGAAGGGCACCACCAG|^TCCACCAACTAAGAACGGCC|^CTGGTGGTGCCCTTCCGTCA|^GGCCGTTCTTAGTTGGTGGA" "^AACGGCTACCACATCCAAGG|^CACCAGACTTGCCCTCCAAT|^CCTTGGATGTGGTAGCCGTT|^ATTGGAGGGCAAGTCTGGTG" "^AAACGATGCCGACTAGCGAT|^TCCACCAACTAAGAACGGCC|^ATCGCTAGTCGGCATCGTTT|^GGCCGTTCTTAGTTGGTGGA" "^GGCCGTTCTTAGTTGGTGGA|^TGCGGCCCAGAACATCTAAG|^TCCACCAACTAAGAACGGCC|^CTTAGATGTTCTGGGCCGCA" "^AACAGGTCTGTGATGCCCTT|^TGTGTACAAAGGGCAGGGAC|^AAGGGCATCACAGACCTGTT|^GTCCCTGCCCTTTGTACACA" "^CTTTGTACACACCGCCCGTC|^CCTTGTTACGACTTTTACTTCCTCT|^GACGGGCGGTGTGTACAAAG|^AGAGGAAGTAAAAGTCGTAACAAGG")
grep_output=$pipeline_dir/"grep_output"
cutadapt_o=$pipeline_dir"/cutadapt_o"

echo "Directories and paths declared!"

for barcode in ${barcodes[@]}; do

  echo Starting QIIME2 for barcode_"$barcode"

  manifest="$pipeline_dir"/manifest_se_"$barcode"
  metadata="$pipeline_dir"/metadata_"$barcode".tsv
  exportpath="$pipeline_dir"/Q2_"$barcode"/export

  #Load data into QIIME2
  
  qiime tools import \
     --type 'SampleData[SequencesWithQuality]' \
     --input-path $manifest \
     --output-path "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
     --input-format SingleEndFastqManifestPhred33

  qiime demux summarize \
    --i-data "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qzv

  echo "Paired reads imported into QIIME2!"

  #DADA2 Denoise and dereplicate
  qiime dada2 denoise-single \
          --i-demultiplexed-seqs "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
          --p-trunc-len 90 \
          --o-representative-sequences "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
          --o-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
          --p-n-threads 24 \
          --output-dir "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"_dada2_out

  qiime feature-table summarize \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qzv \
    --m-sample-metadata-file $metadata

  qiime feature-table tabulate-seqs \
    --i-data "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qzv

  qiime vsearch cluster-features-de-novo \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
    --i-sequences "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
    --p-perc-identity 0.95 \
    --o-clustered-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-095-table.qza \
    --o-clustered-sequences "$pipeline_dir"/Q2_"$barcode"/rep-095-"$runid""$barcode".qza

  qiime tools export \
    --input-path "$pipeline_dir"/Q2_"$barcode"/rep-095-"$runid""$barcode".qza \
    --output-path $exportpath

  qiime tools export \
  --input-path "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-095-table.qza \
  --output-path $exportpath

  echo "Data denoised, merged, and demultiplexed with DADA2 and clustered with VSEARCH!"

done
echo "QIIME2 loop completed!"

echo "End of the pipeline!"
