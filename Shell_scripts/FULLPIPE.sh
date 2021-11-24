#!/bin/bash
#SBATCH --job-name=part2pipe
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=1-00:00:00

### USAGE ###

# Make sure the manifest files and metadata files are formatted adequately
# mkdir pipeline_directory
# Have a copy of the manifest and metadata files in the pipeline_directory
# Manifests paths should be $PWD/$bb_output/filename for each filename
# OTUwrangle.R script should be in the pipeline_directory

#############

#Declare temporary variables
runid="test_3seqs"
pipeline_dir="fullpipelinetrial"
barcodes=("152" "166" "272" "179" "261" "134")

rm -rf $pipelinedir/Q2_*

# *** IMPORTANT *** MANIFEST MUST CONTAIN THE PATHS TO BB_OUTPUT!
source activate qiime2-2018.11
module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

for barcode in ${barcodes[@]}; do

  echo Starting QIIME2 for barcode_"$barcode"

  manifest="$pipeline_dir"/manifest_"$barcode"
  metadata="$pipeline_dir"/metadata_"$barcode".tsv
  exportpath="$pipeline_dir"/Q2_"$barcode"/export

  mkdir "$pipeline_dir"/Q2_"$barcode"
  mkdir $exportpath

  #Load data into QIIME2
  qiime metadata tabulate --m-input-file $metadata --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"_metadata.qzv

  qiime tools import \
     --type 'SampleData[PairedEndSequencesWithQuality]' \
     --input-path $manifest \
     --output-path "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
     --input-format PairedEndFastqManifestPhred33

  qiime demux summarize \
    --i-data "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qzv

  echo "Paired reads imported into QIIME2!"

  #DADA2 Denoise and dereplicate
  qiime dada2 denoise-paired \
          --i-demultiplexed-seqs "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode".qza \
          --p-trunc-len-f 125 \
          --p-trunc-len-r 125 \
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

  qiime tools export \
    --input-path "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
    --output-path $exportpath

  qiime tools export \
  --input-path "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
  --output-path $exportpath

  echo "Data denoised, merged, and demultiplexed with DADA2!"

  #Multiple seqeunce alignment using Mafft
   qiime alignment mafft \
    --i-sequences "$pipeline_dir"/Q2_"$barcode"/rep-"$runid""$barcode".qza \
    --o-alignment "$pipeline_dir"/Q2_"$barcode"/aligned-rep-"$runid""$barcode".qza

  #Mask the alignment to remove positions that are highly variable.
  qiime alignment mask \
    --i-alignment "$pipeline_dir"/Q2_"$barcode"/aligned-rep-"$runid""$barcode".qza \
    --o-masked-alignment "$pipeline_dir"/Q2_"$barcode"/masked-aligned-rep-"$runid""$barcode".qza

  echo "Sequences aligned!"

  #Create the tree using the Fasttree program
  qiime phylogeny fasttree \
    --i-alignment "$pipeline_dir"/Q2_"$barcode"/masked-aligned-rep-"$runid""$barcode".qza \
    --o-tree "$pipeline_dir"/Q2_"$barcode"/unrooted-"$runid""$barcode".qza

  #Root the tree using the longest root
  qiime phylogeny midpoint-root \
    --i-tree "$pipeline_dir"/Q2_"$barcode"/unrooted-"$runid""$barcode".qza \
    --o-rooted-tree "$pipeline_dir"/Q2_"$barcode"/rooted-"$runid""$barcode".qza

  echo "Phylogenetic tree generated!"

  #Alpha rarefaction
  qiime diversity alpha-rarefaction \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
    --i-phylogeny "$pipeline_dir"/Q2_"$barcode"/rooted-"$runid""$barcode".qza \
    --p-max-depth 4000 \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-rarefaction.qzv

  #Alpha diversity metrics
  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$pipeline_dir"/Q2_"$barcode"/rooted-"$runid""$barcode".qza \
    --i-table "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-table.qza \
    --p-sampling-depth 1100 \
    --m-metadata-file $metadata \
    --output-dir "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/faith_pd_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"faith-pd-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/evenness_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"evenness-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/shannon_vector.qza \
    --m-metadata-file $metadata \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"shannon_group-significance.qzv

  #Beta diversity by template
  qiime diversity beta-group-significance \
    --i-distance-matrix "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file $metadata \
    --m-metadata-column template \
    --o-visualization "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-metrics-results/"$runid""$barcode"unweighted-unifrac-template-significance.qzv \
    --p-pairwise

  echo Diversity metrics computed and plotted. End of the QIIME process for barcode "$barcode"

  #Convert feature table to TSV
  biom convert --to-tsv -i "$exportpath"/feature-table.biom -o "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-table.tsv

  #Rename and relocate sequences for all samples
  #cat "$exportpath"/dna-sequences.fasta > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-sequences.fasta
  head -6 "$exportpath"/dna-sequences.fasta > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-sequences.fasta

done

echo "QIIME2 loop completed!"

#### ASSIGN TAXONOMY TO AMPLICONS ####

source deactivate qiime2-2018.11
module purge
module load BLAST+/2.7.1-foss-2016b-Python-2.7.13
module load ClustalW2/2.1-foss-2016b
module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1

"Starting SAP loop!"

for barcode in ${barcodes[@]}; do
  echo SAP started for barcode $barcode !
  ~/anaconda2/bin/sap --project "$pipeline_dir"/Q2_"$barcode"/SAP_"$runid"_"$barcode" \
    --email alejandro.damianserrano@yale.edu \
    --minsignificance 0.3 \
    -s 0.1 \
    --minidentity 0.5 \
    -n 10 \
    --ppcutoff 85 \
    --svg "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-sequences.fasta

  pathtoAssignments="$runid"_"$barcode"/assignments.csv
  pathtoFeatures="$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-feature-table.tsv

  #Create table with OTU IDs and frequency in each sample
  Rscript --vanilla "$pipeline_dir"/OTUwrangle.R $runid $barcode $pathtoAssignments $pathtoFeatures "$pipeline_dir"/Q2_"$barcode"
  tail -n +2 "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""IDtable.tsv" > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""_IDtable.tsv"
  rm "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode""IDtable.tsv"

  echo $barcode Sequences identified by SAP!

done

echo "End of the pipeline!"

