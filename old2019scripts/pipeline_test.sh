#!/bin/bash
#SBATCH --job-name=Q2_pipe  
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=2-23:00:00

### USAGE ###

# Make sure the manifest file and metadata file are formatted adequately
# mkdir pipeline_directory
# Have a copy of the manifest and the metadata in the pipeline_directory
# Manifest paths should be $PWD/$bb_output/filename for each filename


#############

#Load packages
source activate qiime2-2018.11
module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1
module load FastQC/0.11.5-Java-1.8.0_121
module load BBMap/36.62-foss-2016b-Java-1.8.0_121

echo "Modules and environment loaded!"

#Declare temporary variables
metadata="pipelinetrial/metadata_byamplicon.tsv"
manifest="pipelinetrial/manifest_pipe"
runid="TEST"
pathtoTARGZs="/SAY/archive/cwd7-11174224-FASEEB-A/data/sequences/illumina/SIPHWEB_MiSeq_Preliminary_001"
pathtoFASTQs="pipelinetrial/fastq_files"
pathtoHTMLS="pipelinetrial/fastqc_HTMLs"
barcodes=("152" "166" "272" "179" "261" "134")
queries=("^TGACGGAAGGGCACCACCAG|^TCCACCAACTAAGAACGGCC|^CTGGTGGTGCCCTTCCGTCA|^GGCCGTTCTTAGTTGGTGGA" "^AACGGCTACCACATCCAAGG|^CACCAGACTTGCCCTCCAAT|^CCTTGGATGTGGTAGCCGTT|^ATTGGAGGGCAAGTCTGGTG" "^AAACGATGCCGACTAGCGAT|^TCCACCAACTAAGAACGGCC|^ATCGCTAGTCGGCATCGTTT|^GGCCGTTCTTAGTTGGTGGA" "^GGCCGTTCTTAGTTGGTGGA|^TGCGGCCCAGAACATCTAAG|^TCCACCAACTAAGAACGGCC|^CTTAGATGTTCTGGGCCGCA" "^AACAGGTCTGTGATGCCCTT|^TGTGTACAAAGGGCAGGGAC|^AAGGGCATCACAGACCTGTT|^GTCCCTGCCCTTTGTACACA" "^CTTTGTACACACCGCCCGTC|^CCTTGTTACGACTTTTACTTCCTCT|^GACGGGCGGTGTGTACAAAG|^AGAGGAAGTAAAAGTCGTAACAAGG")
grep_output="pipelinetrial/grep_output"
bb_output="pipelinetrial/bb_output"
bb_orphans="pipelinetrial/bb_orphans"
cutadapt_o="pipelinetrial/cutadapt_o"
exportpath="pipelinetrial/dada2_exports"

mkdir $pathtoFASTQs
mkdir $pathtoHTMLS
mkdir $grep_output
mkdir $bb_output
mkdir $bb_orphans
mkdir $cutadapt_o
mkdir $exportpath

echo "Directories and paths declared!"

#Extract GZ balls
gzfiles=($(ls "$pathtoTARGZs"))
for f in ${gzfiles[@]}; do gunzip -c "$pathtoTARGZs"/"$f" > "$pathtoFASTQs"/"${f%.*}" ; done

echo ".gz files extracted!"

#Quality HTMLS
fastqc "$pathtoFASTQs"/* -o $pathtoHTMLS

echo "Fastqc quality html reports generated!"

#Segragate by Barcode/Amplicon
filenames=($(ls "$pathtoFASTQs")) #for the grep by amplicon stage
for file in ${filenames[@]}; do
	for bc in 0 1 2 3 4 5; do
		grep -E "${queries[$bc]}" -C 1 -A 2 pairedendsequences/"$file"  | sed 's\^--$\\g' | sed '/^$/d' >> "$grep_output"/"${barcodes[$bc]}"_"$file"
	done
done

echo "Reads sorted by PCR primers into the 6 different amplicons!"

#Remove PCR primers
files=($(ls "$grep_output")) #for the cutadapt PCR primer removal
for file in ${files[@]}; do
	cutadapt -m 50 -g ^TGACGGAAGGGCACCACCAG -g ^TCCACCAACTAAGAACGGCC -g ^CTGGTGGTGCCCTTCCGTCA -g ^GGCCGTTCTTAGTTGGTGGA  -g ^AACGGCTACCACATCCAAGG -g ^CACCAGACTTGCCCTCCAAT -g ^CCTTGGATGTGGTAGCCGTT -g ^ATTGGAGGGCAAGTCTGGTG  -g ^AAACGATGCCGACTAGCGAT -g ^TCCACCAACTAAGAACGGCC -g ^ATCGCTAGTCGGCATCGTTT -g ^GGCCGTTCTTAGTTGGTGGA  -g ^GGCCGTTCTTAGTTGGTGGA -g ^TGCGGCCCAGAACATCTAAG -g ^TCCACCAACTAAGAACGGCC -g ^CTTAGATGTTCTGGGCCGCA  -g ^AACAGGTCTGTGATGCCCTT -g ^TGTGTACAAAGGGCAGGGAC -g ^AAGGGCATCACAGACCTGTT -g ^GTCCCTGCCCTTTGTACACA  -g ^CTTTGTACACACCGCCCGTC -g ^CCTTGTTACGACTTTTACTTCCTCT -g ^GACGGGCGGTGTGTACAAAG -g ^AGAGGAAGTAAAAGTCGTAACAAGG -o "$cutadapt_o"/"$file" "$grep_output"/"$file"
done

echo "PCR primer sequences removed from the 5' end!"

#Repair non-matching read pairs due to differential sequencing errors in R1 and R2
samplenames=($(grep -oE "^\w+-\w+" $manifest))
for file in ${samplenames[@]}; do
	repair.sh in="$cutadapt_o"/"$file"_R1_001.fastq in2="$cutadapt_o"/"$file"_R2_001.fastq out1="$bb_output"/"$file"_R1_001.fastq out2="$bb_output"/"$file"_R2_001.fastq outsingle="$bb_orphans"/"$file"_orphans.fastq
done

echo "R1 and R2 reads paired up and ready to import!"

# *** IMPORTANT *** MANIFEST MUST CONTAIN THE PATHS TO BB_OUTPUT!

#Load data into QIIME2
qiime metadata tabulate --m-input-file $metadata --o-visualization "$runid"_metadata.qzv

qiime tools import \
   --type 'SampleData[PairedEndSequencesWithQuality]' \
   --input-path $manifest \
   --output-path "$runid".qza \
   --input-format PairedEndFastqManifestPhred33

qiime demux summarize \
  --i-data "$runid".qza \
  --o-visualization "$runid".qzv

echo "Paired reads imported into QIIME2!"

#DADA2 Denoise and dereplicate
qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$runid".qza \
        --p-trim-left-f 20 \
        --p-trim-left-r 20 \
        --p-trunc-len-f 150 \
        --p-trunc-len-r 150 \
        --o-representative-sequences rep-"$runid".qza \
        --o-table "$runid"-table.qza \
        --p-n-threads 24 \
        --output-dir "$runid"_dada2_out

qiime feature-table summarize \
  --i-table "$runid"-table.qza \
  --o-visualization "$runid"-table.qzv \
  --m-sample-metadata-file metadata_"$runid".tsv

qiime feature-table tabulate-seqs \
  --i-data rep-"$runid".qza \
  --o-visualization rep-"$runid".qzv

qiime tools export \
  --input-path rep-"$runid".qza \
  --output-path $exportpath

echo "Data denoised, merged, and demultiplexed with DADA2!"

#Multiple seqeunce alignment using Mafft
 qiime alignment mafft \
  --i-sequences rep-"$runid".qza \
  --o-alignment aligned-rep-"$runid".qza

#Mask the alignment to remove positions that are highly variable.
qiime alignment mask \
  --i-alignment aligned-rep-"$runid".qza \
  --o-masked-alignment masked-aligned-rep-"$runid".qza

echo "Sequences aligned!"

#Create the tree using the Fasttree program
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-"$runid".qza \
  --o-tree unrooted-"$runid".qza

#Root the tree using the longest root
qiime phylogeny midpoint-root \
  --i-tree unrooted-"$runid".qza \
  --o-rooted-tree rooted-"$runid".qza

echo "Phylogenetic tree generated!"

#Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table "$runid"-table.qza \
  --i-phylogeny rooted-"$runid".qza \
  --p-max-depth 4000 \
  --m-metadata-file $metadata \
  --o-visualization "$runid"-rarefaction.qzv

#Diversity metrics
qiime diversity "$runid"-metrics-phylogenetic \
  --i-phylogeny rooted-"$runid".qza \
  --i-table "$runid"-table.qza \
  --p-sampling-depth 1109 \
  --m-metadata-file $metadata \
  --output-dir "$runid"-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$runid"-metrics-byamaplicon/faith_pd_vector.qza \
  --m-metadata-file $metadata.tsv \
  --o-visualization "$runid"-metrics-byamplicon/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$runid"-metrics-byamplicon/evenness_vector.qza \
  --m-metadata-file $metadata.tsv \
  --o-visualization "$runid"-metrics-byamplicon/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity "$runid"-metrics-byamplicon/shannon_vector.qza \
  --m-metadata-file $metadata.tsv \
  --o-visualization "$runid"-metrics-byamplicon/shannon_group-significance.qzv

#Beta diversity by template
qiime diversity beta-group-significance \
  --i-distance-matrix "$runid"-metrics-byamplicon/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $metadata \
  --m-metadata-column template \
  --o-visualization "$runid"-metrics-byamplicon/unweighted-unifrac-template-significance.qzv \
  --p-pairwise

echo "Diversity metrics computed and plotted. End of the pipeline!"
