#!/bin/bash
#SBATCH --job-name=Q2_pipe  
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=2-23:00:00

#Load packages
source activate qiime2-2018.11
module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.1
module load FastQC/0.11.5-Java-1.8.0_121
module load BBMap/36.62-foss-2016b-Java-1.8.0_121

#Declare temporary variables
metadata="pathtometadata.tsv"
manifest="pathtomanifest"
runid="runid"
pathtoTARGZs="/SAY/archive/cwd7-11174224-FASEEB-A/data/sequences/illumina/path"
pathtoFASTQs="path/"
pathtoHTMLS="path/"
barcodes=("152" "166" "272" "179" "261" "134")
queries=("^TGACGGAAGGGCACCACCAG|^TCCACCAACTAAGAACGGCC|^CTGGTGGTGCCCTTCCGTCA|^GGCCGTTCTTAGTTGGTGGA" "^AACGGCTACCACATCCAAGG|^CACCAGACTTGCCCTCCAAT|^CCTTGGATGTGGTAGCCGTT|^ATTGGAGGGCAAGTCTGGTG" "^AAACGATGCCGACTAGCGAT|^TCCACCAACTAAGAACGGCC|^ATCGCTAGTCGGCATCGTTT|^GGCCGTTCTTAGTTGGTGGA" "^GGCCGTTCTTAGTTGGTGGA|^TGCGGCCCAGAACATCTAAG|^TCCACCAACTAAGAACGGCC|^CTTAGATGTTCTGGGCCGCA" "^AACAGGTCTGTGATGCCCTT|^TGTGTACAAAGGGCAGGGAC|^AAGGGCATCACAGACCTGTT|^GTCCCTGCCCTTTGTACACA" "^CTTTGTACACACCGCCCGTC|^CCTTGTTACGACTTTTACTTCCTCT|^GACGGGCGGTGTGTACAAAG|^AGAGGAAGTAAAAGTCGTAACAAGG")
grep_output="pathtogrepoutput"
bb_output="pathtoBBoutput"
bb_orphans="pathtoBBorphans"
cutadapt-o="pathtocutadaptoutput"
exportpath="path/"

#Extract GZ files
gzfiles=($(ls "$pathtoTARGZs"))
for f in ${gzfiles[@]}; do gunzip -c "$pathtoTARGZs"/"$f" > "$pathtoFASTQs"/"${f%.*}" ; done

#Quality HTMLS
fastqc "$pathtoFASTQs"/* -o $pathtoHTMLS

#Segragate by Barcode/Amplicon
filenames=($(ls "$pathtoFASTQs")) #for the grep by amplicon stage
for file in ${filenames[@]}; do
	for bc in 0 1 2 3 4 5; do
		grep -E "${queries[$bc]}" -C 1 -A 2 pairedendsequences/"$file"  | sed 's\^--$\\g' | sed '/^$/d' >> "$grep_output"/"${barcodes[$bc]}"_"$file"
	done
done

#Remove PCR primers
files=($(ls "$grep_output")) #for the cutadapt PCR primer removal
for file in ${files[@]}; do
	cutadapt -m 50 -g ^TGACGGAAGGGCACCACCAG -g ^TCCACCAACTAAGAACGGCC -g ^CTGGTGGTGCCCTTCCGTCA -g ^GGCCGTTCTTAGTTGGTGGA  -g ^AACGGCTACCACATCCAAGG -g ^CACCAGACTTGCCCTCCAAT -g ^CCTTGGATGTGGTAGCCGTT -g ^ATTGGAGGGCAAGTCTGGTG  -g ^AAACGATGCCGACTAGCGAT -g ^TCCACCAACTAAGAACGGCC -g ^ATCGCTAGTCGGCATCGTTT -g ^GGCCGTTCTTAGTTGGTGGA  -g ^GGCCGTTCTTAGTTGGTGGA -g ^TGCGGCCCAGAACATCTAAG -g ^TCCACCAACTAAGAACGGCC -g ^CTTAGATGTTCTGGGCCGCA  -g ^AACAGGTCTGTGATGCCCTT -g ^TGTGTACAAAGGGCAGGGAC -g ^AAGGGCATCACAGACCTGTT -g ^GTCCCTGCCCTTTGTACACA  -g ^CTTTGTACACACCGCCCGTC -g ^CCTTGTTACGACTTTTACTTCCTCT -g ^GACGGGCGGTGTGTACAAAG -g ^AGAGGAAGTAAAAGTCGTAACAAGG -o "$cutadapt-o"/"$file" "$grep_output"/"$file"
done

#Repair non-matching read pairs due to differential sequencing errors in R1 and R2
samplenames=($(grep -oE "^\w+-\w+" $manifest))
for file in ${samplenames[@]}; do
	repair.sh in="$cutadapt-o"/"$file"_R1_001.fastq in2="$cutadapt-o"/"$file"_R2_001.fastq out1="$bb_output"/"$file"_R1_001.fastq out2="$bb_output"/"$file"_R2_001.fastq outsingle="$bb_orphans"/"$file"_orphans.fastq
done

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

#Multiple seqeunce alignment using Mafft
 qiime alignment mafft \
  --i-sequences rep-"$runid".qza \
  --o-alignment aligned-rep-"$runid".qza

#Mask the alignment to remove positions that are highly variable.
qiime alignment mask \
  --i-alignment aligned-rep-"$runid".qza \
  --o-masked-alignment masked-aligned-rep-"$runid".qza

#Create the tree using the Fasttree program
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-"$runid".qza \
  --o-tree unrooted-"$runid".qza

#Root the tree using the longest root
qiime phylogeny midpoint-root \
  --i-tree unrooted-"$runid".qza \
  --o-rooted-tree rooted-"$runid".qza

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
