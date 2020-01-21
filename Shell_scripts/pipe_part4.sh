#!/bin/bash
#SBATCH --job-name=pipe_part4
#SBATCH --partition=general
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --time=1-00:00:00

runid="test"
pipeline_dir="fullpipelinetrial"
barcodes=("152" "166" "272" "179" "261" "134")
module load dSQ

for barcode in ${barcodes[@]}; do
  exportpath="$pipeline_dir"/Q2_"$barcode"/export

  #Split sequences in sets of 100 for SAP
  split -l 100 "$exportpath"/dna-sequences.fasta "$pipeline_dir"/Q2_"$barcode"/splitout"$runid""$barcode"/splitout
  #head -6 "$exportpath"/dna-sequences.fasta > "$pipeline_dir"/Q2_"$barcode"/"$runid""$barcode"-3sequences.fasta
  splitfiles=($(ls "$pipeline_dir"/Q2_"$barcode"/splitout"$runid""$barcode"))

  #Make DeadSimple jobfile for each fasta subfile
  for i in ${splitfiles[@]}; do
    echo "module load BLAST+/2.7.1-foss-2016b-Python-2.7.13; module load ClustalW2/2.1-foss-2016b; ~/anaconda2/bin/sap --project $pipeline_dir/Q2_$barcode/splitSAP$runid$barcode/$i --email alejandro.damianserrano@yale.edu --minsignificance 0.3 -s 0.1 --minidentity 0.5 -n 10 --ppcutoff 85 --svg $pipeline_dir/Q2_$barcode/splitout$runid$barcode/$i" | cat >> "$pipeline_dir"/Q2_"$barcode"/SAP"$runid""$barcode"_jobfile.txt
  done

  dSQ --jobfile "$pipeline_dir"/Q2_"$barcode"/SAP"$runid""$barcode"_jobfile.txt --mem-per-cpu=4g -t 7-10:00:00 > "$pipeline_dir"/Q2_"$barcode"/runSAP_"$runid""$barcode".sh
  sbatch "$pipeline_dir"/Q2_"$barcode"/runSAP_"$runid""$barcode".sh
  echo SAP started for barcode "$barcode"!

done


