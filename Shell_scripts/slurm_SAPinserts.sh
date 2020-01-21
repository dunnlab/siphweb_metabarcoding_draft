#!/bin/bash
#SBATCH --partition=scavenge
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem-per-cpu=5G
#SBATCH -t 2-24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alejandro.damianserrano@yale.edu
#SBATCH --job-name=sap_insilicos

module load BLAST+/2.7.1-foss-2016b-Python-2.7.13
module load ClustalW2/2.1-foss-2016b

~/anaconda2/bin/sap --project 18S134 --email alejandro.damianserrano@yale.edu dereplicated_inserts/dereplicated.134.raw.fasta
~/anaconda2/bin/sap --project 18S166 --email alejandro.damianserrano@yale.edu dereplicated_inserts/dereplicated.166.raw.fasta
~/anaconda2/bin/sap --project 18S261 --email alejandro.damianserrano@yale.edu dereplicated_inserts/dereplicated.261.raw.fasta
~/anaconda2/bin/sap --project 18S153 --email alejandro.damianserrano@yale.edu dereplicated_inserts/dereplicated.153.raw.fasta  
~/anaconda2/bin/sap --project 18S179 --email alejandro.damianserrano@yale.edu dereplicated_inserts/dereplicated.179.raw.fasta
~/anaconda2/bin/sap --project 18S272 --email alejandro.damianserrano@yale.edu dereplicated_inserts/dereplicated.272.raw.fasta