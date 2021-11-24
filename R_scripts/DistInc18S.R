library(tidyverse)
library(magrittr)
library(vegan)
library(seqinr)
library(ape)
library(phangorn)

#Load sequence alignments
setwd('/Volumes/GoogleDrive/My Drive/Primer design/BLASTpulls/18Sdistsandincompatibilities')
ALL18Saln <- read.fasta("../18Sblasts/18Souts.aligned.fasta")
ids18S = seqinr::getName(ALL18Saln)
#ALL16Saln <- read.fasta("../16Sblasts/16Sout.aligned.fasta")
#ids16S = seqinr::getName(ALL16Saln)


#### 18S ####

#Incompatible taxa for each primer pair
incV9 = read.csv('18SV9ADS_incompatible', sep = "|", header = F) %>% .[. != "AAAA"]
inc153 = read.csv('18S153_incompatible', sep = "|", header = F) %>% .[. != "AAAA"]
inc166 = read.csv('18S166_incompatible', sep = "|", header = F) %>% .[. != "AAAA"]
inc179 = read.csv('18S179_incompatible', sep = "|", header = F) %>% .[. != "AAAA"]
inc261 = read.csv('18S261_incompatible', sep = "|", header = F) %>% .[. != "AAAA"]
inc272 = read.csv('18S272_incompatible', sep = "|", header = F) %>% .[. != "AAAA"]
inc_all = list(incV9, inc153, inc166, inc179, inc261, inc272)

#Compatible taxa for each primer pair
comp_all = list()
for(i in 1:length(inc_all)){
  ids18S[which(!(ids18S %in% inc_all[[i]]))] -> comp_all[[i]]
}

#Taxonomic coverage
unlist(comp_all[c(2,6)]) %>% unique()
unlist(comp_all) %>% table() %>% hist()

#Species resolution %
dists = read.csv('18S261_dists.csv')
rownames(dists) = dists$X
dists %<>% .[,-1]
dists[which((rownames(dists) %in% inc)),which((colnames(dists) %in% inc261))] = NA
dists_total = dists[!is.na(dists)] %>% length()
length(dists[(dists < 100) & !is.na(dists)])/dists_total

#Included taxa eyeballing
ALL18Saln %>% getAnnot() %>% .[which(ids18S %in% inc_all[[6]])] %>% str_sub(start=13, end=38) -> excluded_taxa
ALL18Saln %>% getAnnot() %>% .[which(!(ids18S %in% inc_all[[1]]))] %>% str_sub(start=13, end=38) -> included_taxa
ALL18Saln %>% getAnnot() %>% str_sub(start=13, end=38) -> total_taxa

### Analyze inserts ###

setwd("/Volumes/GoogleDrive/My Drive/Primer design/Insert_extractions")

ins_134 <- read.fasta("18S134insert.fasta") %>% .[which(!(getName(.) %in% incV9))] #134 is V9ADS
ins_153 <- read.fasta("18S153insert.fasta") %>% .[which(!(getName(.) %in% inc153))]
ins_166 <- read.fasta("18S166insert.fasta") %>% .[which(!(getName(.) %in% inc166))]
ins_179 <- read.fasta("18S179insert.fasta") %>% .[which(!(getName(.) %in% inc179))]
ins_261 <- read.fasta("18S261insert.fasta") %>% .[which(!(getName(.) %in% inc261))]
ins_272 <- read.fasta("18S272insert.fasta") %>% .[which(!(getName(.) %in% inc272))]

#export pruned fastas
write.fasta(ins_134, names = getName(ins_134), file.out="18S134insert_pruned.fasta")
write.fasta(ins_153, names = getName(ins_153), file.out="18S134insert_pruned.fasta")
write.fasta(ins_166, names = getName(ins_166), file.out="18S166insert_pruned.fasta")
write.fasta(ins_179, names = getName(ins_179), file.out="18S179insert_pruned.fasta")
write.fasta(ins_261, names = getName(ins_261), file.out="18S261insert_pruned.fasta")
write.fasta(ins_272, names = getName(ins_272), file.out="18S272insert_pruned.fasta")


