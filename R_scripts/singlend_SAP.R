library(phyloseq)
library(tidyverse)
library(reshape2)
library(scales)
library(RColorBrewer)
library(forcats)
library(vegan)

runid = "test_3seqs"
barcodes=c("134", "152","166", "179", "261", "272")
Assignments152 <- read.csv(file = "Desktop/test_152_pooled_assignments.csv", header = FALSE, sep = ",")[,c(4:10)]
Assignments166 <- read.csv(file = "Desktop/test_166_pooled_assignments.csv", header = FALSE, sep = ",")[,c(4:10)]
Features152 <- read.csv(file = "Desktop/test_3seqs152-feature-table.tsv", header = FALSE, sep = "\t")[-1,]
Features166 <- read.csv(file = "Desktop/test_3seqs166-feature-table.tsv", header = FALSE, sep = "\t")[-1,]

#Order the OTUs as in the features table
Features152_pruned = Features152[c(1,which(Features152$V1 %in% Assignments152$V4)),]
Taxonomy152 <- Assignments152[match(Features152_pruned[2:nrow(Features152_pruned),1], Assignments152[,1]),2:7]
Features166_pruned = Features166[c(1,which(Features166$V1 %in% Assignments166$V4)),]
Taxonomy166 <- Assignments166[match(Features166_pruned[2:nrow(Features166_pruned),1], Assignments166[,1]),2:7]


#Paste the ranks into a single string
Ranks_pasted152 = apply(Taxonomy152, 1, function(x){paste(x, collapse="_")})
Ranks_pasted166 = apply(Taxonomy166, 1, function(x){paste(x, collapse="_")})

#Combine the taxonomic ranks ID with the feature statistics
IDtable152 <- cbind(c("OTU",Ranks_pasted152), Features152_pruned[,2:ncol(Features152_pruned)])
IDtable166 <- cbind(c("OTU",Ranks_pasted166), Features166_pruned[,2:ncol(Features166_pruned)])

#Write an outfile
IDtable152<-apply(IDtable152, 2, as.character)
colnames(IDtable152) = as.vector(IDtable152[1,])
IDtable152 = as.data.frame(IDtable152[-1,])
write.table(x=IDtable152, file="Desktop/152_IDtable.tsv", row.names=FALSE)

IDtable166<-apply(IDtable166, 2, as.character)
colnames(IDtable166) = as.vector(IDtable166[1,])
IDtable166 = as.data.frame(IDtable166[-1,])
write.table(x=IDtable166, file="Desktop/166_IDtable.tsv", row.names=FALSE)


summary_plot <- function(x, name){
  x  %>% group_by(variable) %>% summarise(N=n()) %>% ggplot(aes(x=fct_reorder(variable, N), y=N)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


IDTables = list()
i=1
for(j in list.files("Desktop/IDtables/")){
  tableJ = read.table(paste("Desktop/IDtables/", j, sep=""), header=T)
  names(tableJ) = str_replace_all(names(tableJ), "\\w+_(\\w+_\\w+)_.+$", "\\1")
  tableJ <- melt(tableJ, id.vars = c('OTU'))
  for(sample in unique(tableJ$variable)){
    aliquot = tableJ[which(tableJ$variable==sample),]
    #aliquot = aliquot[which(aliquot$value<max(aliquot$value) & aliquot$value>100),]
    #tableJ <- tableJ[-which(tableJ$variable==sample),]
    #tableJ <- rbind(tableJ, aliquot)
    #print(aliquot[which(aliquot$value>100),])
  }
  tableJ = tableJ[which(tableJ$value > 1),]
  print(quantile(tableJ$value)[2])
  tableJ = tableJ[which(tableJ$value > 30),]
  #tableJ = tableJ[which(tableJ$value > 1),]
  IDTables[[i]] = tableJ
  i=i+1
}
names(IDTables) = barcodes


i = 1
plotlist = list()
for(t in IDTables){
  plot_T <- summary_plot(t) + ggtitle(names(IDTables)[i])
  plotlist[[i]] <- plot_T
  i=i+1
}
multiplot(plotlist=plotlist, cols=3)

i = 1
plotlist2 = list()
for(t in IDTables){
  plotlist2[[i]] <- ggplot(t, aes(x = variable, y = log(value), fill = OTU)) + geom_bar(position = "fill",stat = "identity") + theme(legend.key.size = unit(0.1,"line"), legend.text = element_text(size=4),legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=rep(c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"),10)) + ggtitle(names(IDTables)[i]) + guides(shape = guide_legend(override.aes = list(size = 1)))
  plotlist2[[i]] %>% print()
  i=i+1
}
multiplot(plotlist=plotlist2[3:6], cols=2)