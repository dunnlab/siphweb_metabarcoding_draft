library(phyloseq)
library(tidyverse)
library(reshape2)
library(scales)
library(RColorBrewer)
library(forcats)
library(vegan)
#getwd()
#list.files("Desktop")
#import_biom("Desktop/feature-table.biom")
feature272 = read.table("Desktop/test272-feature-table.tsv", header=T)
names(feature272) = str_replace_all(names(feature272), "^X", "")
names(feature272) = str_replace_all(names(feature272), "\\w+_(\\w+_\\w+)_.+$", "\\1")
featureMelt272 <- melt(feature272, id.vars = c('OTU'))
featureCompact272 = featureMelt272[which(featureMelt272$value > 1),]
featureRecoded272 = featureCompact272
featureRecoded272$OTU = as.character(1:length(featureRecoded272$OTU))
ggplot(featureRecoded272, aes(x = variable, y = value, fill = OTU)) + geom_bar(position = "fill",stat = "identity") + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))
unique(featureCompact272$OTU)
for(i in 2:ncol(feature272)){
  feature272$`21_mix1`[which(feature272[,i]>0)] %>% length() %>% paste(feature272[which(feature272[,i]>0),1], ., sep="_") %>%  print()
}
feature272$`21_mix1`[which(feature272$`21_mix1`>0)] %>% length()
ggplot(featureRecoded272, aes(x = variable, y = length(unique()))) + geom_bar(position = "fill",stat = "identity")

feature134 = read.table("Desktop/test134-feature-table.tsv", header=T)
names(feature134) = str_replace_all(names(feature134), "^X", "")
names(feature134) = str_replace_all(names(feature134), "\\w+_(\\w+_\\w+)_.+$", "\\1")
featureMelt134 <- melt(feature134, id.vars = c('OTU'))
featureCompact134 = featureMelt134[which(featureMelt134$value > 1),]
featureRecoded134 = featureCompact134
featureRecoded134$OTU = as.character(1:length(featureRecoded134$OTU))
ggplot(featureRecoded134, aes(x = variable, y = value, fill = OTU)) + geom_bar(position = "fill",stat = "identity") + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))

feature152 = read.table("Desktop/test152-feature-table.tsv", header=T)
names(feature152) = str_replace_all(names(feature152), "^X", "")
names(feature152) = str_replace_all(names(feature152), "\\w+_(\\w+_\\w+)_.+$", "\\1")
featureMelt152 <- melt(feature152, id.vars = c('OTU'))
featureCompact152 = featureMelt152[which(featureMelt152$value > 1),]
featureRecoded152 = featureCompact152
featureRecoded152$OTU = as.character(1:length(featureRecoded152$OTU))
ggplot(featureRecoded152, aes(x = variable, y = value, fill = OTU)) + geom_bar(position = "fill",stat = "identity") + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))

feature166 = read.table("Desktop/test166-feature-table.tsv", header=T)
names(feature166) = str_replace_all(names(feature166), "^X", "")
names(feature166) = str_replace_all(names(feature166), "\\w+_(\\w+_\\w+)_.+$", "\\1")
featureMelt166 <- melt(feature166, id.vars = c('OTU'))
featureCompact166 = featureMelt166[which(featureMelt166$value > 1),]
featureRecoded166 = featureCompact166
featureRecoded166$OTU = as.character(1:length(featureRecoded166$OTU))
ggplot(featureRecoded166, aes(x = variable, y = value, fill = OTU)) + geom_bar(position = "fill",stat = "identity") + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))

feature179 = read.table("Desktop/test179-feature-table.tsv", header=T)
names(feature179) = str_replace_all(names(feature179), "^X", "")
names(feature179) = str_replace_all(names(feature179), "\\w+_(\\w+_\\w+)_.+$", "\\1")
featureMelt179 <- melt(feature179, id.vars = c('OTU'))
featureCompact179 = featureMelt179[which(featureMelt179$value > 100),]
featureRecoded179 = featureCompact179
featureRecoded179$OTU = as.character(1:length(featureRecoded179$OTU))
#ggplot(featureRecoded179, aes(x = variable, y = value, fill = OTU)) + geom_bar(position = "fill",stat = "identity") + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))


feature261 = read.table("Desktop/test261-feature-table.tsv", header=T)
names(feature261) = str_replace_all(names(feature261), "^X", "")
names(feature261) = str_replace_all(names(feature261), "\\w+_(\\w+_\\w+)_.+$", "\\1")
featureMelt261 <- melt(feature261, id.vars = c('OTU'))
featureCompact261 = featureMelt261[which(featureMelt261$value > 100),]
featureRecoded261 = featureCompact261
featureRecoded261$OTU = as.character(1:length(featureRecoded261$OTU))
ggplot(featureRecoded261, aes(x = variable, y = value, fill = OTU)) + geom_bar(position = "fill",stat = "identity") + scale_fill_manual(values=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"))

summary_plot <- function(x, name){
  x  %>% group_by(variable) %>% summarise(N=n()) %>% ggplot(aes(x=fct_reorder(variable, N), y=N/2)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
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


featureTables = list()
i=1
for(j in list.files("Desktop/FeatureTables/singlend")){
  tableJ = read.table(paste("Desktop/FeatureTables/singlend/", j, sep=""), header=T)
  names(tableJ) = str_replace_all(names(tableJ), "^X", "")
  names(tableJ) = str_replace_all(names(tableJ), "\\w+_(\\w+_\\w+)_.+$", "\\1")
  tableJ <- melt(tableJ, id.vars = c('OTU'))
  totalJ = length(tableJ$value)
  upper_cutoff = min(sort(tableJ$value)[c(totalJ, totalJ-1)])
  print(upper_cutoff)
  #tableJ = tableJ[which(tableJ$value > 100 & tableJ$value < upper_cutoff),]
  tableJ = tableJ[which(tableJ$value > 100),]
  #tableJ = tableJ[which(tableJ$value > 1),]
  featureTables[[i]] = tableJ
  i=i+1
}
names(featureTables) = str_replace_all(list.files("Desktop/FeatureTables/singlend"), "singlend(\\w+)-.+", "\\1")

i = 1
plotlist = list()
for(t in featureTables){
  plot_T <- summary_plot(t) + ggtitle(names(featureTables)[i])
  plotlist[[i]] <- plot_T
  i=i+1
}
multiplot(plotlist=plotlist, cols=3)

i = 1
plotlist2 = list()
for(t in featureTables){
  plotlist2[[i]] <- ggplot(t, aes(x = variable, y = log(value), fill = OTU)) + geom_bar(position = "fill",stat = "identity") + theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values=rep(c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"),10)) + ggtitle(names(featureTables)[i])
  i=i+1
}
multiplot(plotlist=plotlist2, cols=3)
#plot(plotlist2[[3]])

for(n in 1:length(featureTables)){
  print(names(featureTables)[n])
  featureTables[[n]][which(featureTables[[n]]$variable=="21_mix1"),] %>% .[order(.$value),] %>%  print()
}
featureTables$`272`[which(featureTables$`272`$variable=="25_Physo"),] %>% .[order(.$value),] 

