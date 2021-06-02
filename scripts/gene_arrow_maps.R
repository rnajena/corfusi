###### visualization: gene arrow maps for cORFusi ######

library(rtracklayer)
library(stringr)
library(ggplot2)
library(gggenes)
library(dplyr)
library(scales)
library(grid)
library(gridExtra)
library(tidyverse)


############ input variables ############
###### PLEASE UPDATE FOR YOUR DATA ######

log <- "path/to/log/file"
assembly <- "path/to/hybrid/assembly/gff"
short <- "path/to/short-read/assembly/gff"
corrected <- "path/to/corrected/assembly/gff"
save_path <- "path/to/figures/"





############ path and variables ############
log.data <- read.table(file = log, sep = '\t', header = TRUE)

hybrid.data <- readGFF(assembly)
short_read.data <- readGFF(short)
hybrid_corr.data <- readGFF(corrected)

header <- c("genome", "start", "end", "gene", "strand")
corr_genes <- list()

hybrid.infos <- data.frame()
short_read.infos <- data.frame()
hybrid_corr.infos <- data.frame()



############ extract genes and and positions ############
for(i in 1:nrow(log.data)){
  id <- log.data[i,1]
  start <- log.data[i,2]
  c <- sapply('ID', function(x) grep(id, short_read.data[,x]))
  corrected_gene <- short_read.data[c,]$Name
  contig <- short_read.data[c,]$seqid
  
  if(is.na(corrected_gene) == FALSE){
    corr_genes <- c(corr_genes, list(corrected_gene))
    before <- FALSE
    after <- FALSE
    
    #short read assembly
    for(i in (c-1):(c+1)){
      if(contig != short_read.data[i,]$seqid){
        if(i < c){before <- TRUE}
        if(i > c){after <- TRUE}
        break}
      if(is.na(short_read.data[i,]$Name) == TRUE){
        gene <- "unknown"
      } else{
        gene <- short_read.data[i,]$Name
      }
      if(short_read.data[i,]$strand == "+"){
        fw <- 1
      } else {
        fw <- -1
      }
      h <- data.frame(id, short_read.data[i,]$start, short_read.data[i,]$end, gene, fw)
      short_read.infos <- rbind.data.frame(short_read.infos, h)
    }
    
    #hybrid assembly
    fragments <- hybrid.data[grep(corrected_gene, hybrid.data$Name), ]
    d <- sapply('Name', function(x) grep(fragments[1,]$Name, hybrid.data[,x]))
    
    for(i in (d-1):(d+nrow(fragments)+1)){
      if(i == (d-1) & before == TRUE){break}
      if(i == (d+nrow(fragments)) & after == TRUE){break}
      if(is.na(hybrid.data[i,]$Name) == TRUE){
        gene <- "unknown"
      } else{
        gene <- hybrid.data[i,]$Name
      }
      if(hybrid.data[i,]$strand == "+"){
        fw <- 1
      } else {
        fw <- -1
      }
      h <- data.frame(id, hybrid.data[i,]$start, hybrid.data[i,]$end, gene, fw)
      hybrid.infos <- rbind.data.frame(hybrid.infos, h)
    }
    
    #corrected assembly
    fragments <- hybrid_corr.data[grep(corrected_gene, hybrid_corr.data$Name), ]
    d <- sapply('Name', function(x) grep(fragments[1,]$Name, hybrid_corr.data[,x]))
    
    for(i in (d-1):(d+nrow(fragments)+1)){
      if(i == (d-1) & before == TRUE){break}
      if(i == (d+nrow(fragments)) & after == TRUE){break}
      if(is.na(hybrid_corr.data[i,]$Name) == TRUE){
        gene <- "unknown"
      } else{
        gene <- hybrid_corr.data[i,]$Name
      }
      if(hybrid_corr.data[i,]$strand == "+"){
        fw <- 1
      } else {
        fw <- -1
      }
      h <- data.frame(id, hybrid_corr.data[i,]$start, hybrid_corr.data[i,]$end, gene, fw)
      hybrid_corr.infos <- rbind.data.frame(hybrid_corr.infos, h)
    }
  }
}


colnames(hybrid.infos) <- header
colnames(short_read.infos) <- header
colnames(hybrid_corr.infos) <- header


############ plot genes and neighborhood ############
for(id in log.data[,1]){
  j <- short_read.infos[grep(id, short_read.infos$genome), ]
  if(nrow(j) == 0) next
  k <- hybrid.infos[grep(id, hybrid.infos$genome), ]
  l <- hybrid_corr.infos[grep(id, hybrid_corr.infos$genome), ]
  levels(j$genome) <- sub(id, "Short-read assembly", levels(j$genome))
  levels(k$genome) <- sub(id, "Hybrid assembly", levels(k$genome))
  levels(l$genome) <- sub(id, "Corrected assembly", levels(l$genome))
  
  final <- rbind.data.frame(j,k,l)
  
  c <- sapply('ID', function(x) grep(id, short_read.data[,x]))
  name <- short_read.data[c,]$Name
  
  ggplot(final, aes(xmin = start, xmax = end, y = genome, fill = gene, label = gene, forward = strand)) +
    geom_gene_arrow(arrowhead_width = grid::unit(8, "mm"), arrowhead_height = grid::unit(7, "mm"), arrow_body_height = grid::unit(5, "mm")) +
    geom_gene_label(align = "centre", padding.x = grid::unit(0.7,"mm"), padding.y = grid::unit(0.2,"mm")) +
    facet_wrap(~ genome, scales="free", ncol = 1) +
    scale_fill_brewer(palette = "Set3") +
    theme_genes() +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.0), legend.position = "none") +
    ylab("") #+ ggtitle(paste("ORF correction of ",name, " gene in strain ", s, sep = ""))
  
  ggsave(paste(save_path, "/", name, ".pdf", sep = ""),
         width = 8, height = 3, dpi = 300, units = "in")
  
}