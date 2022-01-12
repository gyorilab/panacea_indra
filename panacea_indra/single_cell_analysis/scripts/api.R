#load libraries
# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
## Functions script
plot_density <- function(wd, all_data, fname, feature, intercept){
  
  png(filename = paste0(wd, 'output/images/', 'QC/',fname, '_',feature,'_density.png'), 
      res = 150, width = 1500, height = 1000)
  
  if(feature == 'percent_mito'){
    plot <- ggplot(all_data, aes(color=orig.ident, x=percent_mito, 
                                 fill= orig.ident))
  }else if(feature == 'nFeature_RNA'){
    plot <- ggplot(all_data, aes(color=orig.ident, x=nFeature_RNA, 
                                 fill= orig.ident))
  }else if(feature == 'nCount_RNA'){
    plot <- ggplot(all_data, aes(color=orig.ident, x=nCount_RNA, 
                                 fill= orig.ident))
  }
  plot <- plot + geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = intercept)
  
  print(plot)
  dev.off()
}

plot_vln <- function(wd, all_data, fname){
  # Visualize QC metrics as a violin plot and save as PNG
  png(filename = paste0(wd, 'output/images/', 'QC/',fname, '_','vln_plot.png'), 
      res = 150, width = 1500, height = 1000)
  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
  print(VlnPlot(all_data, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
          NoLegend())
  dev.off()
}

plot_ncells <- function(wd, fname, obj){
  
  png(paste0(wd, '/output/images/ncells/',fname, '.png'),
      width = 1000, height = 1000, res=150)
  print(ggplot(obj, aes(color=orig.ident, x=unname(obj@active.ident),
                        fill= orig.ident)) + geom_bar() + xlab('') +
          geom_text(stat='count', aes(label=..count..), vjust=-1))
  dev.off()
}

volcano_generate <-function(res, titlename, file_loc){
  keyvals <- ifelse(
    res$log2FoldChange < 0 & res$padj < 0.05, 'blue3',
    ifelse(res$log2FoldChange > 0 & res$padj < 0.05, 'red3',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red3'] <- 'Upregulated'
  names(keyvals)[keyvals == 'black'] <- 'NS'
  names(keyvals)[keyvals == 'blue3'] <- 'Downregulated'
  pngname <- paste(file_loc,".png",sep="")
  
  
  png(pngname, width=2000,height=1500, res=150)
  
  print(EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        FCcutoff = FALSE,
                        colCustom = keyvals,
                        colAlpha = 1,
                        #legendPosition = 'right',
                        title = titlename,
                        pointSize = 3.0,
                        labSize = 6.0,
                        
                        xlim = c(-8, 8)))
  dev.off() 
}