# Venn diagrams
# Libraries
library(tm)
library(dplyr)
library(proustr)
library(stringr)
library(tidyverse)
library(hrbrthemes)
library(VennDiagram)
library(RColorBrewer)

# Functions
# Chart
draw_venn <- function(set_list, file_loc, title, cat_names, mycol, cat_dist){
  
  venn.diagram(
    x = set_list,
    category.names = cat_names,
    filename = paste0(file_loc),
    output=T,
    main = '',
    
    # Output features
    imagetype="tiff" ,
    height = 1200, 
    width = 1500 , 
    resolution = 250,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = mycol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.5,
    cat.fontface = "bold",
    #cat.default.pos = "inner",
    #cat.pos = c(0, 1, 2),
    #cat.dist = cat_dist,
    cat.fontfamily = "sans",
    #rotation = 1
  )
  
  
}



# Read the datasets
setwd("./output/pseudo_bulk_analysis/")
dir.create('./plots/venn_diagrams', showWarnings = F)
identities <- list()
clusters <- c('RM', 'Macs4', 'Dermal Macs')
comparisons <- c('Zymo_Saline', 'Incision_Healthy', 'UVB_Sham')
for(cluster in clusters){
  files <- na.omit(str_extract(list.files('./diff_files/'), paste0(cluster, '.*')))
  gene_sets <- list()
  for(cc in comparisons){
    for(f in files){
    if(str_detect(f, paste0(".*",cc,".*"))){
      df <- read.csv(paste0('./diff_files/', str_extract(f, paste0(".*",cc,".*"))))
      filtered_genes <- df %>% filter(log2FoldChange > 1 & padj < 0.05) 
      gene_sets[[cc]] <- filtered_genes$X
      }
    }
  }
  # Prepare a palette of 3 colors with R colorbrewer:
  myCol <- brewer.pal(3, "Pastel2")
  
  loc = paste0('./plots/venn_diagrams/up_reg_',cluster,'.png') 
  all_genes <- gene_sets
  cat_names <- names(gene_sets)
  draw_venn(all_genes, loc, type, cat_names, myCol[1:3], cat_dist[1:3])
  identities[[cluster]] <- intersect(gene_sets$Zymo_Saline, 
                                   c(gene_sets$Incision_Healthy, 
                                     gene_sets$UVB_Sham))
}





