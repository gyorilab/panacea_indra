# Libraries
library(tm)
library(proustr)
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
    main = title,
    
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
    cat.pos = c(-11, 10),
    cat.dist = cat_dist,
    cat.fontfamily = "sans",
    #rotation = 1
  )
  
  
}

# Read the datasets
setwd("~/gitHub/panacea_indra/panacea_exam/transcriptomics/inputs/")

all_kinases <- read.csv('./kinase.csv', row.names = 1)

df <- readxl::read_xlsx('./List_of_genes_for_heatmaps.xlsx')

# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(2, "Pastel2")

# Creating generalized sets
set_list_1 <- list_1[,c(1,8)]
set_list_2 <- list_2[,c(1,4)]



# common genes in kinases and hDRG kinases
type <- 'All Kinase vs hDRG Kinases'
cat_dist <- rep(0.02,5)
dir.create(paste0('../outputs//venn_diagrams/'),
           showWarnings = F, 
           recursive = T)

loc = paste0('../outputs/venn_diagrams/all_kinase_vs_hDRG_kinase.tiff')
all_genes <- list(all_kinases$x, df$`Gene ste5: hDRG enriched kinases`)
cat_names <- c('all_kinases', 'hDRG_kinases')
draw_venn(all_genes, loc, type, cat_names, myCol[1:2], cat_dist[1:2])



