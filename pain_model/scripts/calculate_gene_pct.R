library(tidyverse)
library(stringr)
library(Seurat)
library(dplyr)

source('~/gitHub/panacea_indra/pain_model/scripts/functions.R')

## Notes

# calculate fraction of the gene in a cluster:
# Number of single nuclei/cells expressing gene in the cluster
# /
# Total number of single nuclei/cells in the cluster

# calculate percentage of the gene expressed in a cluster
# 
#---X---

wd <- ('/Users/sbunga/gitHub/panacea_indra/pain_model/')
setwd(wd)
human_mouse <- read.csv('./data/human_mouse_orth.txt', sep='\t')
colnames(human_mouse) <- c('HUMAN_SYMBOL', 'MOUSE_SYMBOL')

df <- readxl::read_xlsx(paste0(wd, 'data/Primary_mouse/other/mouse_10-fold_high_genes_expression_matrix.xlsx'))
enriched_genes <- df[,c(1,2,3,4)]
colnames(enriched_genes) <- c('HUMAN_SYMBOL', 'Description', 'Subcellular_localization',
                              'Functional_type_of_protein')
enriched_genes  <- merge(enriched_genes, human_mouse, by.y ='HUMAN_SYMBOL')
enriched_genes <- enriched_genes[!enriched_genes[ , 5] == '', ]
enriched_genes <- enriched_genes[!duplicated(enriched_genes$MOUSE_SYMBOL), ]
enriched_genes$HUMAN_SYMBOL <- NULL
enriched_genes <- enriched_genes[, c(4, 1, 2, 3)]



# Visual cortex data
vis_cortex <- read.csv('./data/Primary_mouse/visual_cortex/GSE71585_RefSeq_TPM.csv')
vis_cortex_meta <- read.csv('./data/Primary_mouse/visual_cortex/GSE71585_Clustering_Results.csv')

vis_cortex <- vis_cortex %>% column_to_rownames('gene')
vis_cortex <- CreateSeuratObject(vis_cortex)
vis_cortex$cell_type <- vis_cortex_meta$primary_type
vis_cortex_enriched_mat <- calculate_pct(vis_cortex, enriched_genes)
write.csv(vis_cortex_enriched_mat, './output/gene_pct_visual_cortex_clusters.csv', row.names = F)


# DRG neuron data from Will's lab
mDRG <- readRDS(paste0(wd, 'data/Primary_mouse/will_renthal_mDRG/GSE154659_C57_Raw_counts.RDS'))
mDRG <- CreateSeuratObject(mDRG, project = 'mDRG')
mDRG$cell_type <- NA
mDRG$cell_type <- str_match(rownames(mDRG@meta.data), "(.*rep\\d_)([A-Za-z0-9\\s^_]+)(_.*)")[,3]

mDRG_enriched_mat <- calculate_pct(mDRG, enriched_genes)
write.csv(mDRG_enriched_mat, './output/gene_pct_drg_clusters.csv', row.names = F)



# un-mapped genes
#unmapped_genes <- df$...1[sapply(df$...1, function(x) !x %in% enriched_genes$HUMAN_SYMBOL)]
#print(unmapped_genes)



# Mouse Heart
mheart <- read.csv('./data/Primary_mouse/mHeart/Heart-counts.csv', row.names = 1)
mheart_meta <- read.csv('./data/Primary_mouse/mHeart/metadata_FACS.csv')
mheart_anno <- read.csv('./data/Primary_mouse/mHeart/annotations_FACS.csv')
# Filter annotations
mheart_anno <- mheart_anno[mheart_anno$cell %in% colnames(mheart), ]
mheart <- CreateSeuratObject(mheart)

# Filter out cells which are not in the annotation table
mheart <- mheart[,mheart_anno$cell]

# rearrange the cell ids in metadat
mheart@meta.data <- mheart@meta.data[mheart_anno$cell,]
mheart@meta.data$cell_type <- mheart_anno$free_annotation
mheart_enriched_mat <- calculate_pct(mheart, enriched_genes)
write.csv(mheart_enriched_mat, './output/gene_pct_heart_clusters.csv',
          row.names = F)  

# Motor cortex
mcortex <- Read10X('./data/Primary_mouse/motor_cortex/', 
                   gene.column = 2)
