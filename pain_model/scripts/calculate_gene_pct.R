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

# get enriched genes
enriched_genes <- get_enriched_genes()

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
#mcortex <- Read10X('./data/Primary_mouse/motor_cortex/', 
#                   gene.column = 2)
mcortex <- readRDS('./data/Primary_mouse/motor_cortex/mcortex.Rds')
cluster_memb <- read.csv('./data/Primary_mouse/motor_cortex/cluster.membership.csv')
meta_file <- read.csv('./data/Primary_mouse/motor_cortex/cluster.annotation.csv')
colnames(cluster_memb)[2] <- 'cluster_id'
#meta_file$cluster_label
cluster_memb <- merge(cluster_memb, meta_file, by.x = 'cluster_id') 
# Subset the cells based on the meta file
#mcortex$X <- colnames(mcortex)
#mcortex <- mcortex[, mcortex$X %in% cluster_memb$X]

mcortex@meta.data = merge(mcortex@meta.data, cluster_memb, by = 'X')
colnames(mcortex@meta.data)[6] <- 'cell_type'
rownames(mcortex@meta.data) <- mcortex@meta.data$X
mcortex_enriched_mat <- calculate_pct(mcortex, enriched_genes)
write.csv(mcortex_enriched_mat, './output/gene_pct_cortex_clusters.csv',
          row.names = F)

# Mouse motor neurons
mneurons <- read.csv('./data/Primary_mouse/mMotor_neurons/GSE103892_Expression_Count_Matrix.txt.gz',
                     sep='\t')
annot <- read.csv('./data/Primary_mouse/mMotor_neurons/GSE103892_Sample_Cell_Cluster_Information.txt',
                  sep='\t')
rownames(mneurons) <- mneurons$Gene
mneurons$Gene <- NULL
mneurons <- CreateSeuratObject(mneurons)

View(mneurons@meta.data)
mneurons$sample_cellbarcode <- rownames(mneurons@meta.data)
mneurons@meta.data <- merge(mneurons@meta.data, annot, by='sample_cellbarcode')
rownames(mneurons@meta.data) <- mneurons@meta.data$sample_cellbarcode
colnames(mneurons@meta.data)[9] <- 'cell_type'

mneurons_enriched_mat <- calculate_pct(mneurons, enriched_genes)
write.csv(mneurons_enriched_mat, './output/gene_pct_motor_neuron_clusters.csv',
          row.names = F)


