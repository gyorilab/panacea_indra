library(stringr)
library(Seurat)
library(dplyr)

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
df <- readxl::read_xlsx(paste0(wd, 'data/Primary_mouse/other/mouse_10-fold_high_genes_expression_matrix.xlsx'))
mDRG <- readRDS(paste0(wd, 'data/Primary_mouse/will_renthal_mDRG/GSE154659_C57_Raw_counts.RDS'))
mDRG <- CreateSeuratObject(mDRG, project = 'mDRG')
mDRG$cell_type <- NA
mDRG$cell_type <- str_match(rownames(hDRG@meta.data), "(.*rep\\d_)([A-Za-z0-9\\s^_]+)(_.*)")[,3]

# Subsetting the mDRG object and extracting the cells
# from each cluster
all_cell_types <- names(table(mDRG@meta.data$cell_type))
all_clusters = list()

enriched_genes <- df[,c(1,2,3,4)]
colnames(enriched_genes) <- c('Gene_name', 'Description', 'NA1', 'NA2')
for(n in all_cell_types){
  enriched_genes[, n] <- NA
}

for(n in 1:length(all_cell_types)){
  all_clusters[[all_cell_types[n]]] = subset(mDRG, cell_type == 
                                               all_cell_types[n])
}

for(g in 1:nrow(enriched_genes)){
  for(cluster in all_cell_types){
    exp <- FetchData(all_clusters[[cluster]], enriched_genes$Gene_name[g])
    # Calculating the percentage
    pct <- as.matrix(colMeans(exp  > 0))*100
    enriched_genes[g, cluster] <- pct
  }
}
