library(stringr)
library(Seurat)

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
mouse_genes_enrch <- readxl::read_xlsx(paste0(wd, 'data/Primary_mouse/other/mouse_10-fold_high_genes_expression_matrix.xlsx'))
hDRG <- readRDS(paste0(wd, 'data/Primary_mouse/will_renthal_hDRG/GSE154659_C57_Raw_counts.RDS'))
hDRG <- CreateSeuratObject(hDRG, project = 'hDRG')
hDRG$cell_type <- NA
hDRG$cell_type <- str_match(rownames(hDRG@meta.data), "(.*rep\\d_)([A-Za-z0-9\\s^_]+)(_.*)")[,3]

View(hDRG@meta.data)
cd_genes <- c("Trpv1")
a <- DotPlot(object = hDRG, features = cd_genes, group.by = '')
