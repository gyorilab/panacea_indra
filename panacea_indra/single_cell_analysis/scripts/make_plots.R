wd <- '~/gitHub/panacea_indra/panacea_indra/single_cell_analysis/'
source(paste0(wd, '/scripts/api.R'))

seurat_integrated <- readRDS(paste0(wd, '/RDS/seurat_integrated_0.6.RDS'))
seurat_integrated.markers <- readRDS(paste0(wd, '/RDS/seurat_integrated_markers.Rds'))
seurat_integrated_unknown_removed <- readRDS(paste0(wd,'/RDS/seurat_integrated_unknown_removed.Rds'))

dir.create(paste0(wd, 'output/images/res_0.6'), showWarnings = F)


seurat_integrated.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

gene_list <- c("S100a9", "Csf3r", "H2-Ab1", "Cd74","Lyz2","Csf1r","Il1b", "Cxcl1",
               "F13a1", "Selenop", "Cd163", "Mrc1", "Chil3","Fn1", "Icos", "Thy1", 
               "Mcpt4", "Krtdap")

tiff(paste0(wd,'output/images/res_0.6/dotplot_gene_list_rna_assay.png'), res=200, 
    height = 1200, width = 2900)
print(DotPlot(seurat_integrated, features = gene_list, assay = 'RNA', scale = T))
dev.off()

# UMAP of all the integrated clusters

png(paste0(wd, '/output/images/res_0.6/umap_orig_ident_all_samples.png'), width = 2500,
    height = 1900, res=200)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = TRUE,
              label.size = 4))
dev.off()

# UMAP group by replicates
png(paste0(wd, '/output/images/res_0.6/umap_replicate_all_samples.png'), width = 2500,
    height = 1900, res=200)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = TRUE,
              label.size = 4, split.by = 'replicate'))
dev.off()


# Show the top 10 genes in each cluster
png(paste0(wd, '/output/images/res_0.6/top_10_genes_all_clusters.png'), width = 3750,
    height=4500, res=200)
print(DoHeatmap(object = seurat_integrated, features = top10$gene) + NoLegend())
dev.off()




# Remove unknown cluster and subset injured
injury_subset <- subset(seurat_integrated, type == 'Injured')
injury_subset <- subset(injury_subset, idents='Unknown', invert=T)

# Plot UMAP
png(paste0(wd, '/output/images/res_0.6/umap_model_injured.png'), width = 3500,
    height = 1900, res=200)
print(DimPlot(injury_subset,
              reduction = "umap",
              label = TRUE,
              label.size = 4, split.by = 'model'))
dev.off()



################################################################################


# Highlight interesting genes
png(paste0(wd, '/output/images/heatmap/ditto_highlighted_features.png'), width = 1250,
    height=1000, res=100)
print(dittoHeatmap(seurat_integrated_unknown_removed, gene_list, annot.by = c("cell_type"),
                   #highlight.features = gene_list,
                   complex = TRUE, scaled.to.max = T))

dev.off()



# subset HT
#seurat_HT <- subset(seurat_integrated, type=='HT')

# subet Stim
#seurat_stim <- subset(seurat_integrated, type=='STIM')



# Plot UMAP
dir.create(paste0(wd, '/output/images/UMAP'), showWarnings = F)



# Final images
# ------ X ------


# UMAPs
png(paste0(wd, '/output/images/UMAP/umap_model.png'), width = 5500,
    height = 3500, res=300)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = F,
              label.size = 4, split.by = 'model'))
dev.off()


png(paste0(wd, '/output/images/UMAP/umap_cluster_no_label.png'), width = 4500,
    height = 3500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = F,
              label.size = 4))
dev.off()


png(paste0(wd, '/output/images/UMAP/umap_all_samples_model_condition.png'), width = 3500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = F,
              label.size = 4, split.by = 'model', group.by = 'condition'))
dev.off()





# No. of cells in each condition
table(seurat_integrated@meta.data$condition)

# No. of cells in each cluster
table(seurat_integrated$integrated_snn_res.0.6)

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"


# Markers
marker_genes <- c('Icos', 'Trdc', 'Trbc', 'Csf1r', 'Mrc1', 'Fcgr1a', 'H2-Ab1',
                  'Itgax', 'Cd207', 'Kit', 'Cd14', 'S100a9', 'Cx3cr1', 'Maf')

png(paste0(wd, '/output/images/Feature_plots/new_genes_1_all_samples.png'), width = 4500,
    height = 2000, res=400)
print(FeaturePlot(seurat_integrated, features = marker_genes))
dev.off()


# Select assay back to integrated
DefaultAssay(seurat_integrated) <- "integrated"
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"


# Get stats for each cluster
stats <- seurat_integrated %>% dplyr::select(integrated_snn_res.0.6, condition) %>% 
  dplyr::group_by(integrated_snn_res.0.6,condition) %>% dplyr::summarise(n())

colnames(stats)[3] <- 'count'
write.csv(stats, paste0(wd, 'output/cluster_0.6_condition_stats_all_samples.csv'))
