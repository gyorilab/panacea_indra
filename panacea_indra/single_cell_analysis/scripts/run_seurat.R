library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(dittoSeq)
library(tidyseurat)
library(ComplexHeatmap)

## Notes
# plot barplots of no. of cells of all samples before and after filtering

# Visualize QC as density plots
# nFeature
# trial -> compare cluster 2 (Zymo, UVB, Incision)

wd <- '/Users/sbunga/gitHub/singlecell_model_analysis/'
source(paste0(wd, '/scripts/functions.R'))

dir.create(paste0(wd, '/output/images/ncells/'), showWarnings = F)

input_dirs <- list.dirs(paste0(wd, 'input_files'), full.names = T)
input_files <- na.omit(str_extract(input_dirs, ".*\\/[ATGC]+"))

# Create outputs
paths <- c('images')
for(p in paths){
  dir.create(paste0(wd, '/output/', p), showWarnings = F, recursive = T)
}
dir.create(paste0(wd, 'output/images/', 'QC'), showWarnings = F, recursive = T)
dir.create(paste0(wd, 'RDS'), showWarnings = F)

sample_names_list <- list(
  'CTAGTCGA'='Zymo_1',
  'TCTTACGC'='Zymo_2',
  'AGCTAGAA'='Saline_1',
  'ACTCTAGG'='Saline_2',
  'TACTCCTT'='Sham_1',
  'ATTAGACG'='Sham_2',
  'AGGCTTAG'='UVB_1',
  'CGGAGAGA'='UVB_2',
  'ATAGAGAG'='Healthy_1',
  'TATGCAGT'='Healthy_2',
  'AGAGGATA'='Incision_1',
  'CTCCTTAC'='Incision_2'
)

# Give sample types
HT <- c('Healthy_1', 'Healthy_2', 'Saline_1', 'Saline_2', 'Sham_1', 'Sham_2')
STIM <- c('Zymo_1', 'Zymo_2', 'Incision_1', 'Incision_2', 'UVB_1', 'UVB_2')

models <- list(
  'Incision'='Incision',
  'Healthy'='Incision',
  'Zymo'='Zymo',
  'Saline'='Zymo',
  'UVB'='UVB',
  'Sham'='UVB'
)

condition <- list(
  'Incision_1' = 'Incision',
  'Incision_2' = 'Incision',
  'Healthy_1' = 'Healthy',
  'Healthy_2' = 'Healthy',
  'Zymo_1' = 'Zymo',
  'Zymo_2' = 'Zymo',
  'Saline_1' = 'Saline',
  'Saline_2' = 'Saline',
  'UVB_1' = 'UVB',
  'UVB_2' = 'UVB',
  'Sham_1' = 'Sham',
  'Sham_2' = 'Sham'
)

# Read samples
input_files <- na.omit(input_files[basename(input_files) %in% 
                                     names(sample_names_list)])
infiles <- lapply(input_files, Read10X)

# Get sample names in order from input paths
barcodes <- basename(input_files)
sample_names <- unname(unlist(sample_names_list[barcodes]))

# Create Seurat Objects and add meta data
object_list <- list()
for (s in 1:length(sample_names)) {
  object_list[[s]] <- CreateSeuratObject(infiles[[s]], project = sample_names[s])
  # Check for sample type
  type <- if(sample_names[s] %in% HT) 'HT' else 'STIM'
  object_list[[s]]$type <- type
  # add replicate in the metadata
  object_list[[s]]$replicate <- str_split(sample_names[s], '_', simplify = F)[[1]][2]
  # add sample model
  object_list[[s]]$model <- models[str_split(sample_names[s], '_', simplify = F)[[1]][1]]
  object_list[[s]]$condition <- condition[sample_names[s]][[1]]
  }

# Merge Incision
incision_samples <- c(1, 3, 2, 4)
incision <- object_list[incision_samples]
incision <- merge(incision[[1]], incision[-1],
                  all.cell.ids = sample_names[incision_samples])
incision <- PercentageFeatureSet(incision, "^mt-", col.name = "percent_mito")

plot_density(wd, incision, 'incision', 'nFeature_RNA', 300)
plot_density(wd, incision, 'incision', 'nCount_RNA', 500)
plot_density(wd, incision, 'incision', 'percent_mito', 300)
plot_vln(wd, incision, 'all_incision')

dir.create(paste0(wd, '/output/images/Feature_plots'), showWarnings = F)
png(paste0(wd, '/output/images/Feature_plots/incision.png'),
    width = 1500, height = 1500, res=200)
VlnPlot(incision, features = c("Thbs1","Csf1r","Ptgs2","Il1b"))
dev.off()

saveRDS(incision, paste0(wd, '/RDS/incision.Rds'))

# n cells before filtering
plot_ncells(wd,'incision_before_filter', incision)
incision <- subset(incision, subset = nFeature_RNA > 300 & 
                     nFeature_RNA < 2500 & nCount_RNA > 500 & percent_mito < 5)
plot_ncells(wd, 'incision_after_filter', incision)


# Merge UVB
#uvb_samples <- c(5, 7, 6)
uvb_samples <- c(5, 7, 8, 6)
uvb <- object_list[uvb_samples]
uvb <- merge(uvb[[1]], uvb[-1],
             all.cell.ids = sample_names[uvb_samples])

uvb <- PercentageFeatureSet(uvb, "^mt-", col.name = "percent_mito")

plot_density(wd, uvb, 'UVB', 'nFeature_RNA', 250)
plot_density(wd, uvb, 'UVB', 'nCount_RNA', 300)
plot_density(wd, uvb, 'UVB', 'percent_mito', 300)
plot_vln(wd, uvb, 'all_uvb')

png(paste0(wd, '/output/images/Feature_plots/uvb.png'),
    width = 1500, height = 1500, res=200)
VlnPlot(uvb, features = c("Thbs1","Csf1r","Ptgs2","Il1b"))
dev.off()

saveRDS(uvb, paste0(wd, '/RDS/uvb.Rds'))

plot_ncells(wd, 'uvb_before_filter', uvb)
uvb <- subset(uvb, subset = nFeature_RNA > 250 &
                nFeature_RNA < 2200 & nCount_RNA > 300 & percent_mito < 5)
plot_ncells(wd, 'uvb_after_filter', uvb)

# Merge Zymo
#zymo_sample <- c(9, 8, 10, 11)
zymo_sample <- c(10, 9 , 11, 12)
zymo <- object_list[zymo_sample]
zymo <- merge(zymo[[1]], zymo[-1],
              all.cell.ids = sample_names[zymo_sample])

zymo <- PercentageFeatureSet(zymo, "^mt-", col.name = "percent_mito")

plot_density(wd, zymo, 'ZYMO', 'nFeature_RNA', 300)
plot_density(wd, zymo, 'ZYMO', 'nCount_RNA', 500)
plot_density(wd, zymo, 'ZYMO', 'percent_mito')
plot_vln(wd, zymo, 'all_zymo')

png(paste0(wd, '/output/images/Feature_plots/zymo.png'),
    width = 1500, height = 1500, res=200)
VlnPlot(zymo, features = c("Thbs1","Csf1r","Ptgs2","Il1b"))
dev.off()

saveRDS(zymo, paste0(wd, '/RDS/zymo.Rds'))
plot_ncells(wd, 'zymo_before_filter', zymo)
zymo <- subset(zymo, subset = nFeature_RNA > 300 & 
                     nFeature_RNA < 2200 & nCount_RNA > 500 & percent_mito < 5)
plot_ncells(wd, 'zymo_after_filter', zymo)

# Merge objects
#all_data <- merge(object_list[[1]], object_list[-1], 
#                  add.cell.ids = sample_names )


# Merge all the samples
all_data <- merge(incision, c(uvb, zymo), add.cell.ids = c("incision", "uvb",
                                                           "zymo"))

saveRDS(all_data, paste0(wd, '/RDS/merged_object.Rds'))

healthy_merged <- subset(all_data, type == 'HT')
saveRDS(healthy_merged, paste0(wd, '/RDS/healthy_merged.Rds'))
stim_merged <- subset(all_data, type == 'STIM')
saveRDS(stim_merged, paste0(wd, '/RDS/stim_merged.Rds'))

# SCTransform
split_seurat <- SplitObject(all_data, split.by = "condition")

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent_mito"))
}

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
saveRDS(integ_anchors, 
        paste0(wd, '/output/integrated_anchors.Rds'))



# Integrate across conditions
#seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
#                                   normalization.method = "SCT")
seurat_integrated <- readRDS(paste0(wd, 
                                    './RDS/seurat_integrated_all_samples.Rds'))

#DefaultAssay(seurat_integrated) <- 'RNA'
#seurat_integrated <- NormalizeData(seurat_integrated)
#vfeatures <- FindVariableFeatures(seurat_integrated, selection.method = "vst", nfeatures = 3000)
#top10 <- head(VariableFeatures(vfeatures), 10)
#plot1 <- VariableFeaturePlot(vfeatures)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#png('./gitHub/singlecell_model_analysis/output/images/QC/variable_features.png',
#    height=1000, width=1000, res=150)
#print(plot2)
#dev.off()


# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
png(paste0(wd, '/output/images/QC/pca_split_orig_ident_all_samples.png'), width = 1050,
    height=1000, res=100)
print(PCAPlot(seurat_integrated,
        split.by = "orig.ident"))
dev.off()

png(paste0(wd, '/output/images/QC/pca_all_samples.png'), width = 1050,
    height=1000, res=100)
print(PCAPlot(seurat_integrated))
dev.off()

png(paste0(wd, '/output/images/QC/elbow_plot.png'), width = 1050,
    height=1000, res=100)
print(ElbowPlot(seurat_integrated, ndims = 60))
dev.off()
# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")


# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 50)

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.6)
saveRDS(seurat_integrated, paste0(wd, '/RDS/seurat_integrated_clustered.Rds'))

seurat_integrated <- readRDS(paste0(wd,'/RDS/seurat_integrated_clustered.Rds'))

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"

# Find all the DE genes from all clusters
#seurat_integrated.markers <- FindAllMarkers(seurat_integrated, 
#                                            only.pos = TRUE,
#                                            min.pct = 0.25, 
#                                            logfc.threshold = 0.25)
#saveRDS(seurat_integrated.markers, paste0(wd, '/RDS/seurat_integrated_markers.Rds'))

seurat_integrated.markers <- readRDS(paste0(wd,'/RDS/seurat_integrated_markers.Rds'))


# Annotating clusters
clusters.ids <- c("Neutrophils", "DC2", "Macs4", "Dermal Macs", "Macs3",
                  "T cells", "RM","Macs1", "Mast cells", "DC1", "LCs", "Macs2",
                  "Keratinocytes")
names(clusters.ids) <- levels(seurat_integrated)[c(1:length(clusters.ids))]
seurat_integrated <- RenameIdents(seurat_integrated, clusters.ids)

#cluster_cells <- list(
#  'Neutrophil' = WhichCells(seurat_integrated, idents = c("Neutrophils")),
#  'RM1' = WhichCells(seurat_integrated, idents = c("RM1")),
#  'RM2' = WhichCells(seurat_integrated, idents = c("RM2")),
#  'Dermal Macs' = WhichCells(seurat_integrated, idents = c("Dermal Macs"))
#)

seurat_integrated$cell_type <- Idents(seurat_integrated)

gene_list <- c("S100a9", "Csf3r", "Mgl2", "H2-Ab1", "Cd74", "Chil3", "Fn1", 
               "Fcgr1", "Lyz2", "F13a1", "Selenop", "Cd163", "Mrc1", "Csf1r", 
               "Il2rb", "Thy1", "Ptgs2", "Tnf", "Il1b", "Cxcl1", "Thbs1", "Fcer1a")

seurat_integrated.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


# Highlight certain genes
png(paste0(wd, '/output/images/heatmap/highlighted_features.png'), width = 1250,
    height=1500, res=100)
print(dittoHeatmap(seurat_integrated, top10$gene, annot.by = c("cell_type"),
             highlight.features = gene_list,
             complex = TRUE, assay = 'integrated', slot='data',scaled.to.max=T)
)
dev.off()

DoHeatmap(object = seurat_integrated, features = top10$gene)
dev.off()
# subset HT
#seurat_HT <- subset(seurat_integrated, type=='HT')

# subet Stim
#seurat_stim <- subset(seurat_integrated, type=='STIM')



# Plot the UMAP
dir.create(paste0(wd, '/output/images/UMAP'), showWarnings = F)



# Final images
# ------ X ------


# Heatmaps
dir.create(paste0(wd, '/output/images/heatmap/'), showWarnings = F)
png(paste0(wd, '/output/images/heatmap/top_20_subset_clusters.png'), width = 3500,
    height = 4200, res=400)
cluster_subset.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
print(DoHeatmap(cluster_subset, features = top20$gene) + NoLegend())
dev.off()

png(paste0(wd, '/output/images/heatmap/top_10_all_clusters.png'), width = 3500,
    height = 7200, res=400)
seurat_integrated.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
print(DoHeatmap(seurat_integrated, features = top10$gene) + NoLegend())
dev.off()

# UMAPs
png(paste0(wd, '/output/images/UMAP/umap_HT_model.png'), width = 3500,
    height = 2500, res=400)
print(DimPlot(seurat_HT,
              reduction = "umap",
              label = F,
              label.size = 4, split.by = 'model'))
dev.off()


png(paste0(wd, '/output/images/UMAP/umap_sub_cluster_no_label.png'), width = 3500,
    height = 2500, res=400)
print(DimPlot(cluster_subset,
              reduction = "umap",
              label = F,
              label.size = 4))
dev.off()



png(paste0(wd, '/output/images/UMAP/umap_sub_cluster_type_condition.png'), width = 3500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = F,
              label.size = 4, split.by = 'condition',
              group.by = 'model', cells.highlight = cluster_cells,
              cols.highlight = c('blue', 'green', 'red', 'yellow')))
dev.off()

png(paste0(wd, '/output/images/UMAP/umap_sub_cluster_condition_type.png'), width = 3500,
    height = 2500, res=400)
print(DimPlot(cluster_subset,
              reduction = "umap",
              label = F,
              label.size = 4, split.by = 'condition',
              group.by = 'type'))
dev.off()

png(paste0(wd, '/output/images/UMAP/umap_all_samples_no_label.png'), width = 3500,
    height = 2500, res=400)
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

png(paste0(wd, '/output/images/UMAP/umap_all_samples_model.png'), width = 3500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = F,
              label.size = 4, split.by = 'model'))
dev.off()


png(paste0(wd, '/output/images/UMAP/sample_orig_ident_all_samples.png'), width = 6500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 4,#group.by = 'type', 
        split.by = 'orig.ident'))
dev.off()


# ----- X -----
png(paste0(wd, '/output/images/UMAP/sample_model_all_samples_labelled.png'), width = 6500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = TRUE,
              label.size = 4,#group.by = 'type', 
              split.by = 'model'))
dev.off()

png(paste0(wd, '/output/images/UMAP/sample_type_all_samples.png'), width = 6500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              #label = TRUE,
              label.size = 4, group.by = 'type', 
              split.by = 'model'))
dev.off()


png(paste0(wd, '/output/images/UMAP/sample_condition_all_samples.png'), width = 6500,
    height = 2500, res=400)
print(DimPlot(seurat_integrated,
              reduction = "umap",
              label = TRUE,
              label.size = 4,
              #group.by = 'type', 
              split.by = 'condition'))
dev.off()

# No. of cells in each condition
table(seurat_integrated@meta.data$condition)

# No. of cells in each cluster
table(seurat_integrated$integrated_snn_res.0.6)

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

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



write.csv(seurat_integrated.markers, paste0(wd, 'output/all_markers_0.6_all_samples.csv'))

stats <- seurat_integrated %>% dplyr::select(integrated_snn_res.0.6, condition) %>% 
  dplyr::group_by(integrated_snn_res.0.6,condition) %>% dplyr::summarise(n())

colnames(stats)[3] <- 'count'

write.csv(stats, paste0(wd, 'output/cluster_0.6_condition_stats_all_samples.csv'))
saveRDS(seurat_integrated, paste0(wd,'/RDS/seurat_integrated_all_samples_clustered.Rds'))

png(paste0(wd,'/dotplot_rna_assay.png'), res=200, height = 1200, width = 1000)
print(DotPlot(seurat_integrated, features = c('Fgfr1', 'Mif'), assay = 'RNA'))
dev.off()

