library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(dittoSeq)
library(tidyseurat)
library(ComplexHeatmap)


wd <- 'gitHub/panacea_indra/panacea_indra/single_cell_analysis/'
source(paste0(wd, '/scripts/api.R'))

dir.create(paste0(wd, '/output/images/ncells/'), 
           showWarnings = F, recursive = T)

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
  'AGCTAGAA'='Zymo_contra_1',
  'ACTCTAGG'='Zymo_contra_2',
  'TACTCCTT'='UVB_contra_1',
  'ATTAGACG'='UVB_contra_2',
  'AGGCTTAG'='UVB_1',
  'CGGAGAGA'='UVB_2',
  'ATAGAGAG'='Incision_contra_1',
  'TATGCAGT'='Incision_contra_2',
  'AGAGGATA'='Incision_1',
  'CTCCTTAC'='Incision_2'
)

# Give sample types
HT <- c('Incision_contra_1', 'Incision_contra_2', 'Zymo_contra_1', 
        'Zymo_contra_2', 'UVB_contra_1', 'UVB_contra_2')
STIM <- c('Zymo_1', 'Zymo_2', 'Incision_1', 'Incision_2', 'UVB_1', 'UVB_2')

models <- list(
  'Incision'='Incision',
  'Incision_contra'='Incision',
  'Zymo'='Zymo',
  'Zymo_contra'='Zymo',
  'UVB'='UVB',
  'UVB_contra'='UVB'
)

condition <- list(
  'Incision_1' = 'Incision',
  'Incision_2' = 'Incision',
  'Incision_contra_1' = 'Incision_contra',
  'Incision_contra_2' = 'Incision_contra',
  'Zymo_1' = 'Zymo',
  'Zymo_2' = 'Zymo',
  'Zymo_contra_1' = 'Zymo_contra',
  'Zymo_contra_2' = 'Zymo_contra',
  'UVB_1' = 'UVB',
  'UVB_2' = 'UVB',
  'UVB_contra_1' = 'UVB_contra',
  'UVB_contra_2' = 'UVB_contra'
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
  object_list[[s]] <- CreateSeuratObject(infiles[[s]], project = sample_names[s],
                                         min.cells = 3, min.features = 200)
  # Check for sample type
  type <- if(sample_names[s] %in% HT) 'Healthy' else 'Injured'
  object_list[[s]]$type <- type
  # add replicate in the metadata
  object_list[[s]]$replicate <- str_match(sample_names[6], '(.*_)(\\d+)')[,3]
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

plot_density(wd, incision, 'incision', 'nFeature_RNA', 300, 2800)
plot_density(wd, incision, 'incision', 'nCount_RNA', 500)
plot_density(wd, incision, 'incision', 'percent_mito', 300)
plot_vln(wd, incision, 'all_incision')
saveRDS(incision, paste0(wd, '/RDS/incision.Rds'))

# n cells before filtering
plot_ncells(wd,'incision_before_filter', incision)
incision <- subset(incision, subset = nFeature_RNA > 300 & 
                     nFeature_RNA < 2800 & nCount_RNA > 500 & percent_mito < 5)
plot_ncells(wd, 'incision_after_filter', incision)


# Merge UVB
uvb_samples <- c(5, 7, 8, 6)
uvb <- object_list[uvb_samples]
uvb <- merge(uvb[[1]], uvb[-1],
             all.cell.ids = sample_names[uvb_samples])

uvb <- PercentageFeatureSet(uvb, "^mt-", col.name = "percent_mito")

plot_density(wd, uvb, 'UVB', 'nFeature_RNA', 300, 2800)
plot_density(wd, uvb, 'UVB', 'nCount_RNA', 300)
plot_density(wd, uvb, 'UVB', 'percent_mito', 300)
plot_vln(wd, uvb, 'all_uvb')
saveRDS(uvb, paste0(wd, '/RDS/uvb.Rds'))

plot_ncells(wd, 'uvb_before_filter', uvb)
uvb <- subset(uvb, subset = nFeature_RNA > 300 &
                nFeature_RNA < 2800 & nCount_RNA > 300 & percent_mito < 5)
plot_ncells(wd, 'uvb_after_filter', uvb)


# Merge Zymo
zymo_sample <- c(10, 9 , 11, 12)
zymo <- object_list[zymo_sample]
zymo <- merge(zymo[[1]], zymo[-1],
              all.cell.ids = sample_names[zymo_sample])

zymo <- PercentageFeatureSet(zymo, "^mt-", col.name = "percent_mito")

plot_density(wd, zymo, 'ZYMO', 'nFeature_RNA', 300, 2800)
plot_density(wd, zymo, 'ZYMO', 'nCount_RNA', 500)
plot_density(wd, zymo, 'ZYMO', 'percent_mito', 300)
plot_vln(wd, zymo, 'all_zymo')
saveRDS(zymo, paste0(wd, '/RDS/zymo.Rds'))

plot_ncells(wd, 'zymo_before_filter', zymo)
zymo <- subset(zymo, subset = nFeature_RNA > 300 & 
                     nFeature_RNA < 2800 & nCount_RNA > 500 & percent_mito < 5)
plot_ncells(wd, 'zymo_after_filter', zymo)


# Merge all the samples
all_data <- merge(incision, c(uvb, zymo), 
                  add.cell.ids = c("incision", "uvb",
                                   "zymo"))
saveRDS(all_data, paste0(wd, '/RDS/all_samples_merged_object.Rds'))

# SCTransform
split_seurat <- SplitObject(all_data, split.by = "condition")

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = 
                                     c("percent_mito"))
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
        paste0(wd, '/RDS/all_samples_integrated_anchors.Rds'))



# Integrate across conditions (Run this below step
# on machine with atleast 50GB Mem)
#seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
#                                   normalization.method = "SCT")
seurat_integrated <- readRDS(paste0(wd, 
                                    '/RDS/all_samples_integrated.Rds'))

# Normalize the RNA Assay
DefaultAssay(seurat_integrated) <- 'RNA'
seurat_integrated <- NormalizeData(seurat_integrated)
seurat_integrated <- ScaleData(seurat_integrated)


# Run PCA
DefaultAssay(seurat_integrated) <- 'integrated'
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
tiff(paste0(wd, '/output/images/QC/pca_split_orig_ident_all_samples.tiff'), width = 5050,
    height=2000, res=200)
print(PCAPlot(seurat_integrated,
        split.by = "orig.ident"))
dev.off()

tiff(paste0(wd, '/output/images/QC/pca_all_samples.tiff'), width = 2050,
    height=1500, res=200)
print(PCAPlot(seurat_integrated))
dev.off()

tiff(paste0(wd, '/output/images/QC/elbow_plot.tiff'), width = 1050,
    height=1000, res=200)
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


# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = 0.6)

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"

# Find all the DE genes from all clusters
seurat_integrated.markers <- FindAllMarkers(seurat_integrated, 
                                            only.pos = TRUE,
                                            min.pct = 0.25, 
                                            logfc.threshold = 0.25)

write.csv(seurat_integrated.markers, paste0(wd, 'output/all_markers_0.6_all_samples.csv'))


saveRDS(seurat_integrated.markers, paste0(wd, 
                                          '/RDS/seurat_integrated_markers.Rds'))

# Annotating clusters
clusters.ids <- c("Neutrophils", "DC2", "Macs4", "Dermal Macs", "Macs3",
                  "T cells", "RM", "Macs1", "Mast cells", "DC1", "LCs", "Macs2",
                  "Keratinocytes", "Unknown", "Unknown", "Unknown", "Unknown",
                  "Unknown")

names(clusters.ids) <- levels(seurat_integrated)[c(1:length(clusters.ids))]
seurat_integrated <- RenameIdents(seurat_integrated, clusters.ids)
seurat_integrated$cell_type <- Idents(seurat_integrated)
saveRDS(seurat_integrated, paste0(wd, '/RDS/seurat_integrated_0.6.RDS'))

## End of cluster annotation


## Notes
#1. Remove Unknown clusters and remake the object
seurat_integrated_unknown_removed <- subset(seurat_integrated, idents='Unknown', invert=T)
saveRDS(seurat_integrated_unknown_removed, paste0(wd,'/RDS/seurat_integrated_unknown_removed.Rds'))
#2. Fix the dittoHeatmap color

