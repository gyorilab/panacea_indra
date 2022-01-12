library(dplyr)
library(Seurat)
library(DESeq2)
library(stringr)
library(pheatmap)
library(enrichplot)
library(EnhancedVolcano)
library(SingleCellExperiment)

## Notes
# Collapse the samples into one group 
# Cluster res = 0.8 (25 clusters)
# Cluster res = 0.6 (17 clusters)
# For res(0.6) imp clusters -> 0,5,2,7,13,10,4 (Run analysis on these clusters)
# Comparisons -> make comparisons b/w each samples

wd <- ('~/gitHub/singlecell_model_analysis/')
setwd(wd)
source('./functions.R')
seurat_integrated <- readRDS('./RDS/seurat_integrated_all_samples_clustered.Rds')

#seurat_stim <- subset(seurat_integrated, type=='STIM')


# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat_integrated@assays$RNA@counts 

metadata <- seurat_integrated@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat_integrated@active.ident)

metadata$sample_id <- factor(str_replace(seurat_integrated$orig.ident, '_', ''))

#Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk


# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)
ns


# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select(-"cluster_id")
ei

# Perform QC if not already performed
dim(sce)

# Calculate quality control (QC) metrics
sce_metrics <- scuttle::perCellQCMetrics(sce)

# Get cells w/ few/many detected genes
sce$is_outlier <- scuttle::isOutlier(
  metric = sce_metrics$total,
  nmads = 2, type = "both", log = TRUE)

# Remove outlier cells
sce <- sce[, !sce$is_outlier]
dim(sce)

## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

dim(sce)


# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]


# Aggregate across cluster-sample groups
pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)


dim(pb)

pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) set_colnames(
    t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+"))
    )


# Explore the different components of list
str(pb)

# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)



# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()


# Get cluster IDs for each of the samples
samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "condition", "type")])


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, condition, type)

metadata        

# Generate vector of cluster IDs
clusters <- levels(as.factor(metadata$cluster_id))
clusters

comparisons <- c()
imp_clusters <- c('RM', 'Macs4', 'Dermal Macs')
condition_groups <- list('condition'=c('Incision_Zymo', 
                                       'Incision_UVB', 
                                       'Zymo_UVB',
                                       'Incision_Healthy',
                                       'Zymo_Saline',
                                       'UVB_Sham'),
                         'type'=c('HT_STIM','STIM_HT'))

dir.create('./output/pseudo_bulk_analysis', 
           recursive = T, showWarnings = F)
dir.create(paste0(wd,'/output/pseudo_bulk_analysis/plots/volcano'), 
           showWarnings = F)
dir.create(paste0(wd,'/output/pseudo_bulk_analysis/diff_files'), 
           showWarnings = F)

for(n in imp_clusters){
  cluster_index <- which(clusters == n)
  # Testing the DESEQ with first cluster group
  # Subset the metadata to only the B cells
  cluster_metadata <- metadata[which(
    metadata$cluster_id == clusters[cluster_index]), ]
  head(cluster_metadata)
  
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  head(cluster_metadata)
  
  # Subset the counts to only the B cells
  counts <- pb[[clusters[cluster_index]]]
  
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% 
                                                rownames(cluster_metadata))])
  
  for(condition_name in names(condition_groups)){
    for(comparisons in condition_groups[[condition_name]]){
      # needd to only substitute the comparisons counts
      comparisons <- str_split(comparisons, '_')[[1]]
      active_rows <- which(cluster_metadata[, c(condition_name)] %in% 
                             comparisons)
      active_counts <- cluster_counts[, active_rows]
      dds <- DESeqDataSetFromMatrix(active_counts, 
                                    colData = cluster_metadata[active_rows,], 
                                    design = formula(paste0('~',condition_name)))
      dds <- DESeq(dds)
      
      res <- results(dds, 
                     contrast = c(condition_name, comparisons[1], 
                                  comparisons[2]),
                     alpha = 0.05)
      write.csv(res,
                paste0(wd,'output/pseudo_bulk_analysis/diff_files/',n,'_',comparisons[[1]],
                       '_',comparisons[[2]],'.csv'))
      dir.create(paste0('./output/pseudo_bulk_analysis/plots/volcano'), showWarnings = F,
                 recursive = T)
      volcano_generate(res, titlename = paste0(n,' ',paste0(comparisons[[1]], '_', comparisons[[2]])), 
                       file_loc = paste0(wd,'output/pseudo_bulk_analysis/plots/volcano/',
                                         n,'_',paste0(comparisons[[1]], '_', comparisons[[2]])))
      # Transform counts for data visualization
      rld <- rlog(dds, blind=TRUE)
      
      # Extract the rlog matrix from the object and compute pairwise correlation values
      anno <- cluster_metadata[active_rows, ]
      rld_mat <- assay(rld)
      rld_cor <- cor(rld_mat)
      
      # Plot heatmap
      dir.create('./output/pseudo_bulk_analysis/plots/heatmaps/', showWarnings = F)
      png(paste0('./output/pseudo_bulk_analysis/plots/heatmaps/',
                 n, '_', comparisons[[1]], '_', comparisons[[2]],
                 '_corr.png'),
          width = 1000, height = 1000, res=200)
      print(pheatmap(rld_cor, annotation = anno[, c(condition_name), 
                                                drop=F]))
      dev.off()
      
      
      # Plot PCA
      dir.create('./output/pseudo_bulk_analysis/plots/pca/', showWarnings = F,
                 recursive = T)
      png(paste0('./output/pseudo_bulk_analysis/plots/pca/cluster_',n,'_',
                 comparisons[[1]], '_', comparisons[[2]],'.png'),
          res=150, width = 1500, height = 1200)
      print(DESeq2::plotPCA(rld, intgroup = condition_name))
      dev.off()

      #cluster_metadata$condition <- as.factor(cluster_metadata$condition[active_rows])
      
      significant_genes <- list()
      grouped_comp <- data.frame('gene'='',
                                 'padj'='',
                                 'group'='')
      
      # Turn the results object into a tibble for use with tidyverse functions
      res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
      
      # Significant genes
      # Set thresholds
      padj_cutoff <- 0.05
      
      # Subset the significant results
      sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
        dplyr::arrange(padj)
      up_reg <- dplyr::filter(res_tbl, log2FoldChange > 0) %>%
        dplyr::arrange(desc(log2FoldChange))
      
      dir.create('./output/pseudo_bulk_analysis/up-regulated',showWarnings = F)
      write.csv(up_reg,
                paste0('./output/pseudo_bulk_analysis/up-regulated/cluster',n, "_", 
                       comparisons[[1]], '_', comparisons[[2]],'_up_reg.csv'),
                quote = FALSE, 
                row.names = FALSE)
      
      if(nrow(sig_res) >= 1){
        # Check significant genes output
        #sig_res
        
        # Write significant results to file
        dir.create('./output/pseudo_bulk_analysis/sig_diff_files',showWarnings = F)
        write.csv(sig_res,
                  paste0('./output/pseudo_bulk_analysis/sig_diff_files/cluster',n, "_",
                         comparisons[[1]], '_', comparisons[[2]],'_sig_genes.csv'),
                  quote = FALSE, 
                  row.names = FALSE)
        
        comp <- paste(comparisons[[1]], comparisons[[2]], sep='_')
        significant_genes[[comp]] <- sig_res %>%
          dplyr::select('gene', 'padj') %>%
          dplyr::arrange(padj) %>% head(n=20)
        significant_genes[[comp]]$group <- comp
        grouped_comp <- rbind(grouped_comp, 
                              significant_genes[[comp]])
        
        ## ggplot of top genes
        normalized_counts <- counts(dds, 
                                    normalized = TRUE)
        
        
        ## Order results by padj values
        top20_sig_genes <- sig_res %>%
          dplyr::arrange(padj) %>%
          dplyr::pull(gene) %>%
          head(n=20)
        d <- as.data.frame(sapply(colnames(normalized_counts), 
                                  function(x) str_extract(
                                    colnames(normalized_counts), x)))
        
        cols <- as.numeric(rownames(d[rowSums(is.na(d)) != ncol(d), ]))
        
        top20_sig_norm <- data.frame(normalized_counts)[, cols] %>% 
          rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% top20_sig_genes)
        
        gathered_top20_sig <- top20_sig_norm %>%
          gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))],
                 key = "samplename", value = "normalized_counts")
      
        gathered_top20_sig <- inner_join(ei[, c("sample_id", "condition", "type" )], 
                                         gathered_top20_sig, by = c("sample_id" = "samplename"))
        
        ## plot using ggplot2
        dir.create('./output/pseudo_bulk_analysis/plots/scatter_plot/',
                   showWarnings = F)
        png(paste0('./output/pseudo_bulk_analysis/plots/scatter_plot/',
                   n,'top20_sig_de_',  comparisons[[1]], '_',
                   comparisons[[2]], '.png'),
            width = 1000, height = 1000, res=200)
        
        print(ggplot(gathered_top20_sig) +
                geom_point(aes(x = gene, 
                               y = normalized_counts, 
                               color = type), 
                           position=position_jitter(w=0.1,h=0)) +
                scale_y_log10() +
                xlab("Genes") +
                ylab("log10 Normalized Counts") +
                ggtitle("Top 20 Significant DE Genes") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                theme(plot.title = element_text(hjust = 0.5)))
        dev.off()
        
        # Extract normalized counts for only the significant genes
        sig_norm <- data.frame(normalized_counts)[, cols] %>%
          rownames_to_column(var = "gene") %>%
          dplyr::filter(gene %in% sig_res$gene)
        
        # Set a color palette
        heat_colors <- brewer.pal(6, "YlOrRd")
        
        # Run pheatmap using the metadata data frame for the annotation
        png(paste0('./output/pseudo_bulk_analysis/plots/heatmaps/',
                   n,'heatmaps_', comparisons[[1]], '_',
                   comparisons[[2]],'.png'),
            width = 1000, height = 1000, res=200)
        print(pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
                       color = heat_colors[2:length(colnames(sig_norm))], 
                       cluster_rows = T, 
                       show_rownames = F,
                       annotation = cluster_metadata[, c(condition_name, "cluster_id")][active_rows, ],
                       border_color = NA, 
                       fontsize = 10, 
                       scale = "row", 
                       fontsize_row = 10, 
                       height = 20
        ))
        dev.off()
        
      }
      else{
        print(paste0('Skipping, no significant genes for ', comparisons[[1]],
                     ' VS ',comparisons[[2]], ' in cluster ', n))
        }
      }
    }
}


