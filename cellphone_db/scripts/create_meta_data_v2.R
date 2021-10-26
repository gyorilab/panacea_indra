library(tm)
library(dplyr)
library(Seurat)
library(ggplot2)
library(proustr)
library(stringr)
library(pheatmap)
library(tidyverse)
library(hrbrthemes)
library(VennDiagram)
library(RColorBrewer)
library(AnnotationDbi)
library('org.Hs.eg.db')
library('org.Mm.eg.db')

# Set working directory
setwd('~/gitHub/panacea_indra/cellphone_db/')

# Sourcing functions
source('./scripts/functions.R')

###NOTES

# Normalization -> Immune cells
#               !-> Neuron cells  
# Need to check for lignads/receptors in the multi mapped genes

####

## functions


clusters.ids <- list('0'="Neutrophil", '2'="RM2", '3'="Dermal Macs", '6'="RM1")



# Read in Seurat objects
seurat_integrated <- readRDS('~/gitHub/singlecell_model_analysis/RDS/seurat_integrated_all_samples_clustered.Rds')
seurat_integrated <- subset(seurat_integrated, subset=cell_type %in% c(13,14,15,16),
                            invert = T)

# Make a list of all the conditions
immune_object <- list(
  'Healthy' = subset(seurat_integrated, condition == 'Healthy'),
  'Incision' = subset(seurat_integrated, condition == 'Incision'),
  'Saline' = subset(seurat_integrated, condition == 'Saline'),
  'Sham' = subset(seurat_integrated, condition == 'Sham'),
  'UVB' = subset(seurat_integrated, condition == 'UVB'),
  'Zymo' = subset(seurat_integrated, condition == 'Zymo')
)

# Temp fix for uvb trailing new line in cell names
immune_object[['UVB']] <- RenameCells(immune_object[['UVB']], new.names = 
                                        str_replace(colnames(immune_object[['UVB']]), 
                                                    '\\s+ ', ""))
immune_object[['Sham']] <- RenameCells(immune_object[['Sham']], new.names = 
                                        str_replace(colnames(immune_object[['Sham']]), 
                                                    '\\s+ ', ""))

# Read neuron RDS
neuron_cells <- readRDS('./input/data/GSE154659_C57_Raw_counts.Rds')
neuron_cells <- CreateSeuratObject(neuron_cells)
neuron_cells@meta.data$cells <- rownames(neuron_cells@meta.data)
neuron_cells@meta.data$sample <- NA
neuron_cells$sample[which(str_detect(neuron_cells$cells, "male_C57_Naive_0_"))] <- "male_C57_Naive"
neuron_cells$sample[which(str_detect(neuron_cells$cells, "female_C57_Naive_0_"))] <- "female_C57_Naive"
neuron_cells <- subset(neuron_cells,
                       subset = (sample == "male_C57_Naive" |  sample== "female_C57_Naive"))

# SCTransform
neuron_cells_split <- SplitObject(neuron_cells, split.by = "sample")

for (i in 1:length(neuron_cells_split)) {
  neuron_cells_split[[i]] <- SCTransform(neuron_cells_split[[i]])
}
neuron_cells <- merge(neuron_cells_split[[1]], neuron_cells_split[[2]], merged.data=T)

DefaultAssay(neuron_cells) <- 'SCT'

# Extract the neuron cell types from the meta data
neuron_cell_types <- str_match(neuron_cells@meta.data$cells, 
                               "(.*rep\\d_)([A-Za-z0-9\\s^_]+)(_.*)")[,3]
neuron_cells@meta.data$cellID <- neuron_cell_types


# Subset the object by celltypes
neuron_cells <-     subset(neuron_cells,
                                subset = (
                                  cellID == "cLTMR1" |  cellID == "p_cLTMR2" | 
                                    cellID == "PEP1" |  cellID == "PEP2" | 
                                    cellID == "NP" | cellID == "SST" |
                                    cellID == "NF1" | cellID == "NF2" | 
                                    cellID =="NF3")
                                )

## Below chunk of code is to rename cell idents in neuron data
active.ident_neuron <- data.frame('cells' = names(neuron_cells@active.ident),
                                  'gender' = unname(neuron_cells@active.ident))
neuron_cell_types <- data.frame('cells' = neuron_cells$cells,
                                'cell_types' = neuron_cells$cellID)

active.ident_neuron <- merge(active.ident_neuron, neuron_cell_types, by='cells')
rownames(active.ident_neuron) <- active.ident_neuron$cells
active.ident_neuron$gender <- NULL
renamed_idents <- factor(active.ident_neuron$cell_types)
names(renamed_idents) <- active.ident_neuron$cells
neuron_cells@active.ident <- renamed_idents

# Make dotplots
genes_to_plot <- c('Ptger1', 'Il1rap', 'Tnfrsf1a', 'Osmr', 
                   'Ccr1', 'Itga9', 'Cd47', 'Igf2r')

png('../cellphone_db/output/plots/dotplot_1.png', width=1500,
    height = 1000, res = 150)
print(DotPlot(neuron_cells, features = genes_to_plot))
dev.off()
saveRDS(neuron_cells, './output/pkl/neuron_cells.Rds')

dir.create('output/counts', showWarnings = F)
dir.create('output/meta_data', showWarnings = F)
count = 0
## Runnning loop to merge and save neuro and immune objects
for (immune in immune_object) {
  count = count+1
  xname <- names(immune_object)[count]
  
  print(paste0("Processing ",xname, " sample"))
  
  # Extract meta data 
  immune_meta <- data.frame('Cell' = names(immune@active.ident),
                            'cell_type' = unname(immune@active.ident))
  
  neuron_meta <- data.frame('Cell' = names(neuron_cells@active.ident),
                            'cell_type' = unname(neuron_cells@active.ident))
  
  neuro_immune_meta <- rbind(immune_meta, neuron_meta)
  
  # Write meta table of cells and cell type
  # Separated by tab-space
  write.table(neuro_immune_meta, paste0('./output/meta_data/neuro_',xname,'_meta_table.tsv'), 
              row.names = F, sep='\t', quote = F)
  
  
  # Combining immune incision cells and neuron cells
  DefaultAssay(immune) <- 'SCT'
  neuro_immune.combined <- merge(immune, y = neuron_cells, 
                                 #add.cell.ids = c(xname, "neurons"), 
                                 #project = paste0(xname,"_neurons"), 
                                 merge.data=T)
  
  # Read BioMart Mouse Human Orthologue table
  mart_table <- read.table('./input/mart_export.txt', sep='\t', header = T)
  colnames(mart_table) <- c('HGNC_ID', 'HGNC_SYMBOL', 'MOUSE_SYMBOL', 'MOUSE_ID')
  
  # Removing rows with no value
  mart_table <- mart_table[Reduce(`&`, lapply(mart_table, function(x) !x=="")),]

  # Extracting counts from neuro immune object
  neuro_immune_counts <- as.data.frame(neuro_immune.combined@assays$SCT@data)

  # Filtering the genes in neuro immune counts to the genes in BioMart DB
  neuro_immune_filtered <- neuro_immune_counts[rownames(neuro_immune_counts) %in% 
                                               mart_table$MOUSE_SYMBOL, ]
  
  neuro_immune_filtered$MOUSE_SYMBOL <- rownames(neuro_immune_filtered)
  
  neuro_immune_filtered <- merge(neuro_immune_filtered, mart_table[,c(2,3)], 
                     by='MOUSE_SYMBOL')
  neuro_immune_filtered <- neuro_immune_filtered[, c(ncol(neuro_immune_filtered),1:ncol(neuro_immune_filtered)-1)]
  

  # Removing duplicated genes and only keeping
  # one -> one mapped genes
  neuro_immune_filtered <- neuro_immune_filtered[!neuro_immune_filtered$HGNC_SYMBOL %in% 
                                    neuro_immune_filtered[duplicated(neuro_immune_filtered$HGNC_SYMBOL), c(1)], ]
  
  rownames(neuro_immune_filtered) <- neuro_immune_filtered$HGNC_SYMBOL
  cnum <- which(colnames(neuro_immune_filtered) == 'HGNC_SYMBOL')
  colnames(neuro_immune_filtered)[cnum] <- 'Gene'
  
  neuro_immune_filtered$MOUSE_SYMBOL <- NULL
  
  # Convert HGNC symbol to ENSMBL
  neuro_immune_filtered_ensmbl <- symbol2ensembl(as.matrix(neuro_immune_filtered))
  neuro_immune_filtered_ensmbl[,1] <- rownames(neuro_immune_filtered_ensmbl)
  
  # Save counts in cellphoneDB format
  write.table(neuro_immune_filtered_ensmbl, paste0('./output/counts/neuron_',xname,'_counts.txt'),
              col.names = T, row.names = F, quote = F, sep='\t')

  
}



#### END OF CREATING META DATA and EXTRACTING COUNTS ####







# subsetting only M1a cluster
m1a_cluster <- immune_neuro_immune_filtered_ensmbl[, which(colnames(immune_neuro_immune_filtered_ensmbl)
                                               %in% subset(immune_meta_table, cell_type == 'M1a')[,1] == T)]

# subsetting only Monocyte cluster
monocyte_cluster <- immune_neuro_immune_filtered_ensmbl[, which(colnames(immune_neuro_immune_filtered_ensmbl)
                                               %in% subset(immune_meta_table, cell_type == 'Monocytes')[,1] == T)]











# Read significant means and get custom interactions
sig_means <- read.table('./out/only_immune_customDB/significant_means.txt',
                        sep='\t')
colnames(sig_means) <- sig_means[1,]
sig_means <- sig_means[-1,] 

# Change the cell types manually to get the desired
# list of interactions
t_cols <- which(colnames(sig_means) %in% str_extract(colnames(sig_means), 
                                                     pattern = ".*\\|ab T"))
#m2a_cols <- which(colnames(sig_means) %in% str_extract(colnames(sig_means), 
#                                                       pattern = "Neutrophils.*ab T"))
all_count = matrix(ncol=4)
colnames(all_count) = c('SOURCE','TARGET','Count', 'Pairs')


for (n in 1:length(t_cols)) {
  c_abT <- sig_means[,c(1:12, t_cols[n])]
  c_abT <- c_abT %>% filter(c_abT[, ncol(c_abT)] != "" & receptor_a=='False' & 
                              receptor_b=='True')
  pairs <- colnames(sig_means)[t_cols[n]]
  pairs <- strsplit(pairs, '\\|')[[1]]
  p1 <- pairs[[1]]
  p2 <- pairs[[2]]
  intr_pairs <- paste(c_abT$interacting_pair, collapse=', ')
  new_count = c(p1, p2, length(unique(c_abT$interacting_pair)), intr_pairs)
  names(new_count) = c('SOURCE','TARGET','Count', 'Pairs')
  all_count = rbind(all_count, new_count)
}

all_count <- all_count[-1, ]
write.table(all_count, './abT_interactions.tsv', sep='\t', row.names = F)
#%>%
#  write.csv('./out/only_immune_customDB/T_M2a_interactions.csv')














# Get only interactions to abT
abT_intr <- na.omit(str_extract(pairs1, '.*\\|ab T cells'))

all_count = matrix(ncol=4)
colnames(all_count) = c('SOURCE','TARGET','Count', 'Pairs')
count1 = c()
pvalue=0.05
for(i in 1:length(abT_intr)){
  p1 = strsplit(abT_intr[i], split_sep)[[1]][1]
  p2 = strsplit(abT_intr[i], split_sep)[[1]][2]
  
  n1 = intr_pairs[which(all_intr[,abT_intr[i]]<=pvalue)]

  #pairs_rev = paste(p2, p1, sep=join_sep)
  #n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
  
  #if(p1!=p2)
  #  count1 = length(unique(n1))+length(unique(n2))
  #else
  #  count1 = length(unique(n1))
  
  
  pairs <- paste(n1, collapse=' ,')
  
  new_count = c(p1,p2,length(unique(n1)), pairs)
 
  names(new_count) = c('SOURCE','TARGET','Count', 'Pairs')
  all_count = rbind(all_count, new_count)
}

all_count <- all_count[-1, ]

write.table(all_count[-4, ], './out/celldb_user_custom/abT_interactions.tsv', 
            sep='\t', row.names = F)
mat <- matrix(as.numeric(all_count[,3]), nrow = nrow(all_count), ncol = 1)
rownames(mat) <-  all_count[,1]
colnames(mat) <- 'ab T'
mat <- mat[-4,]
png('./example_plot.png')
print(pheatmap(t(mat), cluster_rows = F, cluster_cols = F, show_rownames = T))
dev.off()


