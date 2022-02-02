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


## functions
symbol2ensembl <-function(cts){
  ens <- rownames(cts)
  ens <- as.character(ens)
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('ENSEMBL'), keytype = 'SYMBOL')
  #symbols
  
  rownames(cts) <-symbols
  keep <- !is.na(rownames(cts))
  cts <- cts[keep,]
  return(cts)
}


# Set working directory
setwd('~/gitHub/panacea_indra/panacea_indra/cellphonedb_processor/')
seurat_integrated <- readRDS('../single_cell_analysis/RDS/seurat_integrated_annotated_0.6.RDS')
incision <- subset(seurat_integrated, condition=='Incision')
zymosan <- subset(seurat_integrated, condition=='Zymosan')
uvb <- subset(seurat_integrated, condition=='UVB')
incision_contra <- subset(seurat_integrated, condition=='Incision_contra')

immune_object <- list('incision' = incision,
                      'zymosan' = zymosan,
                      'uvb' = uvb,
                      'incision_contra' = incision_contra)



# Remove Undetermined cluster
immune_object$incision <- subset(immune_object$incision, cell_type != 'Undetermined')
immune_object$zymosan <- subset(immune_object$zymosan, cell_type != 'Undetermined')
immune_object$uvb <- subset(immune_object$uvb, cell_type != 'Undetermined')
immune_object$incision_contra <- subset(immune_object$incision_contra, 
                                        cell_type != 'Undetermined')


if(file.exists('../single_cell_analysis/RDS/neuron_cells.Rds') != T){
  
  # Read neuron data
  neuron_cells <- readRDS('../single_cell_analysis/input_files/GSE154659_C57_Raw_counts.Rds')
  neuron_cells <- CreateSeuratObject(neuron_cells)
  neuron_cells@meta.data$cells <- rownames(neuron_cells@meta.data)
  neuron_cells@meta.data$sample <- NA
  neuron_cells$sample[which(str_detect(neuron_cells$cells, "male_C57_Naive_0_"))] <- "male_C57_Naive"
  neuron_cells$sample[which(str_detect(neuron_cells$cells, "female_C57_Naive_0_"))] <- "female_C57_Naive"
  neuron_cells <- subset(neuron_cells,
                         subset = (sample == "male_C57_Naive" |
                                   sample== "female_C57_Naive"))
  
  # Normalize and scale data for the two samples (male and female)
  neuron_cells_split <- SplitObject(neuron_cells, split.by = "sample")
  
  for (i in 1:length(neuron_cells_split)) {
    neuron_cells_split[[i]] <- NormalizeData(neuron_cells_split[[i]])
    neuron_cells_split[[i]] <- ScaleData(neuron_cells_split[[i]])
  }
  neuron_cells <- merge(neuron_cells_split[[1]], 
                        neuron_cells_split[[2]], 
                        merged.data=T)
  DefaultAssay(neuron_cells) <- 'RNA'
  
  # Extract the neuron cell types from the meta data
  neuron_cell_types <- str_match(neuron_cells@meta.data$cells, 
                                 "(.*rep\\d_)([A-Za-z0-9\\s^_]+)(_.*)")[,3]
  neuron_cells@meta.data$cellID <- neuron_cell_types
  
  
  # Subset the object by celltypes
  neuron_cells <- subset(neuron_cells,
                         subset = (
                           cellID == "cLTMR1" |  cellID == "p_cLTMR2" | 
                           cellID == "PEP1" |  cellID == "PEP2" | 
                           cellID == "NP" | cellID == "SST" |
                           cellID == "NF1" | cellID == "NF2" | 
                           cellID =="NF3"))
  
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
  
  saveRDS(neuron_cells, '../single_cell_analysis/RDS/neuron_cells.Rds')
}else{
  neuron_cells <- readRDS('../single_cell_analysis/RDS/neuron_cells.Rds')
}

dir.create('./cpdb_input/meta_tables', showWarnings = F)
dir.create('./cpdb_input/counts', showWarnings = F)
count = 0
## Runnning loop to merge and save neuro and immune objects
for (immune in immune_object) {
  count = count+1
  xname <- names(immune_object)[count]
  # Combining immune incision cells and neuron cells
  DefaultAssay(immune) <- 'RNA'
  neuro_immune.combined <- merge(immune, y = neuron_cells, 
                                 add.cell.ids = c(xname, "neurons"), 
                                 project = paste0(xname,"_neuron"))
  

  # Write meta table of cells and cell type
  # Separated by tab-space
  meta_table <- data.frame('Cell' = names(neuro_immune.combined@active.ident),
                           'cell_type' = unname(neuro_immune.combined@active.ident))
  write.table(meta_table, paste0('./cpdb_input/meta_tables/neuro_',xname,'_meta_table.tsv'), 
              row.names = F, sep='\t', quote = F)
  

  # Read BioMart Mouse Human Orthologue table
  mart_table <- read.table('~/.biomart/human_mouse_ortho.tsv', sep='\t', header = T)
  # Removing rows with no value
  mart_table <- mart_table[Reduce(`&`, lapply(mart_table, function(x) !x=="")),]

  # Extracting counts from neuro immune object
  neuro_immune_counts <- as.data.frame(neuro_immune.combined@assays$RNA@data)

  neuro_immune_filtered <- neuro_immune_counts[rownames(neuro_immune_counts) %in% 
                                               mart_table$MOUSE_SYMBOL, ]
  neuro_immune_filtered$MOUSE_SYMBOL <- rownames(neuro_immune_filtered)

  neuro_immune_hgnc <- merge(neuro_immune_filtered, mart_table, 
                             by='MOUSE_SYMBOL')
  neuro_immune_hgnc <- neuro_immune_hgnc[, c(ncol(neuro_immune_hgnc),
                                             1:ncol(neuro_immune_hgnc)-1)]

  neuro_immune_hgnc <- neuro_immune_hgnc[!duplicated(neuro_immune_hgnc$HUMAN_SYMBOL), ]

  rownames(neuro_immune_hgnc) <- neuro_immune_hgnc$HUMAN_SYMBOL
  cnum <- which(colnames(neuro_immune_hgnc) == 'HUMAN_SYMBOL')
  colnames(neuro_immune_hgnc)[cnum] <- 'Gene'
  neuro_immune_hgnc$MOUSE_SYMBOL <- NULL

  # Convert HGNC symbol to ENSMBL
  ensmbl <- symbol2ensembl(as.matrix(neuro_immune_hgnc))
  ensmbl[,1] <- rownames(ensmbl)

  # Save counts in cellphoneDB format
  write.table(ensmbl, paste0('./cpdb_input/counts/neuro_',xname,'_counts.txt'),
              col.names = T, row.names = F, quote = F, sep='\t')
}

################################################################################
