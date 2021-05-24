library(tm)
library(dplyr)
library(Seurat)
library(ggplot2)
library(proustr)
library(tidyverse)
library(hrbrthemes)
library(VennDiagram)
library(RColorBrewer)
library(AnnotationDbi)
library('org.Hs.eg.db')


###NOTES

# Normalization -> Immune cells
#               !-> Neuron cells  
# Need to check for lignads/receptors in the multi mapped genes

####

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
setwd('~/gitHub/panacea_indra/cellphone_db/')
# Read in Immune Cells RDS
immune_cells <- readRDS('./data/rds_IC_cluster.rds')
# Rename cell clusters
immune_cells <- RenameIdents(immune_cells,  `2` = "Neutrophils",
                             `7` = "Neutrophils", `0` = "Monocytes",
                             `5` = "M1a", `1` = "M1b", 
                             `3` = "M2a", `4` = "M2a",`9` = "M2b", 
                             `6` = "Dendritic cells",   
                             `8` = "Dendritic cells",  `10` = "Dendritic cells",
                             `11` = "Mast cells", `12` = "ab T")

# Read neuron RDS
neuron_cells <- readRDS('./data/nucseq_interactome.rds')
# Remove samples with NA
neuron_cells <- subset(neuron_cells, sample!='NA')

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
##

# Combining immune cells and neuron cells
neuro_immune.combined <- merge(immune_cells, y = neuron_cells, 
                       add.cell.ids = c("immune_cells", "neurons"), 
                       project = "Immune_Neuron")

# Write meta table of cells and cell type
# Separated by tab-space
meta_table <- data.frame('Cell' = names(neuro_immune.combined@active.ident),
                         'cell_type' = unname(neuro_immune.combined@active.ident))

write.table(meta_table, './meta_table.tsv', row.names = F, sep='\t', quote = F)


# Read BioMart Mouse Human Orthologue table
mart_table <- read.table('./data/mart_export.txt', sep='\t', header = T)
colnames(mart_table) <- c('HGNC_ID', 'HGNC_SYMBOL', 'MOUSE_SYMBOL', 'MOUSE_ID')
# Removing rows with no value
mart_table <- mart_table[Reduce(`&`, lapply(mart_table, function(x) !x=="")),]

# Extracting counts from neuro immune object
neuro_immune_counts <- as.data.frame(neuro_immune.combined@assays$RNA@counts)

# Converting counts to float values as per to 
# cellphoneDB format
for (coln in 1:ncol(neuro_immune_counts)) {
  neuro_immune_counts[,coln] <- noquote(sprintf("%.2f", 
                                                  neuro_immune_counts[,coln]))
} 


# Filtering the genes in neuro immune counts to the genes in BioMart DB
neuro_immune_filtered <- neuro_immune_counts[rownames(neuro_immune_counts) %in% 
                                       mart_table$MOUSE_SYMBOL, ]
neuro_immune_filtered$MOUSE_SYMBOL <- rownames(neuro_immune_filtered)

hgnc_mart <- merge(neuro_immune_filtered, mart_table[,c(2,3)], 
                   by='MOUSE_SYMBOL')
hgnc_mart <- hgnc_mart[, c(ncol(hgnc_mart),1:ncol(hgnc_mart)-1)]

cat('BioMart:', '\n',
  'genes mapped to HGNC ->',nrow(hgnc_mart), '\n',
  'one to one mapping ->', nrow(hgnc_mart[!hgnc_mart$HGNC_SYMBOL %in% 
                                              hgnc_mart[duplicated(hgnc_mart$HGNC_SYMBOL), c(1)], ]))

# Saving duplicated Genes to check for Ligands/Receptors
duplicated_genes <- unique(hgnc_mart[duplicated(hgnc_mart$HGNC_SYMBOL), 1])

# Read ligands and receptors list
ligands <- read.csv('./data/ligands.txt', col.names = 'genes')
receptors <- read.csv('./data/receptors.txt', col.names = 'genes')

table(duplicated_genes %in% ligands | duplicated_genes %in% receptors)


# Removing duplicated genes and only keeping
# one -> one mapped genes
hgnc_mart_filtered <- hgnc_mart[!hgnc_mart$HGNC_SYMBOL %in% 
                                hgnc_mart[duplicated(hgnc_mart$HGNC_SYMBOL), c(1)], ]
rownames(hgnc_mart_filtered) <- hgnc_mart_filtered$HGNC_SYMBOL
cnum <- which(colnames(hgnc_mart_filtered) == 'HGNC_SYMBOL')
colnames(hgnc_mart_filtered)[cnum] <- 'Gene'

hgnc_mart_filtered$MOUSE_SYMBOL <- NULL
# Convert HGNC symbol to ENSMBL
hgnc_mart_ensmbl <- symbol2ensembl(as.matrix(hgnc_mart_filtered))
hgnc_mart_ensmbl[,1] <- rownames(hgnc_mart_ensmbl)

# Save counts in cellphoneDB format
write.table(hgnc_mart_ensmbl, './counts.txt',
            col.names = T, row.names = F, quote = F, sep='\t')



## ggplots for stat check
stats_df <- data.frame('database' = c('HGNC table',
                                      'HGNC table',
                                      'BioMart',
                                      'BioMart',
                                      'informatics.jax',
                                      'informatics.jax'),
                       'mapping' = c('one->many',
                                     'one->one',
                                     'one->many',
                                     'one->one',
                                     'one->many',
                                     'one->one'),
                       'total_genes' = c(12071,
                                         12014,
                                         14192,
                                         12943,
                                         13133,
                                         12841)
)
p<-ggplot(stats_df, aes(x=database, y=total_genes, fill=mapping)) +
  geom_bar(stat="identity", position=position_dodge())+theme_minimal()+
  geom_text(aes(label=total_genes),position = position_dodge(0.9),
            vjust=1.6)

png('./compare_mapping.png', width = 700, height = 600, res = 100)
print(p)
dev.off()







## Venn diagram to compare HGNC genes
# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(2, "Pastel2")

venn.diagram(
  x = list(indra_hgnc_symbols,
           rownames(ic_hgnc_filtered)),
  category.names = c('INDRA', 'MGI_DB'),
  filename = paste0('./compare_hgnc.png'),
  output=T,
  
  # Output features
  imagetype="tiff" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 250,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c('lightgreen', 'lightblue'),
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(0, 1),
  #cat.dist = c(0.055, 0.010),
  cat.fontfamily = "sans",
  #rotation = 1
)
