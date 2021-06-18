library(tm)
library(dplyr)
library(Seurat)
library(ggplot2)
library(proustr)
library(pheatmap)
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
# Read in Seurat objects
incision <- readRDS('./input/data/incision_all.rds')
inn <- NormalizeData(incision)
healthy <- readRDS('./input/data/healthy_all.rds')

# Rename cell clusters  
#immune_cells <- RenameIdents(immune_cells,  `2` = "Neutrophils",
#                             `7` = "Neutrophils", `0` = "Monocytes",
#                             `5` = "M1a", `1` = "M1b", 
#                             `3` = "M2a", `4` = "M2a",`9` = "M2b", 
#                             `6` = "Dendritic cells",   
#                             `8` = "Dendritic cells",  `10` = "Dendritic cells",
#                             `11` = "Mast cells", `12` = "ab T"
#                             )



# Read neuron RDS
neuron_cells <- readRDS('./input/data/nucseq_interactome.rds')
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
# Combining immune incision cells and neuron cells
neuro_immune.combined <- merge(incision, y = neuron_cells, 
                               add.cell.ids = c("incision", "neurons"), 
                               project = "incisio_neuron")
neuro_immune.combined <- merge(healthy, neuron_cells,
                               add.cell.ids = c('healthy', 'neurons'),
                               project = 'healthy_neuron')

# Write meta table of cells and cell type
# Separated by tab-space
meta_table <- data.frame('Cell' = names(neuro_immune.combined@active.ident),
                         'cell_type' = unname(neuro_immune.combined@active.ident))
write.table(meta_table, './input/neuro_healthy_meta_table.tsv', 
            row.names = F, sep='\t', quote = F)


# Make a meta table for immune cells
#immune_meta_table <- meta_table[which(meta_table$Cell ==
#                                        str_extract(meta_table$Cell, 'incision.*')), ]
#write.table(immune_meta_table, './input/immune_meta_table.tsv', row.names = F, 
#          sep = '\t', quote=F)

# Read BioMart Mouse Human Orthologue table
mart_table <- read.table('./input/mart_export.txt', sep='\t', header = T)
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

# Filtering the genes in neuro immune counts to HGNC eq using INDRA
# modules
#indra_hgnc_genes <- unlist(fastGene::convert_symbols(rownames(neuro_immune_counts)))

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
                                            hgnc_mart[duplicated(hgnc_mart$HGNC_SYMBOL), c(1)], ])
  )

# Saving duplicated Genes to check for Ligands/Receptors
#duplicated_genes <- unique(hgnc_mart[duplicated(hgnc_mart$HGNC_SYMBOL), 1])

# Read ligands and receptors list
#ligands <- read.csv('./data/ligands.txt', col.names = 'genes')
#receptors <- read.csv('./data/receptors.txt', col.names = 'genes')

#table(duplicated_genes %in% ligands | duplicated_genes %in% receptors)


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
write.table(hgnc_mart_ensmbl, './input/neuro_healthy_counts.txt',
            col.names = T, row.names = F, quote = F, sep='\t')

# Save only Immune counts in cellponeDB format
#immune_hgnc_mart_ensmbl <- 
#  hgnc_mart_ensmbl[, c(1, which(colnames(hgnc_mart_ensmbl) == 
#                             str_extract(colnames(hgnc_mart_ensmbl), 'immune.*')))]
#write.table(immune_hgnc_mart_ensmbl, './immune_counts.txt',
#            col.names = T, row.names = F, quote = F, sep='\t')





#### END OF CREATING META DATA ####







# subsetting only M1a cluster
m1a_cluster <- immune_hgnc_mart_ensmbl[, which(colnames(immune_hgnc_mart_ensmbl)
                                               %in% subset(immune_meta_table, cell_type == 'M1a')[,1] == T)]

# subsetting only Monocyte cluster
monocyte_cluster <- immune_hgnc_mart_ensmbl[, which(colnames(immune_hgnc_mart_ensmbl)
                                               %in% subset(immune_meta_table, cell_type == 'Monocytes')[,1] == T)]

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










###########
# Code to make custom interactions file 
# for drawing heatmaps
###########
meta = read.csv(file = './input/immune_meta_table.tsv', comment.char = '', sep='\t')
all_intr = read.table('./out/celldb_user_custom/pvalues.txt', header=T, stringsAsFactors = F, 
                      sep='\t', comment.char = '', check.names = F)
all_intr <- subset(all_intr, receptor_a == 'False' & receptor_b=='True')

intr_pairs = all_intr$interacting_pair
all_intr = all_intr[,-c(1:11)]
split_sep = '\\|'
join_sep = '|'

pairs1_all = unique(meta[,2])


# Creating pairs
pairs1 = c()
for (i in 1:length(pairs1_all))
  for (j in 1:length(pairs1_all))
    pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))

# Get only interactions to abT
abT_intr <- na.omit(str_extract(pairs1, '.*\\|ab T'))

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


