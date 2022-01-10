library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(dittoSeq)
library(tidyseurat)
library(ComplexHeatmap)


setwd('~/gitHub/panacea_indra/cellphone_db/')

# Sourcing functions
source('./scripts/functions.R')
# read human mouse ortho table
mouse_human <- read.csv('./input/mouse_human_ortho.tsv', sep='\t')
colnames(mouse_human) <- c('mouse_id', 'mouse_gene_name', 'human_id', 'human_gene_name')

# Read receptors list
aj_receptors <- readxl::read_xlsx('./input/receptorsonly.xlsx', sheet = 1)
aj_receptors_genes <- aj_receptors$geneSymbol

# Read cellphone db receptors
cpdb_receptors <- read.csv('./input/cpdb_receptors.csv')
colnames(cpdb_receptors)[1] <- 'human_gene_name'
cpdb_receptors_merged <- merge(cpdb_receptors, mouse_human, 
                               by.y='human_gene_name')
cpdb_receptors_merged <- cpdb_receptors_merged[!duplicated(cpdb_receptors_merged$mouse_gene_name), ]


# all receptors
receptor_genes <- c(cpdb_receptors_merged$mouse_gene_name, aj_receptors_genes)

# Read neuron data
neuron_cells <- readRDS('./input/data/GSE154659_C57_Raw_counts.Rds')
neuron_cells <- CreateSeuratObject(neuron_cells)
neuron_cells@meta.data$cells <- rownames(neuron_cells@meta.data)
neuron_cells@meta.data$sample <- NA
neuron_cells$sample[which(str_detect(neuron_cells$cells, "male_C57_Naive_0_"))] <- "male_C57_Naive"
neuron_cells$sample[which(str_detect(neuron_cells$cells, "female_C57_Naive_0_"))] <- "female_C57_Naive"
neuron_cells <- subset(neuron_cells,
                       subset = (sample == "male_C57_Naive" |  
                                 sample== "female_C57_Naive"))

neuron_cell_types <- str_match(neuron_cells@meta.data$cells, 
                               "(.*rep\\d_)([A-Za-z0-9\\s^_]+)(_.*)")[,3]
neuron_cells@meta.data$cellID <- neuron_cell_types
neuron_cells_split <- SplitObject(neuron_cells, split.by = "sample")

female_subset <- neuron_cells_split$female_C57_Naive
male_subset <- neuron_cells_split$male_C57_Naive
# filter the data to neuron subsets
female_subset <- get_subsets(female_subset)
male_subset <- get_subsets(male_subset)

female_subset <- NormalizeData(female_subset)
female_subset <- ScaleData(female_subset, features = rownames(female_subset))
male_subset <- NormalizeData(male_subset)
male_subset <- ScaleData(male_subset, features = rownames(male_subset))


png('./output/plots/Il4ra_female.png', res=200, height = 1500, width=1000)
p1 <- VlnPlot(female_subset, features ='Il4ra', log = T)
print(p1)
dev.off()

png('./output/plots/Il4ra_male.png', res=200, height = 1500, width=1000)
p1 <- VlnPlot(male_subset, features ='Il4ra', log = T)
print(p1)
dev.off()

png('./output/plots/Cd44_Itgav_female.png', res=200, height = 1300, width=1200)
p1 <- VlnPlot(female_subset, features =c('Cd44', 'Itgav'), log = T)
print(p1)
dev.off()

png('./output/plots/Cd44_Itgav_male.png', res=200, height = 1300, width=1200)
p1 <- VlnPlot(male_subset, features =c('Cd44', 'Itgav'), log = T)
print(p1)
dev.off()


# get the differentially expressed receptors from male and female subsets
rg <- receptor_genes[receptor_genes %in% rownames(female_subset)==T]
rg <- rg[!duplicated(rg)]

female_np_pep1_cltmr1_sub <- subset(female_subset, subset=(cellID == 'NP' | cellID == 'PEP1' | cellID == 'cLTMR1'))
female_markers <- FindAllMarkers(female_np_pep1_cltmr1_sub, features = rg)
female_markers <- female_markers %>% group_by(cluster)

female_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05 & (avg_log2FC > 0.5 | avg_log2FC < -0.5)) %>%
  top_n(n = 40, wt = avg_log2FC) -> female_top40

png('./output/plots/heatmap_female_sig_genes_cLTMR1_NP_PEP1.png', res=150, height=2500,
    width=2300)
p <- DoHeatmap(female_np_pep1_cltmr1_sub, features = female_top40$gene, slot = 'scale.data')
print(p)
dev.off()

write.csv(female_markers, './output/de_files/neuron_female_np_pep1_cltmr1_all_receptors.csv')


#######
male_np_pep1_cltmr1_sub <- subset(male_subset, 
                                  subset=(cellID == 'NP' | cellID == 'PEP1' 
                                          | cellID == 'cLTMR1'))

male_markers <- FindAllMarkers(male_np_pep1_cltmr1_sub, features = rg)
male_markers_padj_0.05 <- male_markers %>% filter(p_val_adj < 0.05)

male_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05 & (avg_log2FC > 0.5 | avg_log2FC < -0.5)) %>%
  top_n(n = 40, wt = avg_log2FC) -> top40

png('./output/plots/heatmap_male_sig_genes_cLTMR1_NP_PEP1.png', res=150, height=1800,
    width=1800)
p <- DoHeatmap(male_np_pep1_cltmr1_sub, features = top40$gene, slot = 'scale.data')
print(p)
dev.off()

male_markers <- male_markers %>% group_by(cluster)
write.csv(male_markers, './output/de_files/neuron_male_np_pep1_cltmr1_all_receptors.csv')

            