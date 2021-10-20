library(gplots)
library(pheatmap)

myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
wd = '~/gitHub/panacea_indra/panacea_exam/transcriptomics/'

# Read annotation file
sif <- read.csv(paste0(wd, '/inputs/sif.tsv'), sep='\t', row.names = 'Sample')


# Read input
df <- readxl::read_xlsx(paste0(wd, '/inputs/List_of_genes_for_heatmaps.xlsx'))

tubulin_assc <- readxl::read_xlsx(paste0(wd, '/inputs/TUBULIN_ASSOCIATED_GENES.xlsx'),
                                  col_names = F)

# Read kinase list
kinase_list <- read.csv(paste0(wd, '/inputs/kinase.csv'), row.names = 1)

# Read regenration new lit
new_regeneration_list <- read.csv(paste0(wd, '/inputs/common_upregulated_genes_DRG.csv'))

# Read counts file
tpm <- read.csv(paste0(wd, '/outputs/tpm_counts.csv'))
rownames(tpm) <- tpm$X
tpm$X <- NULL
tpm$Gene.description <- NULL


# Remove read columns
tpm[na.omit(stringr::str_extract(colnames(tpm), '.*Read'))] <- NULL
colnames(tpm) <- stringr::str_replace(colnames(tpm), '_TPM', '')
tpm <- tpm[c("S1.no.drug", "S3.no.drug", "S10.Taxol.1hr", "S12.Taxol.6hr", "S14.Taxol.24hr",
             "S16.Kw.Taxol.1hr", "S18.Kw.Taxol.1hr", "S19.Kw.Taxol.6hr", "S21.Kw.Taxol.6hr",
             "S22.Kw.Taxol.24hr", "S24.Kw.Taxol.24hr")]

# Plot
kw_missing_genes <- c('DCLK3', 'MAPK8', 'MAPK9', 'MAPK10', 'MINK1', 'MAP3K19',
                      'MAP2K5', 'STK10')
key_kw_targets <- df$`Gene set1: Key KW targets`
key_kw_targets <- key_kw_targets[!is.na(key_kw_targets)]
key_kw_targets <- c(key_kw_targets, kw_missing_genes)
key_kw_df <- tpm[rownames(tpm) %in% key_kw_targets,]
key_kw_df <- t(scale(t(key_kw_df)))

png(paste0(wd, '/outputs/heatmaps/key_kw_targets.png'),  width = 800, height = 1000, 
    res=100)

print(pheatmap(key_kw_df, fontsize_row = 9, 
         cluster_rows = F,
         cluster_cols = F, annotation_legend = F))

dev.off()


excitability_targets <- df$`Gene set2: Excitability targets`
excitability_targets <- excitability_targets[!is.na(excitability_targets)]
excitability_targets <- tpm[rownames(tpm) %in% excitability_targets,]
excitability_targets <- t(scale(t(excitability_targets)))

png(paste0(wd, '/outputs/excitability_targets.png'),  width = 1200, height = 1800, 
    res=100)

print(pheatmap(excitability_targets, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

regeneration_targets <- df$`Gene set3: de-/regeneration linked targets`
regeneration_targets <- regeneration_targets[!is.na(regeneration_targets)]
regeneration_targets <- tpm[rownames(tpm) %in% regeneration_targets,]
regeneration_targets <- t(scale(t(regeneration_targets)))

png(paste0(wd, '/outputs/de-_regeneration_linked_targets.png'),  
    width = 800, height = 1200, 
    res=100)

print(pheatmap(regeneration_targets, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()


hDRG_enriched_kinases <- df$`Gene ste5: hDRG enriched kinases`
hDRG_enriched_kinases <- hDRG_enriched_kinases[!is.na(hDRG_enriched_kinases)]
hDRG_enriched_kinases <- tpm[rownames(tpm) %in% hDRG_enriched_kinases,]
hDRG_enriched_kinases <- t(scale(t(hDRG_enriched_kinases)))

hDRG_enriched_kinases_1 <- hDRG_enriched_kinases[1:160, ]
hDRG_enriched_kinases_2 <- hDRG_enriched_kinases[161:320,]
hDRG_enriched_kinases_3 <- hDRG_enriched_kinases[321:480, ]
hDRG_enriched_kinases_4 <- hDRG_enriched_kinases[480:640, ]
hDRG_enriched_kinases_5 <- hDRG_enriched_kinases[641:nrow(hDRG_enriched_kinases),]

png(paste0(wd, '/outputs/heatmaps/hDRG_enriched_kinases.png'),  
    width = 1300, height = 6900, 
    res=400)

print(pheatmap(hDRG_enriched_kinases, fontsize_row = 2, 
               cluster_rows = T,
               cluster_cols = F, annotation_legend = F))

dev.off()


png(paste0(wd, '/outputs/hDRG_enriched_kinases_1.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(hDRG_enriched_kinases_1, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/hDRG_enriched_kinases_2.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(hDRG_enriched_kinases_2, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/hDRG_enriched_kinases_3.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(hDRG_enriched_kinases_3, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/hDRG_enriched_kinases_4.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(hDRG_enriched_kinases_4, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/hDRG_enriched_kinases_5.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(hDRG_enriched_kinases_5, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

# Plot Tubulin associated genes
tubulin <- tpm[rownames(tpm) %in% tubulin_assc$...1,]
tubulin <- t(scale(t(tubulin)))

png(paste0(wd, '/outputs/heatmaps/tublin_associated_genes.png'),  
    width = 800, height = 1000, 
    res=150)

print(pheatmap(tubulin, fontsize_row = 5, 
               cluster_rows = T,
               cluster_cols = F, annotation_legend = F))

dev.off()

# plot regeneration
human_mouse <- read.csv(paste0(wd, 'inputs/Human_Mouse.tsv'), sep = '\t')
colnames(new_regeneration_list) <- 'Mouse.gene.name'
new_regeneration_list <- merge(human_mouse, new_regeneration_list, by='Mouse.gene.name')

regeneration <- tpm[rownames(tpm) %in% new_regeneration_list$Gene.name,]
regeneration <- t(scale(t(regeneration)))

png(paste0(wd, '/outputs/heatmaps/regeneration_full.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(regeneration, fontsize_row = 2, 
               cluster_rows = T,
               cluster_cols = F, annotation_legend = F))

dev.off()

# plot kinase 
All_kinases <- kinase_list
All_kinases <- All_kinases[!is.na(All_kinases)]
All_kinases <- tpm[rownames(tpm) %in% All_kinases,]
All_kinases <- t(scale(t(All_kinases)))

All_kinases_1 <- All_kinases[1:160, ]
All_kinases_2 <- All_kinases[161:320,]
All_kinases_3 <- All_kinases[321:480, ]
All_kinases_4 <- All_kinases[481:640, ]
All_kinases_5 <- All_kinases[641:nrow(All_kinases), ]


png(paste0(wd, '/outputs/heatmaps/kinases.png'),  
    width = 1400, height = 6900, 
    res=400)

print(pheatmap(All_kinases, fontsize_row = 2, 
               cluster_rows = T,
               cluster_cols = F, annotation_legend = F))

dev.off()


png(paste0(wd, '/outputs/kinases_1.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(All_kinases_1, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/kinases_2.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(All_kinases_2, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()


png(paste0(wd, '/outputs/kinases_3.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(All_kinases_3, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/kinases_4.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(All_kinases_4, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()

png(paste0(wd, '/outputs/kinases_5.png'),  
    width = 800, height = 1700, 
    res=100)

print(pheatmap(All_kinases_5, fontsize_row = 9, 
               cluster_rows = F,
               cluster_cols = F, annotation_legend = F))

dev.off()
