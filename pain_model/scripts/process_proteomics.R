library(xlsx)
library(dplyr)


wd <- ('/Users/sbunga/gitHub/panacea_indra/pain_model/')
setwd(wd)

human_mouse <- read.csv('./data/human_mouse_orth.txt', sep='\t')
colnames(human_mouse) <- c('HUMAN_SYMBOL', 'MOUSE_SYMBOL')

df <- readxl::read_xlsx(paste0(wd, 'data/Primary_mouse/other/mouse_10-fold_high_genes_expression_matrix.xlsx'))
enriched_genes <- df[,c(1,2,3,4)]
colnames(enriched_genes) <- c('HUMAN_SYMBOL', 'Description', 'Subcellular_localization',
                              'Functional_type_of_protein')
enriched_genes  <- merge(enriched_genes, human_mouse, by.y ='HUMAN_SYMBOL')
enriched_genes <- enriched_genes[!enriched_genes[ , 5] == '', ]
enriched_genes <- enriched_genes[!duplicated(enriched_genes$MOUSE_SYMBOL), ]
enriched_genes$HUMAN_SYMBOL <- NULL
enriched_genes <- enriched_genes[, c(4, 1, 2, 3)]

panacea_16p <- read.xlsx('./data/Primary_mouse/proteomics/Panacea_16plex_oct2020_quant_data_working_forSam.xlsx',
                         sheetIndex = 1)

panacea_16p <- panacea_16p %>% select(2,3,8,9,10,11,12,13,14)
panacea_16p$DRG <- rowMeans(panacea_16p[, c(6,7,8,9)])
panacea_16p <- panacea_16p[, c(1,2,3,4,5,10)]
colnames(panacea_16p)[1] <- 'MOUSE_SYMBOL'
panacea_16p <- panacea_16p[panacea_16p$MOUSE_SYMBOL %in% enriched_genes$MOUSE_SYMBOL, ]
write.csv(panacea_16p, './output/protein_exp.csv')  

