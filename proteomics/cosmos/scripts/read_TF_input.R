setwd("~/gitHub/panacea_indra/proteomics/cosmos/")
library(dorothea)
#First we import the dorothea regulons (using only confidence A, B, and C), see dorothea publication for information on confidence levels
dorothea_df<- as.data.frame(dorothea_hs[dorothea_hs$confidence 
                                        %in% c("A","B","C"), c(3,1,4)])

# import indra sif
indra_sif <- read.csv('./output/indra_sif.csv')


# import indra statements filtered to human genes and TF only
indra_sif_tf_all <- read.csv('./output/indra_all_tf.csv', row.names = 1)

# import indra TF statements filtered to database sources only
indra_sif_tf_db <- read.csv('./output/indra_db_only_tf.csv', row.names = 1)


#import the RNAseq data. It has entrez gene identifiers, but we need it to have gene symbols 
#to match dorothea database, so we have to do some id conversion as well
RNA_differential_analysis <- as.data.frame(
  read.csv("./Data/CFA_panacea_volcano_data_16plex_data.csv"))

mart_table <- read.table('./input/mouse_human_orthologue.txt',
                         sep='\t', header = T)

#import a mapping table downloaded from uniprot
RNAseq_entrez_to_symbol <- as.data.frame(read_delim("./input/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127
