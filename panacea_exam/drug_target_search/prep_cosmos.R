setwd('~/gitHub/panacea_indra/panacea_exam/drug_target_search/')
library(dorothea)
library(cosmosR)
library(dplyr)
library(viper)

source('./viper_functions.R')

CARNIVAL_options <- cosmosR::default_CARNIVAL_options()
CARNIVAL_options$solverPath <- "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/cplex"
CARNIVAL_options$solver <- "cplex" #or cbc
#CARNIVAL_options$solver <- "lpSolve" #or cbc
CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2

## Estimate TF
meta_pkn <- cosmosR::meta_network
# read dorothea
dorothea_df<- as.data.frame(dorothea_hs[dorothea_hs$confidence 
                                        %in% c("A","B","C"), c(3,1,4)])
dorothea_viper <- df_to_viper_regulon(distinct(dorothea_df))

#import a mapping table downloaded from uniprot
RNAseq_entrez_to_symbol <- as.data.frame(read_delim("./input/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"
names(RNAseq_entrez_to_symbol)[3] <- "HUMAN_SYMBOL"
RNAseq_entrez_to_symbol$HUMAN_SYMBOL <- stringr::str_replace(RNAseq_entrez_to_symbol$HUMAN_SYMBOL,
                                                             '_HUMAN', '')
# Read tpm
tpm_counts <- read.csv('./input/tpm_counts.csv')
rownames(tpm_counts) <- tpm_counts$X

# Keep only TPM counts columns
tpm_counts <- tpm_counts[,na.omit(stringr::str_extract(colnames(tpm_counts), '.*TPM'))]

# Extract kw_taxol and taxol(6hr)
kwx_tx <- tpm_counts[,c(5,6,7,8,9,10,13)]

# take average of kw+taxol(6hr)
kwx_tx$kw_taxol_6hr_avg <- (kwx_tx$S19.Kw.Taxol.6hr_TPM +
                              kwx_tx$S21.Kw.Taxol.6hr_TPM)/2

kwx_tx$kwx_6hr_vs_taxol_6hr_fc <- (kwx_tx$kw_taxol_6hr_avg+1)/(kwx_tx$S6.Kw.6hr_TPM+1)
                              
kwx_tx_6hr <- data.frame(row.names = rownames(kwx_tx),
                         'kwx_6hr_vs_taxol_6hr_FC' = kwx_tx$kwx_6hr_vs_taxol_6hr_fc)

kwx_tx_6hr$HUMAN_SYMBOL <- rownames(kwx_tx_6hr)

RNA_differential_analysis <- merge(kwx_tx_6hr, 
                                   RNAseq_entrez_to_symbol[,c(1,3)], 
                                   by='HUMAN_SYMBOL')
eset <- RNA_differential_analysis$kwx_6hr_vs_taxol_6hr_FC
names(eset) <- RNA_differential_analysis$HUMAN_SYMBOL

#run viper
tf_activity <- run_viper(eset, dorothea_viper, 'kwx_tx_6hr')

tf_activity$ID <- convert_genesymbols_to_entrezid(tf_activity$ID)
tf_activity$ID <- paste0('X', tf_activity$ID)
tf <- tf_activity$NES 
names(tf) <- tf_activity$ID
tf <- tf[names(tf) %in% meta_pkn$source | 
           names(tf) %in% meta_pkn$target]

RNA_differential_analysis$HUMAN_SYMBOL <- convert_genesymbols_to_entrezid(RNA_differential_analysis$HUMAN_SYMBOL)
RNA_differential_analysis$HUMAN_SYMBOL <- paste0('X', RNA_differential_analysis$HUMAN_SYMBOL)
rna <- RNA_differential_analysis$kwx_6hr_vs_taxol_6hr_FC
names(rna) <- RNA_differential_analysis$HUMAN_SYMBOL



kw_hits <- read.csv('./input/210809-kinomescan-100nM-10uM-panacea-prelim-hits-KW.csv')
kw_inhibit <- read.csv('./output/kw_hits.txt', header = F)
kw_hits <- kw_hits[,c(4,5)]
kw_hits <- kw_hits[kw_hits$Entrez.Gene.Symbol %in% kw_inhibit$V1, ]

# Check this step again
kw_hits <- kw_hits[!duplicated(kw_hits$Entrez.Gene.Symbol), ]
kw_hits$Entrez.Gene.Symbol <- convert_genesymbols_to_entrezid(kw_hits$Entrez.Gene.Symbol)
kw_hits <- kw_hits[!is.na(kw_hits$Entrez.Gene.Symbol), ]
kw_hits$Entrez.Gene.Symbol <- paste0('X', kw_hits$Entrez.Gene.Symbol)
kw_hits <- kw_hits[kw_hits$Entrez.Gene.Symbol %in% meta_pkn$source |
                   kw_hits$Entrez.Gene.Symbol %in% meta_pkn$target,]

kw <- kw_hits$Percent.Control
names(kw) <- kw_hits$Entrez.Gene.Symbol

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_pkn,
                                                      signaling_data = tf,
                                                      metabolic_data = kw,
                                                      diff_expression_data = rna,
                                                      maximum_network_depth = 10,
                                                      remove_unexpressed_nodes = TRUE,
                                                      CARNIVAL_options = CARNIVAL_options
)
