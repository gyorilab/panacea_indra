library(readr)
library(viper)
library(dplyr)
library(dorothea)
library(OmnipathR)

setwd("~/gitHub/panacea_indra/proteomics/cosmos/")

source('./scripts/read_TF_input.R')
source('./scripts/viper_functions.R')

## FORMATTING WHOLE INDRA SIF
# Take increase/decrease amount from INDRA and run the TF activity prediction
indra_sif_filtered <- indra_sif %>% filter(stmt_type == 'DecreaseAmount' |
                                             stmt_type == 'IncreaseAmount')

indra_sif_filtered <- indra_sif_filtered %>% select('agA_name', 
                                                    'agB_name', 
                                                    'stmt_type')

indra_sif_filtered <- indra_sif_filtered[, c("agB_name", 
                                             "agA_name",
                                             "stmt_type")]

indra_sif_filtered$pairs <- paste(indra_sif_filtered$agA_name,
                                  indra_sif_filtered$agB_name,
                                  sep='_')
indra_sif_filtered <- indra_sif_filtered[!duplicated(indra_sif_filtered$pairs),]
indra_sif_filtered$sign <- ifelse(indra_sif_filtered$stmt_type == 'IncreaseAmount', 1, -1)
indra_sif_filtered <- indra_sif_filtered[, c("agB_name","agA_name","sign")]




## FORMATTING INDRA SIF FILTERED TO TF
# import indra statements filtered to human genes and TF only
indra_sif_tf_all$sign <- ifelse(indra_sif_tf_all$stmt_type == 'IncreaseAmount', 1, -1)
indra_sif_tf_all <- indra_sif_tf_all[, c("agB","agA","sign")]

# import indra TF statements filtered to database sources only
indra_sif_tf_db$sign <- ifelse(indra_sif_tf_db$stmt_type == 'IncreaseAmount', 1, -1)
indra_sif_tf_db <- indra_sif_tf_db[, c("agB","agA","sign")]


## FORMATTING MOUSE AND HUMAN ORTHOLOGUE TABLE
# Mouse and Human gene orthologue table
colnames(mart_table) <- c('MOUSE_ID', 'MOUSE_SYMBOL', 'HUMAN_ID', 'HUMAN_SYMBOL')

# Removing rows with no value
mart_table <- mart_table[Reduce(`&`, lapply(mart_table, function(x) !x=="")),]
colnames(RNA_differential_analysis)[1] <- 'MOUSE_SYMBOL'



RNA_differential_analysis <- merge(RNA_differential_analysis, 
                                   mart_table[,c(2,4)], 
                                   by='MOUSE_SYMBOL')

remove_rows <- grep('\\w+\\-\\w+', RNA_differential_analysis$HUMAN_SYMBOL)
RNA_differential_analysis <- RNA_differential_analysis[-remove_rows, ]
RNA_differential_analysis$HUMAN_SYMBOL <- paste(RNA_differential_analysis$HUMAN_SYMBOL, 
                                               'HUMAN', sep = '_') 


#import a mapping table downloaded from uniprot
names(RNAseq_entrez_to_symbol)[1] <- "ID"
names(RNAseq_entrez_to_symbol)[3] <- "HUMAN_SYMBOL"

#this part is to merge the mapping table with the differential analysis dataframe, using entrez gene id as a common key between the two
#this way, we will have our gene entrez identifiers mapped to their corresponding symbols in the differential analysis dataframe
#of course, there are many other way to achievethis goal.
RNA_differential_analysis <- merge(RNA_differential_analysis, 
                                   RNAseq_entrez_to_symbol[,c(1,3)], 
                                   by='HUMAN_SYMBOL')

#RNA_differential_analysis <- RNA_differential_analysis[,c(8,2:7)]
#names(RNA_differential_analysis)[1] <- "ID"
RNA_differential_analysis$HUMAN_SYMBOL <- gsub("_.*","",
                                               RNA_differential_analysis$HUMAN_SYMBOL)
RNA_differential_analysis <- unique(RNA_differential_analysis)

#now we just need to format the differential analysis data into a format 
#that is compatible with viper
eset <- RNA_differential_analysis$fold.change
names(eset) <- RNA_differential_analysis$HUMAN_SYMBOL

#we also need to format the SIF dataframe into a viper format
indra_sif_all_viper <- df_to_viper_regulon(distinct(indra_sif_filtered))
indra_sif_tf_all_viper <- df_to_viper_regulon(distinct(indra_sif_tf_all))
indra_sif_tf_db_viper <- df_to_viper_regulon(distinct(indra_sif_tf_db))
dorothea_viper <- df_to_viper_regulon(distinct(dorothea_df))

run_viper <- function(eset, prior, outfile){
  #Now we estimate the TF activities using viper
  TF_activities <- as.data.frame(viper(eset = eset, regulon = prior, 
                                       minsize = 10, adaptive.size = F, 
                                       eset.filter = F, pleiotropy = T))
  TF_activities$TF <- row.names(TF_activities)
  
  #that's just to make the dataframe pretty
  TF_activities <- TF_activities[,c(2,1)]
  names(TF_activities) <- c("ID","NES")
  
  TF_activities$NES_abs <- abs(TF_activities$NES)
  TF_activities$NES_abs[which(TF_activities$NES < 0)] <- TF_activities$NES_abs[which(TF_activities$NES < 0)]*(-1)
  TF_activities$NES <- NULL
  TF_activities <- TF_activities %>% arrange(desc(TF_activities$NES_abs))
  
  write.csv(TF_activities, paste0('./kinase_tf_output/', outfile, '.csv'),
            row.names = F)
}

run_viper(eset, indra_sif_all_viper, 'INDRA_SIF_ALL_CFA_TF')
run_viper(eset, indra_sif_tf_all_viper, 'INDRA_SIF_TF_ALL_CFA_TF')
run_viper(eset, indra_sif_tf_db_viper, 'INDRA_SIF_TF_DB_ONLY_CFA_TF')
run_viper(eset, dorothea_viper, 'DOROTHEA_CFA_TF')

