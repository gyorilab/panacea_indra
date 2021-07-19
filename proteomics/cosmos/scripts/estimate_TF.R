library(dorothea)

#First we import the dorothea regulons (using only confidence A, B, and C), see dorothea publication for information on confidence levels
dorothea_df<- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"),c(3,1,4)])

# Take increase/decrease amount from INDRA and run the TF activity prediction
# Also take TF specific statements by using filter in indra/resources/TF

# import indra sif
indra_sif <- read.csv('./output/indra_sif.csv')
indra_sif_filtered <- indra_sif %>% filter(stmt_type == 'DecreaseAmount' |
                                             stmt_type == 'IncreaseAmount')

#indra_sif_tf <- indra_sif_filtered[indra_sif_filtered$agA_name %in% dorothea_df$tf, ]
indra_sif_filtered <- indra_sif_filtered %>% select('agA_name', 
                                                    'agB_name', 
                                                    'stmt_type')

indra_sif_filtered$pairs <- paste(indra_sif_filtered$agA_name, 
                                  indra_sif_filtered$agB_name, sep='_')
indra_sif_filtered <- indra_sif_filtered[, c("agB_name", "agA_name","stmt_type","pairs" )]
indra_sif_filtered <- indra_sif_filtered[!duplicated(indra_sif_fil$pairs),]
indra_sif_filtered$sign <- ifelse(indra_sif_filtered$stmt_type == 'IncreaseAmount', 1, -1)
indra_sif_filtered <- indra_sif_filtered[, c("agB_name","agA_name","sign")]


# import indra TF statements
indra_sif_tf <- read.csv('./output/indra_db_only_tf.csv', row.names = 1)
indra_sif_tf$sign <- ifelse(indra_sif_tf$stmt_type == 'IncreaseAmount', 1, -1)
indra_sif_tf <- indra_sif_tf[, c("agB","agA","sign")]

#import the RNAseq data. It has entrez gene identifiers, but we need it to have gene symbols to match dorothea database, so we have
#to do some id conversion as well
RNA_differential_analysis <- as.data.frame(
  read.csv("./Data/CFA_panacea_volcano_data_16plex_data.csv"))


# Mouse and Human gene ortholog table
mart_table <- read.table('./input/mart_export.txt',
                         sep='\t', header = T)
colnames(mart_table) <- c('HGNC_ID', 'HGNC_SYMBOL', 'MOUSE_SYMBOL', 'MOUSE_ID')

# Removing rows with no value
mart_table <- mart_table[Reduce(`&`, lapply(mart_table, function(x) !x=="")),]
colnames(RNA_differential_analysis)[1] <- 'MOUSE_SYMBOL'

RNA_differential_analysis <- merge(RNA_differential_analysis, 
                                   mart_table[,c(2,3)], 
                                   by='MOUSE_SYMBOL')

RNA_differential_analysis$HGNC_SYMBOL <- paste(RNA_differential_analysis$HGNC_SYMBOL, 
                                               'HUMAN', sep = '_') 
#import a mapping table downloaded from uniprot
RNAseq_entrez_to_symbol <- as.data.frame(read_delim("./input/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"
names(RNAseq_entrez_to_symbol)[3] <- "HGNC_SYMBOL"

#this part is to merge the mapping table with the differential analysis dataframe, using entrez gene id as a common key between the two
#this way, we will have our gene entrez identifiers mapped to their corresponding symbols in the differential analysis dataframe
#of course, there are many other way to achievethis goal.
RNA_differential_analysis <- merge(RNA_differential_analysis, 
                                   RNAseq_entrez_to_symbol[,c(1,3)], 
                                   by='HGNC_SYMBOL')

#RNA_differential_analysis <- RNA_differential_analysis[,c(8,2:7)]
#names(RNA_differential_analysis)[1] <- "ID"
RNA_differential_analysis$HGNC_SYMBOL <- gsub("_.*","",RNA_differential_analysis$HGNC_SYMBOL)
RNA_differential_analysis <- unique(RNA_differential_analysis)

#now we just need to format the differential analysis data into a format 
#that is compatible with viper
eset <- RNA_differential_analysis$fold.change
names(eset) <- RNA_differential_analysis$HGNC_SYMBOL

#we also need to format the dorothea dataframe into a viper format
dorothea_viper <- df_to_viper_regulon(distinct(indra_sif_tf))

#Now we estimate the TF activities using viper
TF_activities <- as.data.frame(viper(eset = eset, regulon = dorothea_viper, 
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

write.csv(TF_activities, './kinase_tf_output/INDRA_SIF_TF_DB_ONLY_CFA.csv',row.names = F)

########## CONCLUSION

#Now you have succefully estimated kinase and TF activities from phosphoproteomic and transcriptomic
View(kin_activity)
View(TF_activities)

#You can now combined them together and use them as input for COSMOS.
#You may also leave them separated and use them a separated input and measurments in cosmos, if you lack metabolomic data
#See https://github.com/saezlab/cosmosR for more info on how to use cosmos
