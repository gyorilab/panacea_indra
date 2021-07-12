library(readr)
library(viper)
library(dplyr)
library(OmnipathR)

#/!\ /!\ /!\ /!\ /!\ /!\ /!\
# Change the working directory to were you 
# cloned the repository on your machine
#/!\ /!\ /!\ /!\ /!\ /!\ /!\
setwd("~/gitHub/panacea_indra/proteomics/cosmos/")
#/!\ /!\ /!\ /!\ /!\ /!\ /!\

source("./scripts/viper_functions.R")

########## PHOSHO and KINASE part ########

#import the data
phospho_differential_analysis <- 
  read.csv("./Data/panacea_phospho_human_mapped.csv")


#format it properly
phospho_filtered <- phospho_differential_analysis %>% 
  filter(!is.na(human_pos_psp))

phospho_filtered$psite <- paste0(toupper(phospho_filtered$GeneSymbol), '_',
                                 phospho_filtered$Residue, 
                                 phospho_filtered$human_pos_psp)

phospho_filtered <- phospho_filtered[!duplicated(phospho_filtered$psite), ]
row.names(phospho_filtered) <- phospho_filtered$psite
phospho_filtered$X <- NULL
phospho_filtered_c <- data.frame(row.names = rownames(phospho_filtered),
                               'CFA_fc' = phospho_filtered[,which(colnames(phospho_filtered) == "CFA_fc")]
                               )

# import indra sif dump
#indra_sif <- read.csv('./output/indra_sif.csv')
#indra_sif_ksn <- indra_sif[indra_sif$stmt_type %in%
#                             c("Dephosphorylation", "Phosphorylation"),]
#indra_sif_ksn <- indra_sif_ksn[-which(indra_sif_ksn$residue == ''), ]

indra_phospho_df <- read.csv('./output/indra_phospho_df.csv')
indra_phospho_df <- indra_phospho_df %>% dplyr::select('agA_id', 'agA_name',
                                                       'agB_id', 'agB_name',
                                                       'residue', 'position', 'stmt_type')
indra_phospho_df$agB_name <- paste(indra_phospho_df$agB_name,
                                   indra_phospho_df$residue, sep='_')
indra_phospho_df$agB_name <- paste(indra_phospho_df$agB_name,
                                   indra_phospho_df$position, sep='')
indra_phospho_df$sign <- ifelse(indra_phospho_df$stmt_type == "Phosphorylation", 1, -1)
indra_phospho_df <- indra_phospho_df %>% select('agB_name', 'agA_name', 'sign')


#import KSN from omnipath
omnipath_ptm <- get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% 
                               c("dephosphorylation","phosphorylation"), ]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

#Format KSN
KSN_viper <- df_to_viper_regulon(KSN)

#Format indra sif
indra_phospho_ksn <- df_to_viper_regulon(indra_phospho_df)

#run viper to get the TF activities from the phosphoproteomic data
#You can also run that on wour normalised intesity matrix of phosphosites directly,
#as long as it is formatted as a dataframe of similar format as here
#User is strongly encouraged to check the viper publication (PMID: 27322546) for more info on the process
kin_activity <- as.data.frame(viper(eset = phospho_filtered_c, regulon = indra_phospho_ksn, 
                      minsize = 5, adaptive.size = F, eset.filter = F))
kin_activity$ID <- row.names(kin_activity)

kin_activity <- kin_activity[,c(2,1)]
names(kin_activity) <- c("ID","NES")
dir.create('../panacea_indra/other/cosmos/input/', showWarnings = F)
write.csv(kin_activity,'./dorothea_output/CFA_kin_activity.csv',row.names = F)



########## RNA and TF part ########

library(dorothea)

#First we import the dorothea regulons (using only confidence A, B, and C), see dorothea publication for information on confidence levels
dorothea_df<- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"),c(3,1,4)])

#import the RNAseq data. It has entrez gene identifiers, but we need it to have gene symbols to match dorothea database, so we have
#to do some id conversion as well
RNA_differential_analysis <- as.data.frame(
  read.csv("./Data/Proteomics/CFA_panacea_volcano_data_16plex_data.csv"))


# Mouse and Human gene ortholog table
mart_table <- read.table('./input/mart_export.txt',
                         sep='\t', header = T)
colnames(mart_table) <- c('HGNC_ID', 'HGNC_SYMBOL', 'MOUSE_SYMBOL', 'MOUSE_ID')

# Removing rows with no value
mart_table <- mart_table[Reduce(`&`, lapply(mart_table, function(x) !x=="")),]
colnames(RNA_differential_analysis)[1] <- 'MOUSE_SYMBOL'

RNA_differential_analysis <- merge(RNA_differential_analysis, mart_table[,c(2,3)], by='MOUSE_SYMBOL')
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
RNA_differential_analysis <- merge(RNA_differential_analysis, RNAseq_entrez_to_symbol[,c(1,3)], by='HGNC_SYMBOL')
#RNA_differential_analysis <- RNA_differential_analysis[,c(8,2:7)]
#names(RNA_differential_analysis)[1] <- "ID"
RNA_differential_analysis$HGNC_SYMBOL <- gsub("_.*","",RNA_differential_analysis$HGNC_SYMBOL)
RNA_differential_analysis <- unique(RNA_differential_analysis)

#now we just need to format the differential analysis data into a format 
#that is compatible with viper
eset <- RNA_differential_analysis$fold.change
names(eset) <- RNA_differential_analysis$HGNC_SYMBOL

#we also need to format the dorothea dataframe into a viper format
dorothea_viper <- df_to_viper_regulon(dorothea_df)

#Now we estimate the TF activities using viper
TF_activities <- as.data.frame(viper(eset = eset, regulon = dorothea_viper, minsize = 10, adaptive.size = F, eset.filter = F, pleiotropy = T))
TF_activities$TF <- row.names(TF_activities)

#that's just to make the dataframe pretty
TF_activities <- TF_activities[,c(2,1)]
names(TF_activities) <- c("ID","NES")
write.csv(TF_activities, './dorothea_output/CFA_TF_activity.csv',row.names = F)

########## CONCLUSION

#Now you have succefully estimated kinase and TF activities from phosphoproteomic and transcriptomic
View(kin_activity)
View(TF_activities)

#You can now combined them together and use them as input for COSMOS.
#You may also leave them separated and use them a separated input and measurments in cosmos, if you lack metabolomic data
#See https://github.com/saezlab/cosmosR for more info on how to use cosmos
