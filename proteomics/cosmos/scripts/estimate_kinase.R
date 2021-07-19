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
source('./scripts/read_kinase_input.R')

########## PHOSHO and KINASE part ########



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

#Format indra prior
indra_phospho_ksn <- df_to_viper_regulon(indra_phospho_df)

#run viper to get the TF activities from the phosphoproteomic data
#You can also run that on wour normalised intesity matrix of phosphosites directly,
#as long as it is formatted as a dataframe of similar format as here
#User is strongly encouraged to check the viper publication (PMID: 27322546) for more info on the process

run_viper <- function(phospho_data, prior, outfile){
    kin_activity <- as.data.frame(viper(eset = phospho_data, regulon = prior, 
                        minsize = 5, adaptive.size = F, eset.filter = F))
    kin_activity$ID <- row.names(kin_activity)
  
    kin_activity <- kin_activity[,c(2,1)]
    names(kin_activity) <- c("ID","NES")
    dir.create('../panacea_indra/other/cosmos/kinase_tf_output/', 
               showWarnings = F)
    
    kin_activity$NES_abs <- abs(kin_activity$NES)
    kin_activity$NES_abs[which(kin_activity$NES < 0)] <- 
      kin_activity$NES_abs[which(kin_activity$NES < 0)]*(-1)
    kin_activity$NES_abs <- NULL
    kin_activity <- kin_activity %>% arrange(desc(kin_activity$NES))
    write.csv(kin_activity, paste('./kinase_tf_output/', outfile, '.csv', sep=''),row.names = F)
}

run_viper(phospho_filtered_c, indra_phospho_ksn, 'INDRA_CFA_KINASE_ACTIVITY')
run_viper(phospho_filtered_c, KSN_viper, 'OP_CFA_KINASE_ACTIVITY')
