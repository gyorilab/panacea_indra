
## instal COSMOS
library(visNetwork)
library(devtools)
#install_github("samuelbunga/cosmosR")
library(cosmosR)
library(stringr)




CARNIVAL_options <- cosmosR::default_CARNIVAL_options()
CARNIVAL_options$solverPath <- "/Applications/CPLEX_Studio201/cplex/bin/x86-64_osx/cplex"
CARNIVAL_options$solver <- "cplex" #or cbc
#CARNIVAL_options$solver <- "lpSolve" #or cbc
CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2


# load PKN
meta_pkn <- cosmosR::meta_network
tf_dorothea <- cosmosR::load_tf_regulon_dorothea()
indra_sif <- read.csv('./output/indra_sif.csv')
indra_sif <- dplyr::select("agA_ns", "agA_id", "agA_name",
                           "agB_ns", "agB_id", "agB_name")
cosmos_dir <- '~/gitHub/panacea_indra/proteomics/cosmos/kinase_tf_output/'

cfa_kinase <- read.csv(paste0(cosmos_dir, 'kinase_tf_output/OP_CFA_KINASE_ACTIVITY.csv'))
cfa_kinase$ID <- convert_genesymbols_to_entrezid(cfa_kinase$ID)
cfa_kinase$ID <- paste0('X', cfa_kinase$ID)
cfa_kinase_full <- cfa_kinase$NES 
names(cfa_kinase_full) <- cfa_kinase$ID
cfa_kinase_full <- cfa_kinase_full[names(cfa_kinase_full) %in% meta_pkn$source |
                                     names(cfa_kinase_full) %in% meta_pkn$target]

cfa_TF <- read.csv(paste0(cosmos_dir, 
                          'kinase_tf_output/DOROTHEA_CFA_TF.csv'))
cfa_TF$ID <- convert_genesymbols_to_entrezid(cfa_TF$ID)
cfa_TF$ID <- paste0('X', cfa_TF$ID)
cfa_TF_full <- cfa_TF$NES 
names(cfa_TF_full) <- cfa_TF$ID
cfa_TF_full <- cfa_TF_full[names(cfa_TF_full) %in% meta_pkn$source |
                             names(cfa_TF_full) %in% meta_pkn$target]

#entrez_ids <- str_extract(names(signaling_data), '[A-Za-z0-9]+')
#gs <- unlist(fastGene::convert_symbols(entrez_ids))
#cosmosR::convert_genesymbols_to_entrezid(gs)
#head(names(signaling_data))
#cosmosR::convert_genesymbols_to_entrezid(names(signaling_data))

# Load transcript data
CFA <- read.csv(paste0(cosmos_dir, 'Data/CFA_panacea_volcano_data_16plex_data.csv'))
CFA_unique <- CFA[!duplicated(CFA$X),]
HGNC_symbols <- fastGene::convert_symbols(CFA_unique$X)
HGNC_symbols[sapply(HGNC_symbols, is.null)] <- NA
HGNC_symbols <- unlist(HGNC_symbols, use.names = F)
CFA_unique$HGNC <- HGNC_symbols
CFA_unique  <- CFA_unique[-which(CFA_unique$HGNC == 'None'), ]
CFA_unique  <- CFA_unique[!is.na(CFA_unique$HGNC), ]
CFA_unique$ID <- convert_genesymbols_to_entrezid(CFA_unique$HGNC)
CFA_unique <- CFA_unique[!duplicated(CFA_unique$ID), ]
transcript_data <- CFA_unique$fold.change
names(transcript_data) <- paste0('X', CFA_unique$ID)
saveRDS(transcript_data, '../transcript_data.RDS')
save.image('../pre-process.RData')





test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_pkn,
                                                      signaling_data = cfa_TF_full,
                                                      metabolic_data = cfa_kinase_full,
                                                      diff_expression_data = transcript_data,
                                                      maximum_network_depth = 15,
                                                      remove_unexpressed_nodes = TRUE,
                                                      CARNIVAL_options = CARNIVAL_options
)



CARNIVAL_options$timelimit <- 14400
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2





test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for, CARNIVAL_options=CARNIVAL_options)
saveRDS(test_result_for, '../output/test_resuly_for.Rds')




test_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metabolite_to_pubchem,
                                     measured_nodes = unique(c(names(kw),
                                                               names(tf))),
                                     omnipath_ptm = omnipath_ptm)


## Tutorial section: metabolism to signaling 


CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2



test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_pkn,
                                                       signaling_data = tf,
                                                       metabolic_data = kw,
                                                       diff_expression_data = rna,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = FALSE,
                                                       CARNIVAL_options = CARNIVAL_options
                                                       
)



CARNIVAL_options$timelimit <- 28800
t <- test_back
t$meta_network <- meta_pkn
test_result_back <- run_COSMOS_metabolism_to_signaling(data = t,
                                                       CARNIVAL_options = CARNIVAL_options)





test_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metabolite_to_pubchem,
                                      measured_nodes = unique(c(names(toy_metabolic_input),
                                                                names(toy_signaling_input))),
                                      omnipath_ptm = omnipath_ptm)



## Tutorial section: Merge forward and backward networks and visualise network


full_sif <- as.data.frame(rbind(test_result_for[[1]], test_result_back[[1]]))
full_attributes <- as.data.frame(rbind(test_result_for[[2]], test_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)




network_plot <- display_node_neighboorhood(central_node = 'BCAT1', 
                                           sif = full_sif, 
                                           att = full_attributes, 
                                           n = 5)

print(network_plot)
sessionInfo()

full_sif <- readRDS('./output/sep-4/full_sif.Rds')
full_attributes <- readRDS('./output/sep-4/full_attributes.Rds')
write.csv(full_sif, './output/sep-4/full_sif.csv', col.names = T, row.names = F)
write.csv(full_attributes, './output/sep-4/full_attributes.csv', col.names = T, row.names = F)

network_plot <- display_node_neighboorhood(central_node = 'JAK3', 
                                           sif = full_sif, 
                                           att = full_attributes, 
                                           n = 5)
#print(network_plot)
visSave(network_plot, 
        file = "./output/sep-4/network.html", 
        background = "white")


