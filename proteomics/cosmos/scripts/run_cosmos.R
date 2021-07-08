## instal COSMOS
#library(devtools)
#install_github("samuelbunga/cosmosR")
library(cosmosR)
library(stringr)

CARNIVAL_options <- cosmosR::default_CARNIVAL_options()
CARNIVAL_options$solverPath <- "/dart/Panacea/cIBM/cplex/bin/x86-64_linux/cplex"
CARNIVAL_options$solver <- "cplex" #or cbc
#CARNIVAL_options$solver <- "lpSolve" #or cbc
CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2

# load PKN
meta_pkn <- cosmosR::meta_network
tf_dorothea <- cosmosR::load_tf_regulon_dorothea()

cosmos_dir <- '/dart/Panacea/cosmos/'

# Load phospho data
phos_data <- read.csv(paste0(cosmos_dir, '/Data/panacea_phospho_human_mapped.csv'))
signaling_data <- phos_CFA$fold.change
names(signaling_data) <- toupper(phos_CFA$X)


cfa_kinase <- read.csv(paste0(cosmos_dir, 'dorothea_output/CFA_kin_activity.csv'))
cfa_kinase$ID <- convert_genesymbols_to_entrezid(cfa_kinase$ID)
cfa_kinase$ID <- paste0('X', cfa_kinase$ID)
cfa_kinase_full <- cfa_kinase$NES 
names(cfa_kinase_full) <- cfa_kinase$ID
cfa_kinase_full <- cfa_kinase_full[names(cfa_kinase_full) %in% meta_pkn$source |
                                     names(cfa_kinase_full) %in% meta_pkn$target]

cfa_TF <- read.csv(paste0(cosmos_dir, 
                          'dorothea_output/CFA_tf_activity.csv'))
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
CFA <- read.csv('./Data/Proteomics/CFA_panacea_volcano_data_16plex_data.csv')
transcript_data <- CFA$fold.change
names(transcript_data) <- CFA$X

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_pkn,
                                                      signaling_data = cfa_TF_full,
                                                      metabolic_data = cfa_kinase_full,
                                                      diff_expression_data = toy_RNA,
                                                      maximum_network_depth = 15,
                                                      remove_unexpressed_nodes = TRUE,
                                                      CARNIVAL_options = CARNIVAL_options
)

CARNIVAL_options$timelimit <- 14400
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = CARNIVAL_options)

test_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metabolite_to_pubchem,
                                     measured_nodes = unique(c(names(toy_metabolic_input),
                                                               names(toy_signaling_input))),
                                     omnipath_ptm = omnipath_ptm)

CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2


test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_pkn,
                                                       signaling_data = cfa_TF_full,
                                                       metabolic_data = cfa_kinase_full,
                                                       diff_expression_data = toy_RNA,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = FALSE,
                                                       CARNIVAL_options = CARNIVAL_options
                                                       
)


CARNIVAL_options$timelimit <- 28800

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = CARNIVAL_options)


test_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metabolite_to_pubchem,
                                      measured_nodes = unique(c(names(toy_metabolic_input),
                                                                names(toy_signaling_input))),
                                      omnipath_ptm = omnipath_ptm)


full_sif <- as.data.frame(rbind(test_result_for[[1]], test_result_back[[1]]))
full_attributes <- as.data.frame(rbind(test_result_for[[2]], test_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)


network_plot <- display_node_neighboorhood(central_node = 'BCAT1', 
                                           sif = full_sif, 
                                           att = full_attributes, 
                                           n = 5)

#print(network_plot)