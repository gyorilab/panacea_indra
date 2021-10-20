#import the data
setwd("~/gitHub/panacea_indra/proteomics/cosmos/")

phospho_differential_analysis <- 
  read.csv("./Data/panacea_phospho_human_mapped.csv")
indra_tf <- read.csv('./input/transcription_factors.csv')

# Read indra sif dump
indra_phospho_df <- read.csv('./output/indra_phospho_df.csv')
