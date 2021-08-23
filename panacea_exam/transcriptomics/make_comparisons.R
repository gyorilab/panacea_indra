library(dplyr)
library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

source(paste0(wd, 'functions.R'))
wd = '~/gitHub/panacea_indra/panacea_exam/transcriptomics/TPM_calculator/'
files <- na.omit(stringr::str_extract(list.files(paste0(wd)), 
                                      ".*_genes.out"))

for (f in files) {
  if(files[1] == f){
    df <- read.csv(paste0(wd, f), sep='\t')
    df <- df[c('Gene_Id', 'Reads','TPM')]
    sample_name <- stringr::str_match(f, "(trimmed_)(\\w+.\\-[a-zA-Z\\-0-9\\+]+)")[3]
    colnames(df)[2] <- paste0(sample_name, '_Read')
    colnames(df)[3] <- paste0(sample_name,'_TPM')
    print(sample_name)
    merged_df <- df
  }else{
    df <- read.csv(paste0(wd, f), sep='\t')
    df <- df[c('Gene_Id', 'Reads', 'TPM')]
    sample_name <- stringr::str_match(f, "(trimmed_)(\\w+.\\-[a-zA-Z\\-0-9\\+]+)")[3]
    colnames(df)[2] <- paste0(sample_name, '_Read')
    colnames(df)[3] <- paste0(sample_name,'_TPM')
    print(sample_name)
    merged_df <- merge(merged_df, df, by='Gene_Id')
  }
  
}

rownames(merged_df) <- merged_df$Gene_Id
merged_df$Gene_Id <- NULL
mtx <- as.matrix(merged_df)
mtx_symbol <- ensembl2symbol(mtx)

# Get the average of the samples
# (S1, S3), (S16, S18), (S19, S21), (S22, S24), 

# No-drug
S1_S3 <- merged_df[na.omit(stringr::str_extract(colnames(merged_df), 
                                                "S\\d+-no-drug.*"))]
S1_S3$Avg_no_drug <- (S1_S3$`S1-no-drug_TPM`+S1_S3$`S3-no-drug_TPM`)/2

# Kw+Taxol - 1Hr
S16_S18 <- merged_df[na.omit(stringr::str_extract(colnames(merged_df), 
                                                  "S\\d+-Kw\\+Taxol-1hr.*"))]
S16_S18$Avg_Kw_Taxol_1hr <- (S16_S18$`S16-Kw+Taxol-1hr_TPM`+
                               S16_S18$`S18-Kw+Taxol-1hr_TPM`)/2

# Kw+Taxol - 6 Hr
S19_S21 <- merged_df[na.omit(stringr::str_extract(colnames(merged_df), 
                                                  "S\\d+-Kw\\+Taxol-6hr.*"))]
S19_S21$Avg_Kw_Taxol_6hr <- (S19_S21$`S19-Kw+Taxol-6hr_TPM`+
                               S19_S21$`S21-Kw+Taxol-6hr_TPM`)/2

# Kw+Taxol - 24 Hr
S22_S24 <- merged_df[na.omit(stringr::str_extract(colnames(merged_df), 
                                                  "S\\d+-Kw\\+Taxol-24hr.*"))]
S22_S24$Avg_Kw_Taxol_24hr <- (S22_S24$`S22-Kw+Taxol-24hr_TPM`+ 
                                S22_S24$`S24-Kw+Taxol-24hr_TPM`)/2


## Make Comparisons
# Control (samples 1&3 -avg TPM) vs KW 1 h

