library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

source(paste0(wd, 'functions.R'))
wd = '~/gitHub/panacea_indra/panacea_exam/transcriptomics/'
files <- na.omit(stringr::str_extract(list.files(paste0(wd, 'TPM_calculator')), 
                                      ".*_genes.out"))

for (f in files) {
  if(files[1] == f){
    df <- read.csv(paste0(wd, f), sep='\t')
    df <- df[c('Gene_Id', 'TPM')]
    sample_name <- stringr::str_match(f, "(trimmed_)(\\w+.\\-[a-zA-Z\\-0-9\\+]+)")[3]
    colnames(df)[2] <- sample_name
    print(sample_name)
    merged_df <- df
  }else{
    df <- read.csv(paste0(wd, f), sep='\t')
    df <- df[c('Gene_Id', 'TPM')]
    sample_name <- stringr::str_match(f, "(trimmed_)(\\w+.\\-[a-zA-Z\\-0-9\\+]+)")[3]
    colnames(df)[2] <- sample_name
    print(sample_name)
    merged_df <- merge(merged_df, df, by='Gene_Id')
  }
  
}

rownames(merged_df) <- merged_df$Gene_Id
merged_df$Gene_Id <- NULL
mtx <- as.matrix(merged_df)
mtx_symbol <- ensembl2symbol(mtx)

