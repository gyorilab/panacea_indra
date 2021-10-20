library(dplyr)
library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)


wd = '~/gitHub/panacea_indra/panacea_exam/transcriptomics/TPM_calculator/'
source(paste0(wd, '../functions.R'))
files <- na.omit(stringr::str_extract(list.files(paste0(wd)), 
                                      ".*_genes.out"))

gene_table <- read.table(paste0(wd, '../inputs/gene_table.tsv'), sep='\t',
                         header = T)
gene_table$Gene.stable.ID <- NULL
gene_table$Gene.stable.ID.version <- NULL
gene_table <- gene_table[!duplicated(gene_table$Gene.name), ]

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
mtx_symbol <- as.data.frame(mtx_symbol)
mtx_symbol$Gene.name <- rownames(mtx_symbol)
mtx_symbol <- merge(mtx_symbol, gene_table, by='Gene.name', all.x=T)
mtx_symbol <- as.data.frame(mtx_symbol)
rownames(mtx_symbol) <- mtx_symbol$Gene.name
mtx_symbol$Gene.name <- NULL
mtx_symbol <- mtx_symbol[, c(ncol(mtx_symbol), 1:ncol(mtx_symbol)-1)]

write.csv(mtx_symbol, paste0(wd, '../outputs/tpm_counts.csv'), 
          col.names = T, row.names = T)

# Save baseline expression
write.csv(mtx_symbol[na.omit(stringr::str_extract(colnames(mtx_symbol), 
                                                  ".*no-drug.*"))],
          paste0(wd, '../outputs/baseline_expression.csv'),
          col.names = T, row.names = T)
# Get the average of the samples
# (S1, S3), (S16, S18), (S19, S21), (S22, S24), 

# No-drug
S1_S3 <- mtx_symbol[na.omit(stringr::str_extract(colnames(merged_df), 
                                                "S\\d+-no-drug.*"))]
S1_S3$Avg_no_drug <- (S1_S3$`S1-no-drug_TPM`+S1_S3$`S3-no-drug_TPM`)/2
S1_S3$Gene_desc <- mtx_symbol$Gene.description
  
# Kw+Taxol - 1Hr
S16_S18 <- mtx_symbol[na.omit(stringr::str_extract(colnames(merged_df), 
                                                  "S\\d+-Kw\\+Taxol-1hr.*"))]
S16_S18$Avg_Kw_Taxol_1hr <- (S16_S18$`S16-Kw+Taxol-1hr_TPM`+
                               S16_S18$`S18-Kw+Taxol-1hr_TPM`)/2
S16_S18$Gene_desc <- mtx_symbol$Gene.description

# Kw+Taxol - 6 Hr
S19_S21 <- mtx_symbol[na.omit(stringr::str_extract(colnames(merged_df), 
                                                  "S\\d+-Kw\\+Taxol-6hr.*"))]
S19_S21$Avg_Kw_Taxol_6hr <- (S19_S21$`S19-Kw+Taxol-6hr_TPM`+
                               S19_S21$`S21-Kw+Taxol-6hr_TPM`)/2
S19_S21$Gene_desc <- mtx_symbol$Gene.description

# Kw+Taxol - 24 Hr
S22_S24 <- mtx_symbol[na.omit(stringr::str_extract(colnames(merged_df), 
                                                  "S\\d+-Kw\\+Taxol-24hr.*"))]
S22_S24$Avg_Kw_Taxol_24hr <- (S22_S24$`S22-Kw+Taxol-24hr_TPM`+ 
                                S22_S24$`S24-Kw+Taxol-24hr_TPM`)/2

S22_S24$Gene_desc <- mtx_symbol$Gene.description

## Make Comparisons
dir.create(paste0(wd, '../outputs/'), showWarnings = F)

# Control (samples 1&3 -avg TPM) vs KW 1 h
Kw_1hr <- as.data.frame(mtx_symbol[,c("S4-Kw-1hr_Read", "S4-Kw-1hr_TPM")])
S1_S3_vs_Kw_1hr <- cbind(S1_S3, Kw_1hr)
S1_S3_vs_Kw_1hr$foldchange <- (S1_S3$Avg_no_drug+1)/(Kw_1hr$`S4-Kw-1hr_TPM`+1)
write.csv(S1_S3_vs_Kw_1hr, paste0(wd, '../outputs/', 'S1_S3_vs_Kw_1hr.csv'))

# Control (samples 1&3 -avg TPM) vs KW 6 h
Kw_6hr <- as.data.frame(mtx_symbol[,c("S6-Kw-6hr_Read", "S6-Kw-6hr_TPM" )])
S1_S3_vs_Kw_6hr <- cbind(S1_S3, Kw_6hr)
S1_S3_vs_Kw_6hr$foldchange <- (S1_S3$Avg_no_drug+1)/(Kw_6hr$`S6-Kw-6hr_TPM`+1)
write.csv(S1_S3_vs_Kw_6hr, paste0(wd, '../outputs/', 'S1_S3_vs_Kw_6hr.csv'))

# Control (samples 1&3 -avg TPM) vs KW 24 h
Kw_24hr <- as.data.frame(mtx_symbol[,c("S8-Kw-24hr_Read", "S8-Kw-24hr_TPM" )])
S1_S3_vs_Kw_24hr <- cbind(S1_S3, Kw_24hr)
S1_S3_vs_Kw_24hr$foldchange <- (S1_S3$Avg_no_drug+1)/(Kw_24hr$`S8-Kw-24hr_TPM`+1)
write.csv(S1_S3_vs_Kw_24hr, paste0(wd, '../outputs/', 'S1_S3_vs_Kw_24hr.csv'))

# Control (samples 1&3 -avg TPM) vs Taxol 1 h
Taxol_1hr <- as.data.frame(mtx_symbol[,c("S10-Taxol-1hr_Read", "S10-Taxol-1hr_TPM" )])
S1_S3_vs_Taxol_1hr <- cbind(S1_S3, Taxol_1hr)
S1_S3_vs_Taxol_1hr$foldchange <- (S1_S3$Avg_no_drug+1)/(Taxol_1hr$`S10-Taxol-1hr_TPM`+1)
write.csv(S1_S3_vs_Taxol_1hr, paste0(wd, '../outputs/', 'S1_S3_vs_Taxol_1hr.csv'))

# Control (samples 1&3 -avg TPM) vs Taxol 6 h
Taxol_6hr <- as.data.frame(mtx_symbol[,c("S12-Taxol-6hr_Read", "S12-Taxol-6hr_TPM")])
S1_S3_vs_Taxol_6hr <- cbind(S1_S3, Taxol_6hr)
S1_S3_vs_Taxol_6hr$foldchange <- (S1_S3$Avg_no_drug+1)/(Taxol_6hr$`S12-Taxol-6hr_TPM`+1)
write.csv(S1_S3_vs_Taxol_6hr, paste0(wd, '../outputs/', 'S1_S3_vs_Taxol_6hr.csv'))

# Control (samples 1&3 -avg TPM) vs Taxol 24 h
Taxol_24hr <- as.data.frame(mtx_symbol[,c("S14-Taxol-24hr_Read", "S14-Taxol-24hr_TPM")])
S1_S3_vs_Taxol_24hr <- cbind(S1_S3, Taxol_24hr)
S1_S3_vs_Taxol_24hr$foldchange <- (S1_S3$Avg_no_drug+1)/(Taxol_24hr$`S14-Taxol-24hr_TPM`+1)
write.csv(S1_S3_vs_Taxol_24hr, paste0(wd, '../outputs/', 'S1_S3_vs_Taxol_24hr.csv'))

# Control (samples 1&3 -avg TPM) vs KW+Taxol 1 h (avgTPM)
S1_S3_vs_Kw_Taxol_1hr <- cbind(S1_S3, S16_S18)
S1_S3_vs_Kw_Taxol_1hr$foldchange <- (S1_S3$Avg_no_drug+1)/(S16_S18$Avg_Kw_Taxol_1hr+1)
write.csv(S1_S3_vs_Kw_Taxol_1hr, paste0(wd, '../outputs/', 'S1_S3_vs_Kw_Taxol_1hr.csv'))


# Control (samples 1&3 -avg TPM) vs KW+Taxol 6 h (avgTPM)
S1_S3_vs_Kw_Taxol_6hr <- cbind(S1_S3, S19_S21)
S1_S3_vs_Kw_Taxol_6hr$foldchange <- (S1_S3$Avg_no_drug+1)/(S19_S21$Avg_Kw_Taxol_6hr+1)
write.csv(S1_S3_vs_Kw_Taxol_6hr, paste0(wd, '../outputs/', 'S1_S3_vs_Kw_Taxol_6hr.csv'))


# Control (samples 1&3 -avg TPM) vs KW+Taxol 24 h (avgTPM)
S1_S3_vs_Kw_Taxol_24hr <- cbind(S1_S3, S22_S24)
S1_S3_vs_Kw_Taxol_24hr$foldchange <- (S1_S3$Avg_no_drug+1)/(S22_S24$Avg_Kw_Taxol_24hr+1)
write.csv(S1_S3_vs_Kw_Taxol_24hr, paste0(wd, '../outputs/', 'S1_S3_vs_Kw_Taxol_24hr.csv'))


# Taxol 1 h  vs KW+Taxol 1 h (avgTPM)
Taxol_1hr_vs_Kw_Taxol_1hr <- cbind(Taxol_1hr, S16_S18)
Taxol_1hr_vs_Kw_Taxol_1hr$foldchange <- (Taxol_1hr$`S10-Taxol-1hr_TPM`+1)/(S16_S18$Avg_Kw_Taxol_1hr)
write.csv(Taxol_1hr_vs_Kw_Taxol_1hr, paste0(wd, '../outputs/', 'Taxol_1hr_vs_Kw_Taxol_1hr.csv'))


# Taxol 1 h  vs KW+Taxol 6 h (avgTPM)
Taxol_1hr_vs_Kw_Taxol_6hr <- cbind(Taxol_1hr, S19_S21)
Taxol_1hr_vs_Kw_Taxol_6hr$foldchange <- (Taxol_1hr$`S10-Taxol-1hr_TPM`)/(S19_S21$Avg_Kw_Taxol_6hr)
write.csv(Taxol_1hr_vs_Kw_Taxol_6hr, paste0(wd, '../outputs/', 'Taxol_1hr_vs_Kw_Taxol_6hr.csv'))


# Taxol 1 h  vs KW+Taxol 24 h (avgTPM)
Taxol_1hr_vs_Kw_Taxol_24hr <- cbind(Taxol_1hr, S22_S24)
Taxol_1hr_vs_Kw_Taxol_24hr$foldchange <- (Taxol_1hr$`S10-Taxol-1hr_TPM`)/(S22_S24$Avg_Kw_Taxol_24hr)
write.csv(Taxol_1hr_vs_Kw_Taxol_24hr, paste0(wd, '../outputs/', 'Taxol_1hr_vs_Kw_Taxol_24hr.csv'))

