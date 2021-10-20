library(dplyr)

### READ INPUT FILES

setwd('./kinase_tf_output/')

# Read TF estimation inputs
dorothea_cfa_tf <- read.csv('./DOROTHEA_CFA_TF.csv')
indra_sif_all_cfa <- read.csv('./INDRA_SIF_TF_ALL_CFA_TF.csv')
indra_sif_db_cfa<- read.csv('./INDRA_SIF_TF_DB_ONLY_CFA_TF.csv')


merged_dorothea_indra_all_df <- merge(dorothea_cfa_tf, 
                                      indra_sif_all_cfa, by='ID')

merged_dorothea_indra_db <- merge(dorothea_cfa_tf,
                                  indra_sif_db_cfa, by='ID')

#merged_indra_all_db <- merge(indra_sif_all_cfa,
#                             indra_sif_db_cfa, by ='ID')

#merged_dorothea_indra <- merge(merged_indra_all_db, dorothea_cfa_tf,
#                               by='ID')

cor(merged_dorothea_indra_all_df$NES_abs.x,
    merged_dorothea_indra_all_df$NES_abs.y)

cor(merged_dorothea_indra_db$NES_abs.x,
    merged_dorothea_indra_db$NES_abs.y)

make_cor_plot <- function(x,y,xlab,ylab,outfile){
  png(outfile, res = 300, height = 1200, width = 1200)
    # Creating the plot
    plot(x, y, pch = 19, col = "lightblue", xlab = xlab,
         ylab = ylab, main = paste("Correlation:", round(cor(x, y), 2)))
    
    # Regression line
    abline(lm(y ~ x), col = "red", lwd = 3)
    
  dev.off()

}

make_cor_plot(merged_dorothea_indra_db$NES_abs.x, merged_dorothea_indra_db$NES_abs.y,
              'dorothea', 'indra_db_only','../kinase_tf_output/corr_dorothea_indra_db_TF.png')

make_cor_plot(merged_dorothea_indra_all_df$NES_abs.x, merged_dorothea_indra_all_df$NES_abs.y,
              'dorothea', 'indra','../kinase_tf_output/corr_dorothea_indra_TF.png')

