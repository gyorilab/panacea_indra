
get_enriched_genes <- function(){
  wd <- ('/Users/sbunga/gitHub/panacea_indra/pain_model/')
  setwd(wd)
  
  human_mouse <- read.csv('./data/human_mouse_orth.txt', sep='\t')
  colnames(human_mouse) <- c('HUMAN_SYMBOL', 'MOUSE_SYMBOL')
  
  df <- readxl::read_xlsx(paste0(wd, 'data/Primary_mouse/other/mouse_10-fold_high_genes_expression_matrix.xlsx'))
  enriched_genes <- df[,c(1,2,3,4)]
  colnames(enriched_genes) <- c('HUMAN_SYMBOL', 'Description', 'Subcellular_localization',
                                'Functional_type_of_protein')
  enriched_genes  <- merge(enriched_genes, human_mouse, by.y ='HUMAN_SYMBOL')
  enriched_genes <- enriched_genes[!enriched_genes[ , 5] == '', ]
  enriched_genes <- enriched_genes[!duplicated(enriched_genes$MOUSE_SYMBOL), ]
  enriched_genes$HUMAN_SYMBOL <- NULL
  enriched_genes <- enriched_genes[, c(4, 1, 2, 3)]
  return(enriched_genes)
}


calculate_pct <- function(obj, enriched_genes){
  
  # Subsetting the obj object and extracting the cells
  # from each cluster
  obj@meta.data$cell_type[obj@meta.data$cell_type == ''] = 'Unknown'
  all_cell_types <- names(table(obj@meta.data$cell_type))
  all_clusters = list()
  
  for(n in all_cell_types){
    enriched_genes[, n] <- NA
  }
  
  for(n in 1:length(all_cell_types)){
    all_clusters[[all_cell_types[n]]] = subset(obj, cell_type == 
                                                 all_cell_types[n])
  }
  
  for(g in 1:nrow(enriched_genes)){
    for(cluster in all_cell_types){
      if(enriched_genes$MOUSE_SYMBOL[g] %in% rownames(all_clusters[[cluster]])){
        exp <- FetchData(all_clusters[[cluster]], enriched_genes$MOUSE_SYMBOL[g])
        # Calculating the percentage
        pct <- as.matrix(colMeans(exp  > 0))*100
        enriched_genes[g, cluster] <- round(pct, 2)
      }else{
        enriched_genes[g, cluster] <- NA
      }
    }
  }
  return(enriched_genes)
}
