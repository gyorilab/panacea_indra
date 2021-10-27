calculate_pct <- function(){
  
  for(n in all_cell_types){
    enriched_genes[, n] <- NA
  }
  
  for(n in 1:length(all_cell_types)){
    all_clusters[[all_cell_types[n]]] = subset(mDRG, cell_type == 
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
  
  write.csv(enriched_genes, './output/gene_pct_drg_clusters.csv', row.names = F)
  
  
  
  
}