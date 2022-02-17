symbol2ensembl <-function(cts){
  ens <- rownames(cts)
  ens <- as.character(ens)
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('ENSEMBL'), keytype = 'SYMBOL')
  #symbols
  
  rownames(cts) <-symbols
  keep <- !is.na(rownames(cts))
  cts <- cts[keep,]
  return(cts)
}

get_subsets <- function(obj){
  
  obj <- subset(obj,
                        subset = (
                          cellID == "cLTMR1" |  cellID == "p_cLTMR2" | 
                            cellID == "PEP1" |  cellID == "PEP2" | 
                            cellID == "NP" | cellID == "SST" |
                            cellID == "NF1" | cellID == "NF2" | 
                            cellID =="NF3"))
  
  
  active.ident_neuron <- data.frame('cells' = names(obj@active.ident),
                                    'gender' = unname(obj@active.ident))
  neuron_cell_types <- data.frame('cells' = obj$cells,
                                  'cell_types' = obj$cellID)
  
  active.ident_neuron <- merge(active.ident_neuron, neuron_cell_types, by='cells')
  rownames(active.ident_neuron) <- active.ident_neuron$cells
  active.ident_neuron$gender <- NULL
  renamed_idents <- factor(active.ident_neuron$cell_types)
  names(renamed_idents) <- active.ident_neuron$cells
  obj@active.ident <- renamed_idents
  return(obj)
}




