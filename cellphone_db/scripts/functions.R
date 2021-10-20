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
