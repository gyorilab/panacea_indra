ensembl2symbol <-function(cts){
  ens <- rownames(cts)
  ens %<>% as.character %>% gsub("\\.\\d+", "", .)
  
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  rownames(cts) <-symbols
  keep <- !is.na(rownames(cts))
  cts <- cts[keep,]
  return(cts)
}

