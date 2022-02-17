library(pheatmap)
library(RColorBrewer)
setwd('~/gitHub/panacea_indra/cellphone_db/output/')
meta_file <- '../output/meta_data/neuro_Healthy_meta_table.tsv'
pvalues_file <- './cellphonedb_output/neuro_healthy/pvalues.txt'
meta_sep='\t'
pvalues_sep = '\t'
pvalue = 0.05
col1 = "dodgerblue4"
col2 = 'peachpuff'
col3 = 'deeppink4'

heatmaps_plot = function(meta_file, pvalues_file, count_filename, log_filename, count_network_filename, interaction_count_filename, count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network
  cell_types <- c("Keratinocytes", "T cells", "Mast cells", "LCs", "DC2",
                  "DC1", "Dermal Macs", "Macs4", "Macs3", "Macs2", "Macs1",
                  "RM", "Neutrophils", "NP", "PEP1", "PEP2", "cLTMR1",
                  "p_cLTMR2", "SST", "NF1", "NF2", "NF3")
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)

  breaksList <- seq(0, 250, by = 25)
  colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList))
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]


  split_sep = '\\|'
  join_sep = '|'

  pairs1_all = unique(meta[,2])

  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
        pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))

  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]

    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))

    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }

  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)

  #######   count interactions

  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]

    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]

    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))

  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    count_matrix <- as.data.frame(count_matrix)
    count_matrix <- as.matrix(count_matrix[cell_types, cell_types])
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)

    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )

    pheatmap(count_matrix, gaps_row = 13, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = F, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename,
             breaks = breaksList)

    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, gaps_row = 13, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = F, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}