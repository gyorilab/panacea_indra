###########
# Code to make custom interactions file 
# for drawing heatmaps
###########
all_samples <- c('healthy', 'incision', 'zymo', 'saline')

for (s in all_samples) {
  
  meta = read.csv(file = paste0('./input/neuro_',s,'_meta_table.tsv'), comment.char = '', sep='\t')
  all_intr = read.table(paste0('./Outputs/neuro_',s,'_cellphonedb_nature/pvalues.txt'), header=T, stringsAsFactors = F, 
                        sep='\t', comment.char = '', check.names = F)
  main_cols <- colnames(all_intr)[1:11]
  
  neuron_cells <- c('cLTMR1', 'cLTMR2', 'NP', 'PEP1', 'PEP2', 'SST')
  
  
  
  intr_pairs = all_intr$interacting_pair
  cell_types = all_intr[,-c(1:11)]
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  
  
  # Creating pairs from immune to neurons
  immune_to_neurons <- c()
  neuron_to_immune <- c()
  neuron_to_neuron <- c()
  
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      if (pairs1_all[j] %in% neuron_cells & !pairs1_all[i] %in% neuron_cells){
        immune_to_neurons = c(immune_to_neurons, 
                              paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
      }
  # Creating neuron -> immune
  for (i in 1:length(neuron_cells))
    for (j in 1:length(pairs1_all))
      if(!pairs1_all[j] %in% neuron_cells){
        neuron_to_immune = c(neuron_to_immune, 
                             paste(neuron_cells[i],pairs1_all[j],sep=join_sep))
      }
  
  immune_neuron_b_df <- all_intr[, c(main_cols, immune_to_neurons)]
  immune_neuron_b_df <- subset(immune_neuron_b_df,  receptor_b=='True')
  #write.csv(immune_neuron_b_df, './Outputs/neuro_healthy_cellphonedb_nature/immune_neuron.csv')
  
  immune_neuron_a_df <- all_intr[, c(main_cols, neuron_to_immune)]
  immune_neuron_a_df <- subset(immune_neuron_a_df, receptor_a=='True')
  #write.csv(neuron_immune_df, './Outputs/neuro_healthy_cellphonedb_nature//neuro_immune.csv')
  
  immune_to_neuron_merged_df <- base::merge(immune_neuron_b_df, immune_neuron_a_df, all=T)
  x <- immune_to_neuron_merged_df
  x[, 11:ncol(x)][is.na(x[, 11:ncol(x)])] <- 1
  write.table(x, paste0('./Outputs/neuro_',s,'_cellphonedb_nature/filtered_pvalues.txt'),
              sep='\t', row.names = F)
  
}