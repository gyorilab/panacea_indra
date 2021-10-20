df_to_viper_regulon <- function(df)
{
  names(df) <- c("feature","pathway","sign")
  df <- df[complete.cases(df),]
  
  pathway_regulon <- list(0)
  i <- 1
  for(pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1,length(features))
    names(pathway_feature_list) <- c("tfmode","likelihood")
    
    pathway_regulon[[i]] <- pathway_feature_list
    i <- i+1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}

run_viper <- function(eset, prior, outfile){
  #Now we estimate the TF activities using viper
  TF_activities <- as.data.frame(viper(eset = eset, regulon = prior, 
                                       minsize = 10, adaptive.size = F, 
                                       eset.filter = F, pleiotropy = T))
  TF_activities$TF <- row.names(TF_activities)
  
  #that's just to make the dataframe pretty
  TF_activities <- TF_activities[,c(2,1)]
  names(TF_activities) <- c("ID","NES")
  
  TF_activities$NES_abs <- abs(TF_activities$NES)
  TF_activities$NES_abs[which(TF_activities$NES < 0)] <- TF_activities$NES_abs[which(TF_activities$NES < 0)]*(-1)
  TF_activities$NES <- NULL
  TF_activities <- TF_activities %>% arrange(desc(TF_activities$NES_abs))
  
  write.csv(TF_activities, paste0('./output/', outfile, '_TF.csv'),
            row.names = F)
  return(TF_activities)
}
