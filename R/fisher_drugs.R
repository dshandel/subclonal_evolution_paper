compute_fisher_drugs <- function(top5, all, pattern) {
  #' @param top5 is a subset of dataframe all than contains info on drugs that are in 
  #' the top5 most correlated to radiation resistancy drugs
  #' @param all see top5
  #' @param pattern the pattern to look for in column moa of top5 and rest
  #' @return fisher's exact output that test if the pattern is more often found in 
  #' top 5 or not. 
  
  # extract frequency of pattern in top5
  pattern_top5 <- nrow(top5_df[grepl(pattern = pattern , top5_df$moa),])
  
  # extract frequency of pattern in all
  pattern_all <- nrow(all[grepl(pattern = pattern , all$moa),])
  # extract frequency of pattern in non top5
  pattern_rest <- pattern_all - pattern_top5
  
  # extract total number of drugs in top5
  total_top5 <- nrow(top5)
  
  # extract total number of drugs in rest
  total_rest <- nrow(all) - total_top5
  
  # make contigency table
  ct <- matrix(c(pattern_top5,pattern_rest, total_top5 - pattern_top5, total_rest - pattern_rest ), nrow =2 , ncol = 2)
  colnames(ct) <- c(pattern, paste0('non-',pattern))
  rownames(ct) <- c('top5', 'rest')
  
  # computing statistics
  fisher <- fisher.test(ct)
  
  # return 
  fisher$p.value
  
}
