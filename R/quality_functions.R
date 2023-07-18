quality_check <- function(inputdir, model) {
  #' @param inputdir is directory containing the standard aneufinder output for each plate
  #' model is which model to use for the quality check: either dnacopy, edivisive or HMM
  #' quality_check will make quality metrics for each plate
  #' quality_check assumes that rdaBaseDirectory contains folder with models for each plate
  #' each plate folder has to contain a MODEL folder and a method-'dnacopy, edivisive, or HMM' folder
  #' @example quality_check(inputdir = rdaBaseDirectory, model = 'edivisive')
  
  # go to folder with appropiate files
  rda_folder <- paste0(inputdir, '/MODELS', '/method-', model)
  
  # extract files
  rda_files <- list.files(rda_folder, full.names = TRUE)
  
  # run quality check on each file
  cl <-
    clusterByQuality(
      rda_files,
      measures = c(
        'spikiness',
        'num.segments',
        'entropy',
        'bhattacharyya',
        'sos'
      )
    )
  
  return(cl)
  
}

quality_select <- function(cl, spik, bhat) {
  #' selects all files that meets the defined quality requirements
  #' @param cl output of the clusterByQuality function in AneuFinder
  #' @param spik the spikiness treshhold. All clusters with spikiness above
  #' the treshold will be removed
  #' @param bhat the bhattacharrya score. All clusters with bhattacharrya below
  #' the treshold will be removed
  #' @return returns a vector of selected files paths
  
  # convert to df to allow easy wrangling
  cl_df <- as.data.frame(cl$parameters)
  
  # remove higher than spik
  spik <- subset(cl_df, spikiness < spik)
  
  # remove lower than bhat
  spik_bhat <- subset(spik,  bhattacharyya > bhat)
  
  cluster_n <- nrow(spik_bhat)
  
  # create vector of selected files
  selected.files <- unlist(cl$classification[0:cluster_n])
  
  return(selected.files)
  
}

check_quality <- function(cl) {
  return(cl$parameters)
}