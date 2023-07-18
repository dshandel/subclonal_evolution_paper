# returns the AUC of the actual data.
AUC_fun <- function(df) {
  
  # defining variables to work with
  d <- df
  
  lines <- character()
  # for loop, checks if CF$cline is in vector lines, if so, do nothing, else append to vector. 
  for (line in d$cline) (
    ifelse(line %in% lines, NA, lines <- c(lines, line)))
  
  # separate into dataframes according to cline
  per_line <- split(d, f = d$cline, drop = TRUE)
  
  # define dose_vetor
  dose_vector <- unique(d$dose)
  
  
  # computes AUC for given predicted survival and dose given by dose_vector
  AUC_calc <- function(df) {
    # extract actual data
    rel_mean <- df$relative_mean
    
    # extracts number of experiemnts 
    n <- length(unique(d$Exp))
    #remove 2th and 3th number (are duplicates)
    rel_mean <- rel_mean[seq(1, length(rel_mean), n)]
    
    # compute maximal 
    max_AUC <- max(dose_vector)*1
    
    #compute AUC
    act_AUC <- AUC(x=dose_vector, y = rel_mean)
    
    # make dataframe
    datafr <- data.frame(max_AUC = max_AUC,
                         act_AUC = act_AUC, 
                         rel_AUC = act_AUC/max_AUC,
                         expcode = unique(df$expcode),
                         cline= unique(df$cline),
                         treatment = unique(df$treatment ) )
    
    
    return(datafr)
  }
  
  # Store in vector
  AUC_vector <- lapply(per_line, AUC_calc) 
  
  #merge dataframes
  df = do.call(rbind, AUC_vector)
  
  return(df) 
  
}