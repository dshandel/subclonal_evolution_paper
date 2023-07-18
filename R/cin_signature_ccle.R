run_cin_ccle <- function(cnv, d) { 
  #' @param cnv is segmental data to input in CINSignatureQuantification
  #' @param d is a dataframe with sample metadata that is joined with cnv
  #' @return outputs CIN signature values.
  
  #############
  # Libraries #
  #############
  require(CINSignatureQuantification)
  require(tidyverse)
  
  ###########
  # Wrangle #
  ###########
  wrangler_cin_signature_ccle <- function(df) {
    #' @param  df: A dataframe containing all the copy number calls of all the files.
    #' @return a dataframe in the format https://github.com/markowetzlab/CINSignatureQuantification
    
    # remove unwanted column names
    keep_columns <-
      c('Chromosome',
        'Start',
        'End',
        'Segment_Mean',
        'DepMap_ID')
    
    new_df <- subset(df, select = keep_columns)
    
    # reorder columns
    col_order <- c('Chromosome',
                   'Start',
                   'End',
                   'Segment_Mean',
                   'DepMap_ID')
    
    new_df <- new_df[, col_order]
    
    # define new column names
    column_names <-
      c('chromosome',
        'start',
        'end',
        'segVal',
        'sample')
    
    colnames(new_df) <- column_names
    
    return(new_df)
    
  }
  
  # Subset cnv to only include samples where radiation sensitivity is known
  cnv <-
    cnv[cnv$DepMap_ID %in% d$DepMap_ID,]
  
  # Wrangle
  cnv <-
    wrangler_cin_signature_ccle(cnv)
  
  # The CIN algorithm (run with quantifyCNSignatures) does not support Y chromosome, we have to get rid of those.
  cnv <-
    cnv[cnv$chromosome != 'Y', ]
  
  #######
  # Run #
  #######
  # Running the CIN scores
  cin_scores <-
    quantifyCNSignatures(cnv, build = 'hg38')
  
  # Extracting the values.
  cin_values <- getActivities(cin_scores)

  ###########
  # Wrangle #
  ###########
  cin_df <- as.data.frame(cin_values)
  
  cin_df$DepMap_ID <- rownames(cin_df)
  
  cin_df <-
    gather(cin_df, cx_signature, cx_value, CX1:CX17, factor_key = T)
  
  # add some metadata
  cin_df <- inner_join(cin_df, d)

  ##########
  # Return #
  ##########
  return(cin_df)
  
}





