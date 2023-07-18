get_cytoband_coverage <- function(Grange) {
  #' @param Grange a grange with seqnames and granges.
  #' @returns returns a GRange with cytogenetic locations including
  #' percentage of overlap and a dataframe containing information on how much of
  #' the total arm was affected.
  
  #############
  # Libraries #
  #############
  require(tidyverse)
  require(plyranges)
  
  ########
  # Data #
  ########
  # Uncomment me, this is to download cytogenetic coordinates
  # library(biovizBase)
  # hg38IdeogramCyto <- getIdeogram("hg38", cytobands = TRUE)
  # save(hg38IdeogramCyto, file = 'rda/cna_analysis/hg38IdeogramCyto.rda')
  
  # loading cytogentic location
  load('rda/cna_analysis/hg38IdeogramCyto.rda')
  
  ###########
  # Wrangle #
  ###########
  # Convert both GRanges to similar karyogram style
  seqlevelsStyle(Grange) <- "NCBI"
  seqlevelsStyle(hg38IdeogramCyto) <- "NCBI"
  
  Grange_list_cytoband <- list()
  Grange_list_pq_affected <- list()
  
  # return empty Grange and empty data.frame if Grange is empty
  if (length(Grange) == 0) {
    # return empty Grange
    return(list(GRanges(
      seqnames = character(0),
      IRanges(start = integer(0), end = integer(0)),
      cytobands_cov =  character(0)
    ),
    data.frame(seqnames = character(), 
               start_of_arm = numeric(), 
               end_of_arm = numeric(), 
               width_of_arm = numeric(), 
               strand = character(), 
               chromosome = character(), 
               arm = character(), 
               percentage_arm_affected = numeric(), 
               stringsAsFactors = FALSE)
    ))
  }
  
  for (row in 1:length(Grange)) {
    for_range <- Grange[row]
    
    ##########
    # Part 1 #
    ##########
    ## First find specific cytoband locations
    # finding overlaps
    hits <- findOverlaps(for_range, hg38IdeogramCyto)
    
    # computing percentage overlap
    overlaps <-
      pintersect(for_range[queryHits(hits)], hg38IdeogramCyto[subjectHits(hits)])
    percentOverlap <-
      width(overlaps) / width(hg38IdeogramCyto[subjectHits(hits)])
    
    # extracting cytobands
    cytobands <-
      paste0(seqnames(hg38IdeogramCyto[subjectHits(hits)]), hg38IdeogramCyto[subjectHits(hits)]$name)
    # replace 'chr'
    cytobands <-
      str_replace(cytobands, pattern = 'chr', replacement = '')
    
    cytobands_cov <-
      list(paste0(cytobands, ';', round(percentOverlap, 2) * 100, '%'))
    
    for_range$cytobands_cov <- cytobands_cov
    
    Grange_list_cytoband[[row]] <- for_range
    
    ##########
    # Part 2 #
    ##########
    # make dataframe that shows how much of an arm is affected 
    
    # Grange with p and q IRanges
    p_q <- hg38IdeogramCyto
    
    p_q$arm <- str_extract(p_q$name, 'p|q')
    
    p_q <- subset(p_q, seqnames %in% c(seq(1:22), 'X', 'Y'))
    
    p_q <- p_q %>% group_by(seqnames, arm) %>% reduce_ranges()
    
    # finding overlaps
    hits_pq <- findOverlaps(for_range, p_q)
    
    # computing percentage overlap
    overlaps <-
      pintersect(for_range[queryHits(hits_pq)], p_q[subjectHits(hits_pq)])
    percentOverlap <-
      width(overlaps) / width(p_q[subjectHits(hits_pq)])
    
    p_q_percentage <- as.data.frame(p_q[subjectHits(hits_pq)])
    
    p_q_percentage$percentage_arm_affected <- percentOverlap * 100
    
    colnames(p_q_percentage) <-
      c(
        'seqnames',
        'start_of_arm',
        'end_of_arm',
        'width_of_arm',
        'strand',
        'chromosome',
        'arm',
        'percentage_arm_affected'
      )
    
    
    
    
    Grange_list_pq_affected[[row]] <- p_q_percentage
    
    
  }
  
  Grange <- bind_ranges(Grange_list_cytoband)
  pq_arm_affect <- bind_rows(Grange_list_pq_affected)
  
  return(list(Grange, pq_arm_affect))
  
  ### END of function ###
  
}


# Grange = common_cnas$deletions$hub005_prerad_a.a
#
#
# test <- hg38IdeogramCyto
# seqlevelsStyle(test) <- "NCBI"
#
# test$arm <- str_extract(test$name, 'p|q')
#
# test <- subset(test, seqnames %in% c(seq(1:22), 'X','Y'))
#
#
# test %>% group_by(seqnames, arm) %>% reduce_ranges()
