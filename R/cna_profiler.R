cna_profiler <- function(df) {
  #' @param df is a datafame containing all the copy number calls of all the files.
  #'                       Columns: 01) seqnames
  #'                                02) start
  #'                                03) end
  #'                                04) width
  #'                                05) strand
  #'                                06) [copy.number cell 1]
  #'                                07) [copy.number cell 2], etc.
  #' @return Outputs a dataframe that states the ratio (normalized to total cells) of
  #' small, medium, large and whole chromosome CNAs

  require(tidyverse)
  require(GWASTools)
  
  # defining some parameters
  # number of cells 
  n_cells <- ncol(df) - 5
  
  # extract centromere positions
  data(centromeres.hg38)
  
  centromeres.hg38 <- subset(centromeres.hg38, chrom != 'Y')
  names(centromeres.hg38) <- c('seqnames', 'start', 'end')
  
  
  # extract total length chromosomes
  total_length <- df %>% group_by(seqnames) %>% summarize(sum(width))
  
  # construct dataframe containing chromosome name, centromere position (start and end) and total length chromosome
  chrom_centro_total <- inner_join(total_length, centromeres.hg38)
  
  # compute arm lengths
  chrom_centro_total$p_arm <- chrom_centro_total$start + (chrom_centro_total$end-chrom_centro_total$start)
  chrom_centro_total$q_arm <- chrom_centro_total$`sum(width)` -  chrom_centro_total$p_arm
  
  ##### !!!!!!!! ######
  # Check if this is still correct when using new assemblies !
  
  # now all p_arms except for chromosome 1 are correct (for some reason the p-arm is the larger arm here). 
  # I am going to swithc it. 
  p_arm_1 <- chrom_centro_total$q_arm[1]
  q_arm_1 <- chrom_centro_total$p_arm[1]
  
  chrom_centro_total[1,5] <- p_arm_1
  chrom_centro_total[1,6] <- q_arm_1
  
  
  # function to extract CNA for each chromosome, run this function for each chromosome
  cna_extract <- function(chromosome) {
    # @parama: chromosome is the chromosome (i.e., 1, 2,3, X, etc) 
    df_chrom <- subset(df, seqnames == chromosome)
    
    # run this function for each column in df_chrom
    cna_extract_cell <- function(column_i) {
      # @ param column i is the integer of the column in df_chrom 
      # from which to extract CNAs
      
      # extract rows that switch CNA value (for example from 2 to 3)
      CNA_change <- which(c(FALSE, tail(df_chrom[,column_i],-1) != head(df_chrom[,column_i],-1)))
      
      # extract indices
      CNA_change_indices <- unique(sort(c(1,CNA_change, (CNA_change-1), nrow(df_chrom))))
      
      # construct dataframe containing end and start of each CNA
      CNA_change_df <- df_chrom[CNA_change_indices,c(1,2,3,4,5, column_i)]
      names(CNA_change_df) <- c('seqnames', 'start', 'end', 'width', 'strand', 'cna')
      
      start_list <- list()
      end_list <- list()
      
      i <- 1
      
      while(i<=nrow(CNA_change_df)){
        
        try(if(CNA_change_df[i+1, 6] == CNA_change_df[i, 6]){
          start_list[[(length(start_list) + 1)]] <- CNA_change_df[i,2]
          end_list[[(length(end_list) + 1)]] <- CNA_change_df[i+1, 3]
          i <- i+2
        })
        
        if(CNA_change_df[i+1, 6] != CNA_change_df[i, 6] || is.na(CNA_change_df[i+1, 6])){
          start_list[[(length(start_list) + 1)]] <- CNA_change_df[i,2]
          end_list[[(length(end_list) + 1)]] <- CNA_change_df[i, 3]
          i <- i+1
        }
        
      }
      
      # constructing a dataframe
      cna_dataframe <- data.frame(seqnames = chromosome, 
                                  start = unlist(start_list), 
                                  end = unlist(end_list))
      
      # remove NA
      cna_dataframe <- na.omit(cna_dataframe)
      
      # add column CNA
      cna_dataframe$cna <- (CNA_change_df %>% filter(cna!=lag(cna) | is.na(lag(cna))))$cna
      
      # remove columsn with cna == 2 (this is no CNA but normal diploid)
      cna_dataframe <- subset(cna_dataframe, cna != 2)
      
      # add total number of CNAs in this cell as colunn num_cna
      try(cna_dataframe$num_cna <- length(cna_dataframe$cna), silent = T)
      
      # add a width column of the CNA
      cna_dataframe$width <- cna_dataframe$end - cna_dataframe$start + 1
      
      # extract chromosome specifics
      chrom_df <- subset(chrom_centro_total, seqnames == chromosome)
      
      # in the next of code, we are going to compute the fraction of SCNA
      # Important: using this measure, a fraction of 1 means that a region the size of one arm
      # was affected. A fraction of 2 means the whole chromosome is affected. Note that these CNA
      # can cross centromeres, so we are unable to easily assign if it is the p-arm or q-arm that is 
      # affected. 
      cna_dataframe$fraction <- ifelse(cna_dataframe$end < chrom_df$p_arm, (cna_dataframe$width / chrom_df$p_arm), 
                                       ifelse(cna_dataframe$start > chrom_df$p_arm, cna_dataframe$width / chrom_df$q_arm,
                                              (((chrom_df$p_arm - cna_dataframe$start) / chrom_df$p_arm) + ((cna_dataframe$end - chrom_df$p_arm) / chrom_df$q_arm))))
      
      # adding cellname
      try(cna_dataframe$cell_id <- colnames(df_chrom)[column_i], silent = T)
      
      return(cna_dataframe)
      
    }
    
    # running for each cell. 
    cna_list <- lapply(6:ncol(df), cna_extract_cell)
    
    # constructing single dataframe 
    cna_df <- bind_rows(cna_list, .id = "cell_num_scar")
    
    # removing cell_num scar
    cna_df <- cna_df[, 2:ncol(cna_df)]
    
    return(cna_df)
    
  }
  
  cna_df_per_chrom_list <- list()
  
  for(chrom in unique(df$seqnames)) {
    cna_df_per_chrom_list[[chrom]] <- cna_extract(chrom)
  }
  
  # constructing single dataframe 
  cna_df_per_chrom_df <- bind_rows(cna_df_per_chrom_list, .id = "chrom_scar")
  
  # removing cell_num scar
  cna_df_per_chrom_df <- cna_df_per_chrom_df[, 2:ncol(cna_df_per_chrom_df)]
}