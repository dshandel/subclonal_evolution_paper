find_specific_cna_sc <-
  function(list_resistant,
           list_sensitive,
           perc_cutoff_within_subclone_resistant,
           perc_cutoff_within_subclone_sensitive,
           per_cutoff_across_subclone_resistant,
           per_cutoff_across_subclone_sensitive,
           return_cnas_per_subclone = F) {
    #' @param list_resistant a list subclone names which itself contains a list of paths of paths to aneuHMM .Rda files. 
    #' @param list_sensitive a list subclone names which itself contains a list of paths of paths to aneuHMM .Rda files. 
    #' @param perc_cutoff_within_subclone_resistant user defined percentage cutoff to use (e.g. 0.6 means a certain CNA is shared by 60% of the cells
    #' within a subclone).
    #' @param perc_cutoff_within_subclone_senstive same as in perc_cutoff_within_subclone_resistant but for the sensitive group
    #' @param per_cutoff_across_subclone_resistant user defined percentage cutoff to use (e.g. 0.7 means a certain CNA is shared in 70% of all subclones)
    #' @param return_cnas_per_subclone boolean. If TRUE will return the shared by perc_cutoff_within_subclone CNAs for each list of each item within 
    #' list_resistant. Note list_sensitive should be set to NULL for proper results!!
    #' @param per_cutoff_across_subclone_sensitive same as resistant but for the sensitive group
    #' @example perc_cutoff_within_subclone = 0.6 and per_cutoff_across_subclone = 0.7 means the algorithm will find CNAs that are shared 
    #' by 60% of single cells within each subclone, and by 70% of all subclones. Important, a cutoff of 50% is the same as computing the median.
    #' @return a list of GRanges showing unique to list_resistant amplifications and deletions and unique to list_sensitive
    #' amplifications and deletions or, if return_cnas_per_subclone see param return_cnas_per_subclone above.
    #' @details Note that this algorithm is NOT sensitive to the number of single cells within each subclone. 
    
    #############
    # Libraries #
    #############
    require(plyranges)
    require(AneuFinder)
    
    ###########
    # Part 1. #
    ###########
    print('Computing (this may take some minutes)...')
    
    resistant_sc <-
      lapply(list_resistant,
             loadFromFiles,
             check.class = c("aneuHMM", "aneuBiHMM"))
    sensitive_sc <-
      lapply(list_sensitive,
             loadFromFiles,
             check.class = c("aneuHMM", "aneuBiHMM"))
    
    # find specific CNAs needs a list of GRanges where each GRange has columns
    # seqnames, ranges, strand (not necessary) and cn.
    wrangle_sc <- function(list_lff) {
      # segments is the reduced bins GRange (so all equal CN bins put together)
      wrangle_within <- function(list_lff_within) {
        segments <- list_lff_within$segments
        # ugly way of adding a cn column
        segments$cn <- segments$copy.number
        
        return(segments)
        
        # END of wrange_within function
        
      }
      within_subcl <- lapply(list_lff, wrangle_within)
      
      # END OF wrangle_sc function
      
    }

    list_resistant <- lapply(resistant_sc, wrangle_sc)
    list_sensitive <- lapply(sensitive_sc, wrangle_sc)
    
    print('Extracted segment info from files')
    Sys.sleep(1) 
    
    ###########
    # Part 2. #
    ###########
    # The Granges are now still arranged in bins. We are first going to
    # combine bins with similar copy number.
    print('Computing (this may take some minutes)...')
    combine_indels <- function(subclone_consensus) {
      #' @param subclone_consensus_list is a GRanges with seqnames, ranges, strand and cn
      #' @returns list of Granges where each Granges has only CNA of cn2 or cn>2
      
      combine_within <- function(subclone_consensus_within) {
        # subset subclone_consensus in cn <2, cn = 2 and cn > 2
        cn_less_2 <-
          subclone_consensus_within[subclone_consensus_within$cn < 2,]
        cn_more_2 <-
          subclone_consensus_within[subclone_consensus_within$cn > 2,]
        
        # reduce each GRanges (joins consecutive ranges together)
        cn_less_2 <- reduce(cn_less_2)
        cn_more_2 <- reduce(cn_more_2)
        
        # reduce deletes metadata.
        # adding metadata
        cn_less_2$cn <- rep('<2', length(cn_less_2))
        
        cn_more_2$cn <- rep('<2', length(cn_more_2))
        
        # make list and rename items
        GR_list <- GRangesList(cn_less_2, cn_more_2)
        names(GR_list) <- c('<2', '>2')
        
        # Return #
        return(GR_list)
      }
      
      combine_subcl <- lapply(subclone_consensus, combine_within)
      
    }
    
    # run an apply function here on each Granges in list_resistant versus list_sensitive
    resistant_split <- lapply(list_resistant, combine_indels)
    sensitive_split <- lapply(list_sensitive, combine_indels)
    
    # splits deletions together within each subclone
    grouper_del <- function(gr) {
      grouper_within <- function(gr_within) {
        return(gr_within$'<2')
      }
      
      within_grouper <- lapply(gr, grouper_within)
    }
    
    # splits amplifciations together within each subclone
    grouper_ampl <- function(gr) {
      grouper_within <- function(gr_within) {
        return(gr_within$'>2')
      }
      
      within_grouper <- lapply(gr, grouper_within)
    }
    
    # run apply on each subclone
    resistant_split_del <- lapply(resistant_split, grouper_del)
    resistant_split_ampl <- lapply(resistant_split, grouper_ampl)
    
    sensitive_split_del <- lapply(sensitive_split, grouper_del)
    sensitive_split_ampl <- lapply(sensitive_split, grouper_ampl)
    
    ###########
    # Part 3. #
    ###########
    # Part 3A #
    ###########
    # Finding common CNAs within each subclone
    shared_cnas <-
      function(list_GRanges,
               perc_cutoff_within_subclone) {
        #' @param list_GRanges is a list of GRanges with seqnames, ranges, strand and cn
        #' @param perc_cutoff_within_subclone user defined percentage cutoff to use (e.g. 0.6)
        #' @returns Granges where each Granges is shared by => perc_cutoff_within_subclone
        
        # finding overlaps in x% of regions: wait for response community:
        # https://stackoverflow.com/questions/74900085/find-ranges-that-are-shared-by-80-or-more-of-10-granges-objects
        # https://support.bioconductor.org/p/9148540/
        
        shared_cnas_within <- function(list_GRanges_within) {
          n <- length(list_GRanges_within) # number of range sets
          shared_cna <-
            bind_ranges(list_GRanges_within, .id = "origin") %>%
            compute_coverage() %>%
            mutate(fraction_cov = score / n) %>%
            filter(fraction_cov >= perc_cutoff_within_subclone) %>%
            reduce_ranges()
          
          return(shared_cna)
        }
        
        lapply(list_GRanges, shared_cnas_within)
      }
    
    resistant_shared_del <-
      shared_cnas(resistant_split_del, perc_cutoff_within_subclone = perc_cutoff_within_subclone_resistant)
    resistant_shared_ampl <-
      shared_cnas(resistant_split_ampl, perc_cutoff_within_subclone = perc_cutoff_within_subclone_resistant)
    sensitive_shared_del <-
      shared_cnas(sensitive_split_del, perc_cutoff_within_subclone = perc_cutoff_within_subclone_sensitive)
    sensitive_shared_ampl <-
      shared_cnas(sensitive_split_ampl, perc_cutoff_within_subclone = perc_cutoff_within_subclone_sensitive)
    
    # END of function if return_cnas_per_subclone == T
    if (return_cnas_per_subclone == T) {
      cna_list <- list (deletions = c(resistant_shared_del, sensitive_shared_del), 
                        amplifications = c(resistant_shared_ampl, sensitive_shared_del))
      
      print(
        paste0(
          'Extracted CNAs that are shared by ',
          as.numeric(perc_cutoff_within_subclone_resistant) * 100,
          '% of single cells within each resistant subclone, and by ',
          as.numeric(perc_cutoff_within_subclone_sensitive) * 100,
          '% of single cells within each sensitive subclone'
        )
      )
      ##########
      # Unload #
      ##########
      suppressWarnings(invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
                                        detach,
                                        character.only = TRUE, unload = TRUE)))
      
      return(cna_list)
    }
    
    # Part 3B #
    ###########
    # Finding common CNAs across subclones
    shared_cnas_across_subclones <-
      function(list_GRanges,
               per_cutoff_across_subclone) {
        #' @param list_GRanges is a list of GRanges with seqnames, ranges, strand and cn
        #' @param perc_cutoff user defined percentage cutoff to use (e.g. 0.6)
        #' @returns Granges where each Granges is shared by => perc_cutoff_within_subclone
        
        # finding overlaps in x% of regions: wait for response community:
        # https://stackoverflow.com/questions/74900085/find-ranges-that-are-shared-by-80-or-more-of-10-granges-objects
        # https://support.bioconductor.org/p/9148540/
        
        n <- length(list_GRanges) # number of range sets
        shared_cna <-
          bind_ranges(list_GRanges, .id = "origin") %>%
          compute_coverage() %>%
          mutate(fraction_cov = score / n) %>%
          filter(fraction_cov >= per_cutoff_across_subclone) %>%
          reduce_ranges()
        
        return(shared_cna)
        
      }
    
    resistant_shared_del_across <-
      shared_cnas_across_subclones(resistant_shared_del, per_cutoff_across_subclone = per_cutoff_across_subclone_resistant)
    resistant_shared_ampl_across <-
      shared_cnas_across_subclones(resistant_shared_ampl, per_cutoff_across_subclone = per_cutoff_across_subclone_resistant)
    sensitive_shared_del_across <-
      shared_cnas_across_subclones(sensitive_shared_del, per_cutoff_across_subclone = per_cutoff_across_subclone_sensitive)
    sensitive_shared_ampl_across <-
      shared_cnas_across_subclones(sensitive_shared_ampl, per_cutoff_across_subclone = per_cutoff_across_subclone_sensitive)
    
    print(
      paste0(
        'Extracted CNAs that are shared by ',
        as.numeric(perc_cutoff_within_subclone_resistant) * 100, ' and ', as.numeric(perc_cutoff_within_subclone_sensitive) * 100,
        '% of single cells within resistant and sensitive subclones, respectively, and that are shared by ',
        as.numeric(per_cutoff_across_subclone_resistant) * 100, ' and ', as.numeric(per_cutoff_across_subclone_sensitive) * 100, 
        '% across resistant and sensitive subclones, respectively.'
      )
    )
    
    Sys.sleep(1) 
    print('Computing (few seconds)...')
    Sys.sleep(2) 
    
    ###########
    # Part 4. #
    ###########
    # Now that we have shared CNAs across subclones with the resistant and sensitive group, we need to extract
    # only those CNAs that are unique to the resistant group.
    
    # unique to resistant lines deletions:
    unique_resistant_del <-
      setdiff(resistant_shared_del_across, sensitive_shared_del_across)
    # unique to resistant lines amplifications
    unique_resistant_ampl <-
      setdiff(resistant_shared_ampl_across,
              sensitive_shared_ampl_across)
    
    # unique to sensitive lines deletions:
    unique_sensitive_del <-
      setdiff(sensitive_shared_del_across, resistant_shared_del_across)
    # unique to sensitive lines amplifications:
    unique_sensitive_ampl <-
      setdiff(sensitive_shared_ampl_across,
              resistant_shared_ampl_across)
    
    print(
      'Identified amplifications and deletions that are unique to resistant or sensitive subclones'
    )
    
    ##########
    # Return #
    ##########
    returner <-
      list(
        unique_resistant_del,
        unique_resistant_ampl,
        unique_sensitive_del,
        unique_sensitive_ampl
      )
    names(returner) <-
      c(
        'unique_resistant_del',
        'unique_resistant_ampl',
        'unique_sensitive_del',
        'unique_sensitive_ampl'
      )
    
    ##########
    # Unload #
    ##########
    suppressWarnings(invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
                     detach,
                     character.only = TRUE, unload = TRUE)))
    
    return(returner)
  }


