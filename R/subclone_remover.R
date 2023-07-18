remove_selection <- function(files, ns_to_remove) {
  #' @param files is a vector of paths of single cells .Rdata (AneuFinder output)
  #' @param ns_to_remove are the bam numbers to remove. Please make ns_to_remove of type character.
  #' @example if you want to remove this bam file: single.SCC−scKaryo−UMC−KRA−003_372.bam
  #' than ns_to_remove should contain c('372', ....etc)
  #' @return returns vector of paths to bam files with ns_to_remove removes
  
  # make file object
  file <- files
  # make empty list (this will contain the index of files that need to be removed)
  remove <- list()
  # get platename
  name <- str_extract(files[1],  'KRA-0\\d+')
  
  # extract index of files that need to be remove and put in a list (remove)
  for (i in ns_to_remove) {
    bool <- which(str_detect(files, paste0('UMC-', name, '_', i)))
    remove[[paste0('element', i)]] <- bool
    
  }
  
  # make list numeric and unlist.
  remove <- as.numeric(unlist(remove))
  
  # remove is an integer of vectors we want to get rid of.
  # construct a boolean.
  
  # first get all index of original files object
  total <- seq(1:length(files))
  
  # bool has the index that needs to be remove set to FALSE, and the rest to TRUE
  bool <- !(total %in% remove)
  
  # if we do this
  clean_file <- file[bool]
  
  # than all TRUEs in bool (integers that need to stay) will stay.
  
  return(clean_file)
  
}