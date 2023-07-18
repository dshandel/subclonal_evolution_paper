
perform_fisher <- function(dataframe,
                           population1,
                           population2) {
  #' @param 
  #' dataframe is a df with columns paths and class with values 'hub015_rad_a', 'hub015_prerad_a' etc.
  #' population 1 and population two are the two pops you want to compare, for example prerad vs rad
  #' or prerad versus rad cycle 2
  #' @return
  #' returns the p value of a Fisher's exact test. If the p value is above 0.05, the distributions
  #' pre and post rad are the same. 
  #' 
  
  dataframe$rad <- ifelse(dataframe$class %in% population1, 'pop1',
                          ifelse(dataframe$class %in% population2, 'pop2',
                                 NA))
  
  org_id <- unique(str_extract(unique(dataframe$class), 'hub\\d{3}'))
  
  # remove na
  dataframe <- dataframe[which(!is.na(dataframe$rad)),]
  
  simple_fun <- function(string) {
    return(tail(string,1))
  }
  
  dataframe$clone <- unlist(lapply(stringr::str_split(dataframe$class, pattern = '_'), simple_fun))
  
  # in hub015 and hub005, we defined subclones. These are actually part of the same clone, so should be 
  # analysed together
  if (org_id %in% c('hub005', 'hub015')){
    dataframe$clone <- str_replace_all(dataframe$clone, pattern = 'a.a|a.b', replacement = 'a')
    
  }
  
  tab <- table(dataframe$clone, dataframe$rad)
  
  
  # df <- data.frame(cline <- org_id, p_val <- fisher.test(tab)$p.value)
  # 
  # return(df) 
  return(fisher.test(tab)$p.value)
  
}
