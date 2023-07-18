extract_other_cells <- function(kra_name, list_of_subclones) {
  i <- which(names(man_select_files_edivisive) == kra_name)
  all <-  str_extract(man_select_files_edivisive[[i]],  '_\\d+\\.') %>% str_remove('_') %>% str_remove('\\.') 
  rest <- all[which(!all %in% list_of_subclones)]
  return(rest)
}


extract_square <- function(list_files1, list_files2, square, x, y) {
  # extracts all cells that lie within a square of a pca plot defined by 
  # user. For example square c(0.1, 0.0) will extract all cells with
  # pc1 lower or higher that 0.1 and pc2 lower or higher than 0.0. User has to rearange functon
  # if he wants x lower or higher than 0.1
  cells <- c(list_files1, list_files2)
  
  df12 <- plot_pca(cells,
                   colorBy = classes, 
                   PC1=1, 
                   PC2=2,
                   plot = F)
  
  colnames(df12) <- c('PC1', 'PC2')
  
  
  # extracting dataframe with correct square
  if(x == '<' & y == '<') {
    df12 <- subset(df12, PC1 < square[1] & PC2 < square[2])
  }
  
  if(x == '<' & y == '>') {
    df12 <- subset(df12, PC1 < square[1] & PC2 > square[2])
  }
  
  if(x == '>' & y == '<') {
    df12 <- subset(df12, PC1 > square[1] & PC2 < square[2])
  }
  
  if(x == '>' & y == '>') {
    df12 <- subset(df12, PC1 > square[1] & PC2 > square[2])
  }
  
  # numbers is: 
  cell_id <-  str_extract(rownames(df12),  '_\\d+\\.') %>% str_remove('_') %>% str_remove('\\.')  
  
  return(cell_id)
  
}

l_g <- function(x, kra_num, list_cells) {
  select <- grepl(paste0(kra_num,'_',x,"\\."), list_cells)
  i <- which(select)
  return(list_cells[i])
  
  
}


subclone_classes_for_genomeheatmap <- function(list_of_subclones, cline) {
  #' @param 
  #' list of subclones is the full list of rda files that have been assigned a subclone_classes_for_genomeheatmap
  #' cline is the line (hub005, hub183 etc.)
  #' rad is either prerad or rad
  #' @return 
  #' returns a vector with that assigns a subclone
  
  extract_subclones <- list_of_subclones[grepl(names(list_of_subclones), pattern = paste0(cline))]
  
  size_subclones <- lapply(extract_subclones, length)
  
  classes <- rep(names(size_subclones), size_subclones)
  
  df <- data.frame(path = unlist(extract_subclones),
                   class = classes)
  
  return(df)
  
  
}
