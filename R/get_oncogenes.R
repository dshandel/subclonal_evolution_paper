get_oncogenes <- function(GRange) {
  #' @param Grange a grange with seqnames and granges. 
  #' @returns returns df with all genes encoded within the GRange regions. 
  source('R/get_genes.R')
  
  # first finding all genes within GRange
  all_genes <- get_genes(GRange) 
  
  # change colnames so dataframes match
  colnames(all_genes)[which(names(all_genes) == 'hgnc_symbol')] <- "Gene Symbol"   
  
  # downloading oncogenes
  require(readr)
  oncogenes <- read_csv("../cna_analysis/data/cna_analysis/Census_allMon May 30 12_43_01 2022.csv")
  
  
  # subset df_genes so to only have oncogenes
  require(tidyverse)
  oncogenes <- inner_join(oncogenes, all_genes)
  
  return(oncogenes)
  
}
