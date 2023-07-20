get_genes <- function(GRange) { 
  #' @param Grange a grange with seqnames and granges. 
  #' @returns returns df with all genes encoded within the GRange regions. 
  
  require(GenomicRanges)
  require(biomaRt)
  
  # This code should be run only once and is therefore commented. 
  # code from: https://www.biostars.org/p/311199/
  
  # Set up an gene annotation template to use
  # mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  # mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  # genes <- getBM(attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"), mart=mart)
  # genes <- genes[genes[,1]!="" & genes[,2] %in% c(1:22,"X","Y"),]
  # xidx <- which(genes[,2]=="X")
  # yidx <- which(genes[,2]=="Y")
  # genes[xidx, 2] <- 23
  # genes[yidx, 2] <- 24
  # genes[,2] <- sapply(genes[,2],as.integer)
  # genes <- genes[order(genes[,3]),]
  # genes <- genes[order(genes[,2]),]
  # save(genes, file = 'rda/cna_analysis/genes.rda')
  # colnames(genes) <- c("GeneSymbol","Chr","Start","End")
  # genes_GR <- makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)
  # save(genes_GR, file = 'rda/cna_analysis/genes_GR.rda')
  
  # loading genelist and genomic location
  load("../cna_analysis/rda/cna_analysis/genes_GR.rda", envir = .GlobalEnv)
  load('../cna_analysis/rda/cna_analysis/genes.rda', envir = .GlobalEnv )
  
  # Convert both GRanges to similar karyogram style
  seqlevelsStyle(GRange) <- "UCSC"
  GRange <- renameSeqlevels(GRange, paste0('chr', c(1:23)))
  
  seqlevelsStyle(genes_GR) <- "UCSC"
  
  # Finding hits with genes_GR and GRanges
  hits <- findOverlaps(genes_GR, GRange, type="within")
  
  # constructing
  df <- cbind(as.data.frame(GRange)[subjectHits(hits),],genes[queryHits(hits),])
  
  return(df)
  
}



