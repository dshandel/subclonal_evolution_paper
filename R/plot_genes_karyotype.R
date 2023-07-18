plot_genes_karyotype <-
  function(GRange,
           dist = -30,
           show_genes,
           show_chromosomes = c(
             "chr1",
             "chr2",
             "chr3",
             "chr4",
             "chr5",
             "chr6",
             "chr7",
             "chr8",
             "chr9",
             "chr10",
             "chr11",
             "chr12",
             "chr13",
             "chr14",
             "chr15",
             "chr16",
             "chr17",
             "chr18",
             "chr19",
             "chr20",
             "chr21",
             "chr22",
             "chrX"
           )) {
    #' @param GRange a GRange with seqnames, ranges, strand and hgnc_symbol
    #' @param dist How far away the genes should be plotted from the lines.
    #' Between -40 and -70 seems to do the trick most of the times.
    #' @param show_genes is a vector of gene symbols to show.
    #' @param show_chromosomes is a vector indicating which chromosomes to show.
    #' @return a karyotype with genes labeled
    #' @details. Throws an error when the GRange does not contain either
    #' ampflications or deletions. This can be ignored, the resulting plot
    #' is correct
    
    #############
    # Libraries #
    #############
    require(karyoploteR)
    require(stringr)
    
    ###########
    # Wrangle #
    ###########
    # making sure seqlevels match
    seqlevelsStyle(GRange) <- 'UCSC'
    
    kp <- plotKaryotype(
      "hg38",
      plot.type = 5,
      labels.plotter = NULL,
      main = "",
      cex = 4,
      chromosomes = show_chromosomes
    )
    
    #renaming seqlevels
    seqlevels_grange_change <-
      str_replace_all(seqlevels(GRange),
                      pattern = 'chr23',
                      replacement = 'chrX')
    
    GRange <- renameSeqlevels(GRange, value = seqlevels_grange_change)
    
    kpAddChromosomeNames(
      kp,
      chr.names = str_remove_all(show_chromosomes, 'chr'),
      cex = 1.3
    )
    
    # split into amplifications and deletions
    ampl <- GRange[GRange$cn == 'amplified']
    ampl <- ampl[ampl$hgnc_symbol %in% show_genes]
    
    del <- GRange[GRange$cn == 'deleted']
    del <- del[del$hgnc_symbol %in% show_genes]
    
    tryCatch(
      cor <-
        cor.test(df$AUC, df$dependency, use = "complete.obs"),
      error = function(e)
        NULL
    )
    
    tryCatch(
      kpPlotMarkers(
        kp,
        data = ampl,
        labels = ampl$hgnc_symbol,
        ignore.chromosome.ends = T,
        r0 = 1,
        r1 = 0.85,
        label.dist = 0.003,
        label.margin = dist,
        marker.parts = c(0.3, 0.1, 0.1),
        line.color = 'firebrick',
        label.color = 'firebrick'
      ),
      error = function(e)
        NULL
    )
    
    
    tryCatch(
      kpPlotMarkers(
        kp,
        data = del,
        labels = del$hgnc_symbol,
        ignore.chromosome.ends = T,
        r0 = 1,
        r1 = 0.4,
        label.dist = 0.003,
        cex = 8,
        label.margin = dist,
        marker.parts = c(0.3, 0.1, 0.1),
        line.color = 'steelblue',
        label.color = 'steelblue'
      ),
      error = function(e)
        NULL
    )
    
    
  }
