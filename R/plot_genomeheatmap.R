genomeheatmap <- function(selected.files, 
                          path, 
                          classes_daan = NULL, 
                          class.col = NULL, 
                          dendogram = F,
                          daan_colours = T,
                          name) {
  # This code was adapted from the AneuFinder package
  # plots and saves genome heatmap of selected.files
  # @param 
  # selected.files: vector of location of selected files.
  # output path: path were plots should be saved. 
  # name: how you want to name the plot
  
  # get platename
  platename <- str_extract(selected.files[1],  'KRA-0\\d+')
  
  # AneuFinder functions: 
  startTimedMessage <- function(...) {
    
    x <- paste0(..., collapse='')
    message(x, appendLF=FALSE)
    ptm <- proc.time()
    return(ptm)
    
  }
  
  transCoord <- function (gr) 
  {
    cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
    cum.seqlengths.0 <- c(0, cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(gr)
    gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    return(gr)
  }
  
  stopTimedMessage <- function(ptm) {
    
    time <- proc.time() - ptm
    message(" ", round(time[3],2), "s")
    
  }
  
  initializeStates <- function (states) 
  {
    somy.states <- grep("somy", states, value = TRUE)
    somy.numbers <- as.integer(sapply(strsplit(somy.states, "-somy"), 
                                      "[[", 1))
    names(somy.numbers) <- somy.states
    if ("zero-inflation" %in% states) {
      multiplicity <- c(`zero-inflation` = 0, somy.numbers)
    }
    else {
      multiplicity <- somy.numbers
    }
    levels.distributions <- c("delta", "dgeom", "dnbinom", "dbinom")
    distributions <- rep(NA, length(states))
    names(distributions) <- states
    distributions[states == "zero-inflation"] <- "delta"
    distributions[states == "0-somy"] <- "dgeom"
    distributions[(states != "zero-inflation") & (states != "0-somy")] <- "dnbinom"
    states <- factor(states, levels = states)
    distributions <- factor(distributions, levels = levels.distributions)
    l <- list(states = states, distributions = distributions, 
              multiplicity = multiplicity)
    return(l)
  }
  
  #  colours
  if (daan_colours) {
    stateColors <- function(states = c("zero-inflation", paste0(0:10, "-somy"),
                                       "total")) {
      
      state.colors <- c(`zero-inflation` = "#1d4661", `0-somy` = "#1d4661",
                        `1-somy` = "#3787BA", `2-somy` = "#95B8C5",
                        `3-somy` = "#F0ECEB", `4-somy` = "#D7A290", `5-somy` = "#BF583B",
                        `6-somy` = "#8D1128", `7-somy` = "#3C0912", `8-somy` = "black",
                        total = "black")
      
      states.with.color <- intersect(states, names(state.colors))
      cols <- rep("black", length(states))
      names(cols) <- states
      cols[states.with.color] <- state.colors[states.with.color]
      return(cols)
    }
    
  }

  
  
  # heatmapgenomewide adapted form AneuFinder package
  heatmapGenomewide_daan <- function (hmms, ylabels = NULL, classes, reorder.by.class = TRUE, 
                                      classes.color, file = NULL, cluster = TRUE, plot.breakpoints = FALSE, 
                                      hotspots = NULL, exclude.regions = NULL) 
  {
    if (!is.null(ylabels)) {
      if (length(ylabels) != length(hmms)) {
        stop("length(ylabels) must equal length(hmms)")
      }
    }
    if (!is.null(classes)) {
      if (length(classes) != length(hmms)) {
        stop("length(classes) must equal length(hmms)")
      }
    }
    if (length(classes.color) != length(unique(classes))) {
      stop("'classes.color' must have the same length as unique(classes)")
    }
    if (is.null(names(classes.color))) {
      names(classes.color) <- unique(classes)
    }
    if (!setequal(names(classes.color), unique(classes))) {
      stop("The names of 'classes.color' must be equal to the unique elements in 'classes'")
    }
    if (length(hmms) == 1 & cluster == TRUE) {
      cluster <- FALSE
      warning("Cannot do clustering because only one object was given.")
    }
    hmms <- loadFromFiles(hmms, check.class = c("aneuHMM", "aneuBiHMM"))
    class.data <- data.frame(ID = sapply(hmms, "[[", "ID"))
    class.data$ID <- factor(class.data$ID, levels = class.data$ID)
    if (is.null(ylabels)) {
      class.data$ylabel <- as.character(class.data$ID)
    }
    else {
      class.data$ylabel <- as.character(ylabels)
    }
    class.data$class <- classes
    mapping <- class.data$ylabel
    names(mapping) <- class.data$ID
    if (reorder.by.class) {
      cl <- clusterHMMs(hmms, cluster = cluster, classes = classes, 
                        exclude.regions = exclude.regions)
    }
    else {
      cl <- clusterHMMs(hmms, cluster = cluster, exclude.regions = exclude.regions)
    }
    hmms <- hmms[cl$IDorder]
    class.data <- class.data[cl$IDorder, ]
    class.data$ID <- factor(class.data$ID, levels = class.data$ID)
    segments.list <- GRangesList()
    for (i1 in 1:length(hmms)) {
      hmm <- hmms[[i1]]
      if (is.null(hmm$segments)) {
        segments.list[[hmm$ID]] <- GRanges()
      }
      else {
        segments.list[[hmm$ID]] <- hmm$segments
      }
    }
    if (plot.breakpoints) {
      breakpoints <- GRangesList()
      for (i1 in 1:length(hmms)) {
        hmm <- hmms[[i1]]
        if (is.null(hmm$breakpoints)) {
          breakpoints[[hmm$ID]] <- GRanges()
        }
        else {
          breakpoints[[hmm$ID]] <- hmm$breakpoints
        }
      }
      if (length(breakpoints) == 0) {
        plot.breakpoints <- FALSE
      }
    }
    ptm <- startTimedMessage("Transforming coordinates ...")
    segments.list <- endoapply(segments.list, transCoord)
    if (plot.breakpoints) {
      breakpoints <- endoapply(breakpoints, transCoord)
    }
    stopTimedMessage(ptm)
    ptm <- startTimedMessage("Making the plot ...")
    df <- list()
    for (i1 in 1:length(segments.list)) {
      df[[length(df) + 1]] <- data.frame(start = segments.list[[i1]]$start.genome, 
                                         end = segments.list[[i1]]$end.genome, seqnames = seqnames(segments.list[[i1]]), 
                                         ID = names(segments.list)[i1], state = segments.list[[i1]]$state)
    }
    df <- do.call(rbind, df)
    df$ID <- factor(df$ID, levels = levels(class.data$ID))
    df$ylabel <- mapping[as.character(df$ID)]
    if (plot.breakpoints) {
      df.breakpoints <- list()
      for (i1 in 1:length(breakpoints)) {
        if (length(breakpoints[[i1]]) > 0) {
          df.breakpoints[[length(df.breakpoints) + 1]] <- data.frame(start = breakpoints[[i1]]$start.genome, 
                                                                     end = breakpoints[[i1]]$end.genome, seqnames = seqnames(breakpoints[[i1]]), 
                                                                     ID = names(segments.list)[i1], mid = (breakpoints[[i1]]$start.genome + 
                                                                                                             breakpoints[[i1]]$end.genome)/2)
        }
        else {
          df.breakpoints[[length(df.breakpoints) + 1]] <- data.frame(start = numeric(), 
                                                                     end = numeric(), seqnames = character(), ID = character(), 
                                                                     mid = numeric())
        }
      }
      df.breakpoints <- do.call(rbind, df.breakpoints)
      df.breakpoints$ID <- factor(df.breakpoints$ID, levels = levels(class.data$ID))
      df.breakpoints$ylabel <- mapping[as.character(df.breakpoints$ID)]
    }
    cum.seqlengths <- cumsum(as.numeric(seqlengths(segments.list[[1]])))
    names(cum.seqlengths) <- seqlevels(segments.list[[1]])
    cum.seqlengths.0 <- c(0, cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(segments.list[[1]])
    label.pos <- round(cum.seqlengths.0 + 0.5 * seqlengths(segments.list[[1]]))
    df.chroms <- data.frame(y = c(0, cum.seqlengths), x = 1, 
                            xend = length(segments.list))
    pltlist <- list()
    widths <- vector()
    df$state <- factor(df$state, levels = names(sort(initializeStates(levels(df$state))$multiplicity)))
    df$x <- as.numeric(df$ID)
    ggplt <- ggplot(df) + geom_linerange(aes_string(ymin = "start", 
                                                    ymax = "end", x = "x", col = "state"), size = 5) + scale_y_continuous(breaks = label.pos, 
                                                                                                                          labels = names(label.pos))
    

    # ggplt <- ggplt + scale_x_continuous(name = "",
    #                      breaks = 1:length(unique(df$ylabel)),
    #                      labels = unique(df$ylabel))
    
    ggplt <- ggplt + scale_color_manual(values = stateColors(levels(df$state))) # adding custom colours
    # adjusintg x axis
    ggplt <- ggplt + theme(panel.background = element_blank(), 
                           axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                           axis.line = element_blank(), axis.title.x = element_blank())
    
    ggplt <- ggplt + geom_segment(aes_string(x = "x", xend = "xend",
                                             y = "y", yend = "y"), data = df.chroms, col = "grey13")
    
    
    
    # adjusting y axis
    ggplt <- ggplt + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                           axis.line = element_blank(), axis.title.x = element_blank())
    
    # removing legend
    ggplt <- ggplt + theme(legend.position="none")
    
    ggplt <- ggplt + coord_flip()
    
    # removing all axis names
    ggplt <- ggplt + ylab("") + xlab("")
    
    # add numeber of cells sequenced text
    ggplt <- ggplt + annotate("text", 
                              label = paste(length(hmms), 'cells'), 
                              x = length(hmms)/2, 
                              y = 3100000000,
                              angle = 270,
                              size = 50)
    
    # decreasing plot margin
    ggplt <- ggplt + theme(plot.margin = unit(c(0,0,0,0), "cm"))
    
    
    if (plot.breakpoints) {
      df.breakpoints$x <- as.numeric(df.breakpoints$ID)
      ggplt <- ggplt + geom_linerange(data = df.breakpoints, 
                                      mapping = aes_string(x = "x", ymin = "start", ymax = "end"), 
                                      size = 2) + ylab("") + geom_point(data = df.breakpoints, 
                                                                        mapping = aes_string(x = "x", y = "mid"))
    }
    if (!is.null(hotspots)) {
      if (length(hotspots) > 0) {
        df.hot <- as.data.frame(transCoord(hotspots))
        df.hot$xmin <- 0
        df.hot$xmax <- length(class.data$ID) + 1
        ggplt <- ggplt + geom_rect(data = df.hot, mapping = aes_string(xmin = "xmin", 
                                                                       xmax = "xmax", ymin = "start.genome", ymax = "end.genome", 
                                                                       alpha = "num.events"), fill = "hotpink4") + scale_alpha_continuous(name = "breakpoints", 
                                                                                                                                          range = c(0.4, 0.8))
      }
    }
    
    
    width.heatmap <- sum(as.numeric(seqlengths(hmms[[1]]$bins)))/3e+09 * 
      150
    height <- max(length(hmms) * 0.5, 2)
    pltlist[["heatmap"]] <- ggplt
    widths["heatmap"] <- width.heatmap
    # adding class colors
    if (!is.null(classes)) {
      width.classes <- 5
      class.data$x <- as.numeric(class.data$ID)
      ggclass <- ggplot(class.data) + geom_linerange(aes_string(ymin = 0, 
                                                                ymax = 1, x = "x", col = "class"), size = 5) + guides(col = FALSE) + 
        xlab("")
      ggclass <- ggclass + theme(panel.background = element_blank(), 
                                 axis.ticks = element_blank(), axis.text = element_blank(), 
                                 axis.line = element_blank(), axis.title.x = element_blank())
      ggclass <- ggclass + coord_flip()
      
      # decreasing plot margin
      ggclass <- ggclass + theme(plot.margin = unit(c(0,0,0,0), "cm"))
      
      if (!is.null(classes.color)) {
        ggclass <- ggclass + scale_color_manual(breaks = names(classes.color), 
                                                values = classes.color)
      }
      pltlist[["classbar"]] <- ggclass
      widths["classbar"] <- width.classes
    }
    # adding dendogram
    if (!is.null(cl$hclust) & dendogram) {
      dhc <- stats::as.dendrogram(cl$hclust)
      ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
      ggdndr <- ggplot(ddata$segments) + geom_segment(aes_string(x = "x", 
                                                                 xend = "xend", y = "y", yend = "yend")) + scale_y_reverse()
      ggdndr <- ggdndr + coord_flip()
      ggdndr <- ggdndr + theme(panel.background = element_blank(), 
                               axis.ticks = element_blank(), axis.text = element_blank(), 
                               axis.line = element_blank(), axis.title = element_blank())
      width.dendro <- 20
      pltlist[["dendro"]] <- ggdndr
      widths["dendro"] <- width.dendro
    }
    # alligning ggpllt with dendogram and classes
    cowplt <- cowplot::plot_grid(plotlist = rev(pltlist), align = "h", 
                                 ncol = length(pltlist), rel_widths = rev(widths))
    
    
    stopTimedMessage(ptm)
    if (!is.null(file)) {
      ptm <- startTimedMessage("Plotting to file ", file, " ...")
      ggsave(file, cowplt, width = sum(widths), height = height, 
             units = "cm", limitsize = FALSE)
      stopTimedMessage(ptm)
    }
    else {
      return(cowplt)
    }
  }
  
  #make heatmap and safe
  suppressWarnings(heatmapGenomewide_daan(selected.files, classes = classes_daan, classes.color = class.col,
                         file = paste0(path, '/', name, '_', platename, '.pdf')))
  
  return(print(paste0(name, ' done.')))
  
}