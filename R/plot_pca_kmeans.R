draw_pca <-
  function(list_files1,
           list_files2,
           size,
           legend_position,
           baseline_kra,
           recurrence_kra,
           ssdna004 = F) {
    cells <- c(list_files1, list_files2)
    
    df <- plot_pca(
      cells,
      colorBy = classes,
      PC1 = 1,
      PC2 = 2,
      plot = F
    )
    
    df$class <- str_extract(rownames(df),  'KRA-0\\d+')
    num <- str_extract(rownames(df),  '_\\d+\\.')
    num <- str_remove_all(num, pattern = '\\.')
    df$label <- paste0(df$class, num)
    pc1 <- colnames(df)[1]
    pc2 <- colnames(df)[2]
    
    colnames(df) <- c('PC1', 'PC2', 'class', 'label')
    
    col_vals <- ifelse(ssdna004 == T, list(c("#D69C4E", "#999999")), list(c("#999999", "#D69C4E")))
    label_vals <- ifelse(ssdna004 == T, list(c('Recurrence', 'Baseline')), list(c('Baseline', 'Recurrence')))
    
    p <- ggplot(data = df, aes(label = label)) +
      geom_point(size = size, aes(x = PC1,
                                  y = PC2,
                                  col = class)) +
      # scale_color_manual(values = c("#FFDB6D", "#00AFBB")) +
      scale_color_manual(
        name = '',
        values = col_vals[[1]],
        labels = label_vals[[1]]
      ) +
      
      # or consider '#9e1303' for second cycle of recurrence
      
      theme_cowplot() +
      theme(legend.position = legend_position,
            text = element_text(size = 17),
            axis.text = element_text(size = 15)) +
      
      labs(x = pc1, y = pc2)
    p
    
  
    
    # ggplotly(p)
    
  }


k_cluster <- function(list_files1,
                      list_files2,
                      cluster_n,
                      cols,
                      labels,
                      legend_position,
                      normalize,
                      return_elbow = F) {
  
  require(factoextra)
  
  cells <- c(list_files1, list_files2)
  
  df12 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 1,
    PC2 = 2,
    plot = F
  )

  
  df34 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 3,
    PC2 = 4,
    plot = F
  )
  
  df56 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 5,
    PC2 = 6,
    plot = F
  )
  
  df78 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 7,
    PC2 = 8,
    plot = F
  )
  
  df910 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 9,
    PC2 = 10,
    plot = F
  )
  
  
  df <- cbind(df12, df34, df56, df78, df910)
  
  
  if (normalize) {
    # extract variance explained for each pca
    pca_var <-
      as.numeric(str_match(colnames(df), "\\(\\s*(.*?)\\s*%")[, 2])
    
    for (i in 1:length(pca_var)) {
      # calculate normalized pca values
      norm_pca <- df[, i] / 100 * pca_var[i]
      #repopulate the dataframe
      df[, i] <- norm_pca
    }
    
    
    
    
    
  }
  
  #check best k for k-means clustering method
  if (return_elbow) {
    return(fviz_nbclust(df, kmeans, method = "silhouette") +
             labs(subtitle = "Elbow method"))
    
  }
  
  # #
  # fviz_nbclust(df, kmeans, method = "silhouette")+
  #   labs(subtitle = "Silhouette method")
  #
  # fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  #   labs(subtitle = "Gap statistic method")
  
  # compute k-means clustering
  set.seed(1)
  res.km <- kmeans(df, cluster_n, nstart = 1000)
  
  # change colnames so it works for our plot
  df12$class <- str_extract(rownames(df12),  'KRA-0\\d+')
  num <- str_extract(rownames(df),  '_\\d+\\.')
  num <- str_remove_all(num, pattern = '\\.')
  df12$label <- paste0(df12$class, num)
  pc1 <- colnames(df12)[1]
  pc2 <- colnames(df12)[2]
  
  colnames(df12) <- c('PC1', 'PC2', 'class', 'label')
  
  # plot nicely
  p <- ggplot(data = df12, aes(label = label)) +
    geom_point(alpha = 1.0,
               aes(
                 x = PC1,
                 y = PC2,
                 col = as.factor(res.km$cluster)
               )) +
    
    
    theme_cowplot() +
    
    labs(x = pc1, y = pc2) +
    scale_color_manual(labels = labels,
                       values = cols,
                       name = "Clone") +
    theme(legend.position = legend_position,
          text = element_text(size = 17),
          axis.text = element_text(size = 15)) 
  
  return(p)
  
  #ggplotly(p)
  
}

draw_pca_double <-
  function(list_files1,
           list_files2,
           list_files3,
           size,
           legend_position) {
    cells <- c(list_files1, list_files2, list_files3)
    
    
    df <- plot_pca(
      cells,
      colorBy = classes,
      PC1 = 1,
      PC2 = 2,
      plot = F
    )
    
    df$class <- str_extract(rownames(df),  'KRA-0\\d+')
    num <- str_extract(rownames(df),  '_\\d+\\.')
    num <- str_remove_all(num, pattern = '\\.')
    df$label <- paste0(df$class, num)
    pc1 <- colnames(df)[1]
    pc2 <- colnames(df)[2]
    
    colnames(df) <- c('PC1', 'PC2', 'class', 'label')
    
    
    p <- ggplot(data = df, aes(label = label)) +
      geom_point(size = size, aes(x = PC1,
                                  y = PC2,
                                  col = class)) +
      # scale_color_manual(values = c("#FFDB6D", "#00AFBB")) +
      scale_color_manual(
        name = '',
        values = c("#999999", "#D69C4E", '#d65c4e'),
        labels = c('Baseline', 'Recurrence Cycle 1', 'Recurrence Cycle 2')
      ) +
      
      # or consider '#9e1303' for second cycle of recurrence
      
      theme_cowplot() +
      theme(legend.position = legend_position,
            text = element_text(size = 17),
            axis.text = element_text(size = 15)) +
      
      
      labs(x = pc1, y = pc2)
    return(p)
    
    # ggplotly(p)
    
  }


k_cluster_double <- function(list_files1,
                             list_files2,
                             list_files3,
                             
                             cluster_n,
                             cols,
                             labels,
                             legend_position,
                             normalize,
                             return_elbow = F) {
  
  require(factoextra)
  
  cells <- c(list_files1, list_files2, list_files3)
  
  
  df12 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 1,
    PC2 = 2,
    plot = F
  )
  
  
  df34 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 3,
    PC2 = 4,
    plot = F
  )
  
  df56 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 5,
    PC2 = 6,
    plot = F
  )
  
  df78 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 7,
    PC2 = 8,
    plot = F
  )
  
  df910 <- plot_pca(
    cells,
    colorBy = classes,
    PC1 = 9,
    PC2 = 10,
    plot = F
  )
  
  
  df <- cbind(df12, df34, df56, df78, df910)
  
  
  if (normalize) {
    # extract variance explained for each pca
    pca_var <-
      as.numeric(str_match(colnames(df), "\\(\\s*(.*?)\\s*%")[, 2])
    
    for (i in 1:length(pca_var)) {
      # calculate normalized pca values
      norm_pca <- df[, i] / 100 * pca_var[i]
      #repopulate the dataframe
      df[, i] <- norm_pca
    }
    
    
    
    
    
  }
  
  #check best k for k-means clustering method
  if (return_elbow) {
    return(fviz_nbclust(df, kmeans, method = "silhouette") +
             labs(subtitle = "Silhouette method"))
    
  }
  
  # #
  # fviz_nbclust(df, kmeans, method = "silhouette")+
  #   labs(subtitle = "Silhouette method")
  #
  # fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  #   labs(subtitle = "Gap statistic method")
  
  # compute k-means clustering
  set.seed(1)
  res.km <- kmeans(df, cluster_n, nstart = 1000)
  
  # change colnames so it works for our plot
  df12$class <- str_extract(rownames(df12),  'KRA-0\\d+')
  num <- str_extract(rownames(df),  '_\\d+\\.')
  num <- str_remove_all(num, pattern = '\\.')
  df12$label <- paste0(df12$class, num)
  pc1 <- colnames(df12)[1]
  pc2 <- colnames(df12)[2]
  
  colnames(df12) <- c('PC1', 'PC2', 'class', 'label')
  
  # plot nicely
  p <- ggplot(data = df12, aes(label = label)) +
    geom_point(alpha = 1.0,
               aes(
                 x = PC1,
                 y = PC2,
                 col = as.factor(res.km$cluster)
               )) +
    
    
    theme_cowplot() +
    
    labs(x = pc1, y = pc2) +
    scale_color_manual(labels = labels,
                       values = cols,
                       name = "Clone") +
    theme(legend.position = legend_position,
          text = element_text(size = 17),
          axis.text = element_text(size = 15)) 
  
  p
  
  #ggplotly(p)
  
}
