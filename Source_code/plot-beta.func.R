plot_beta <- function(combs, times, metrics, df){
  
  ## print time-measure comb
  print(combs)
  
  ## set to characters
  metric <- as.character(metrics)
  time <- as.character(times)
  
  ## subset to specific timepoints
  df.sub <- 
    prune_samples(sample_data(df)$species!= c("no_plant"), df)
  df.f <- prune_samples(sample_data(df.sub)$timepoint %in% times, df.sub)
  
  ## plot with Bray Curtis distances
  ## ordinate
  ord_data <- ordinate(df.f, method = metrics, 
                       distance = "bray")
  
  var_explained <- ord_data$values$Relative_eig * 100
  
  if (metrics == "PCoA"){
  
  # Extract ordination coordinates for samples
  ord_df <- plot_ordination(df.f, 
                            ord_data, 
                            axes = c(1, 2), justDF = TRUE)
  
  }
  
  else if (metrics == "NMDS"){
    
    
    # Extract sample scores (sites)
    ord_df <- as.data.frame(vegan::scores(ord_data, display = "sites"))
    ord_df$SampleID <- rownames(ord_df)
    
    # Merge with sample metadata
    ord_df <- merge(ord_df, as(sample_data(df.f), "data.frame"), 
                    by.x = "SampleID", by.y = "row.names")
    
  }
  
  ## ensure they are consistent
  names(ord_df)[names(ord_df) %in% c("NMDS1", "Axis.1")] <- "Axis.1"
  names(ord_df)[names(ord_df) %in% c("NMDS2", "Axis.2")] <- "Axis.2"
  
  ## plot
  
  ## formatting
  species_labs <- c(
    L = 'Lettuce',
    R = 'Radish',
    P = 'Pea'
  )
  
  p <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2, 
                          color = treat.rn)) +
    geom_point(aes(colour = treat.rn, 
                   shape = species),
               size = 4, 
               stroke = 2,
               alpha = 0.7) + 
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    scale_shape_manual(values = c(15, 16, 17),
                       labels = species_labs) +
    labs(x = paste0(metric, " Axis 1 [", round(var_explained[1], 1), "%]"),
         y = paste0(metric, " Axis 2 [", round(var_explained[2], 1), "%]"),
         shape = "Crop species",
         colour = "Amendment- \n contaminant level") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size= 12), 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =10, face = "bold"),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if (metrics == "PCoA"){
  
  #Extract Eigenvalues
  eigenvalues <- ord_data$values$Eigenvalues
  
  # Number of components
  n <- length(eigenvalues)
  
  # Calculate broken stick model
  broken_stick <- sapply(1:n, function(i) {
    sum(1 / (i:n)) / n * 100
  })
  
  # Create a dataframe for plotting
  scree_data <- data.frame(
    PrincipalCoordinate = 1:n,
    VarianceExplained = eigenvalues / sum(eigenvalues) * 100,
    BrokenStickModel = broken_stick
  )
  
  ggplot(scree_data, 
                  aes(x = PrincipalCoordinate)) + 
    geom_line(aes(y = VarianceExplained),
                  color = "blue", 
                  linewidth = 1, 
                  linetype = "solid") + 
    geom_line(aes(y = BrokenStickModel), 
              color = "red", 
              linewidth = 1, 
              linetype = "dashed") +
    geom_point(aes(y = VarianceExplained)) + 
    labs(title = paste0("Scree Plot for PCoA ", times), 
         x = "Principal Coordinate", 
         y = "Variance Explained (%)") + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  ## set output dir
  ggsave(paste0("./16S_analyses_files/figs/plot-scree_", times, ".png"),
         width = 6, height = 5, units = "in",
         dpi = 400)
  }
  
  return(p)
}