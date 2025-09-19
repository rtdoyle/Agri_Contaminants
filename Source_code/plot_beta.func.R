plot_beta.func <- function(combs, times, metrics, df){
  
  ## print time-measure comb
  print(combs)
  
  ## set to characters
  metric <- as.character(metrics)
  time <- as.character(times)
  
  ## subset to specific timepoints
  df.f1 <- 
    prune_samples(sample_data(df)$species!= c("no_plant"), df)
  df.f <- prune_samples(sample_data(df.f1)$timepoint %in% time, df.f1)
  
  ## ordinate
  RA.ord <- ordinate(df.f, method = metric, 
                              distance = "bray")
  
  
  ## plot
  ### add in ggplot aesthetics:
  p <-  plot_ordination(df.f, RA.ord, axes = c(1,2))
  
  p$layers <- p$layers[-1]
  # To remove the extra jitter layer
  
  p <- p +
    geom_point(aes(colour = treat.rn, shape = Species), alpha = 0.7, size = 5) + 
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    guides(shape="none", colour="none") + 
    scale_shape_discrete(labels=c("Lettuce", "Pea", "Radish")) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
          axis.title.x = element_text(size=14), 
          axis.title.y = element_text(size=14), 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size= 12), 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size =10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if (metrics == "PCoA"){
  
  #Extract Eigenvalues
  eigenvalues <- RA.ord$values$Eigenvalues
  
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
                  size = 1, 
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
  ggsave(paste0("./16S_outputs/figs/scree_plot_", times, ".png"),
         width = 6, height = 5, units = "in",
         dpi = 400)
  }
  
  return(list(p))
}