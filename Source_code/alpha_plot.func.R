plot_alpha_div.alt.func <- function(combs, times, metrics, df){
  
  ## print time-measure comb
  print(combs)
  
  ## set to characters
  metric <- as.character(metrics)
  time <- as.character(times)
  
  ## subset to specific timepoints
  df.f <- prune_samples(sample_data(df)$timepoint %in% time, df)
  
  df.f <- prune_samples(sample_data(df.f)$species != "no_plant", df.f)
  
  ### specify species names
  species_names <- c(
    L = 'Lettuce',
    R = 'Radish',
    P = 'Pea'
  )
  
  ## plot
  
  p1 <- df.f %>%
    plot_richness(x="treat.rn", 
                  color = "treat.rn", 
                  shape = "Species",
                  measures = metric) + 
    geom_point(size = 3) + 
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    guides(color = "none", shape = "none") +
    scale_shape_manual(values = c(16,17,15),
                       labels = c("Lettuce", "Pea", "Radish")) + 
    ggtitle(NULL) + 
    labs(x = NULL, y = NULL) +
    facet_wrap(~Species, 
               labeller = labeller("Species" = species_names)) +  
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(p1)
}