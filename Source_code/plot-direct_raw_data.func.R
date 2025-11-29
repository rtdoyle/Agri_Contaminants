plot_raw_data <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### for pea-specific traits 
  if (traits %in% c("wet_pod_weight",
                     "dry_pod_weight",
                    "per_pod_weight",
                     "infected_peas_perc",
                     "NDVI_mean",
                     "pods",
                     "flowers")){
    
    df <- df %>%
      filter(species == "pea") %>%
      droplevels(.)
  }
  
  ### means
  df.means <- df %>% 
    group_by(species, amend_spike, amend, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  controls <- df %>% 
    filter(amend == "control") %>%
    group_by(species, amend_spike, amend, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  
  ### specify species names
  species_names <- c(
    lettuce = 'Lettuce',
    pea = 'Pea',
    radish = 'Radish'
  )
  
  ### plots
  p <- ggplot(df.means,
              aes(x = amend_spike, 
                  y = mean)) +
    geom_bar(aes(fill = amend_spike),
             stat = "identity") +
    geom_beeswarm(data = df,
                  aes(y = get(traits))) +
    geom_hline(data = controls,
               aes(yintercept = mean),
               linetype = 2) +
    scale_fill_manual(name = "Amendment- \n Spiking level",
      values = c("gray",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500"),
      guide = guide_legend(nrow = 2)) +
    facet_wrap(~species, scales = "free_y",
               ncol = 3,
               labeller = 
                 labeller(species = species_names)) +
    labs(y = paste0(traits), 
         x = "Amendment type + contaminant level") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),  # removes tick marks
      axis.title.y = element_text(size = 14),
      strip.text = element_blank(),
      legend.position = "none",
      legend.key.size = unit(0.5, "cm"),   # Shrinks the size of legend keys
      legend.text = element_text(size = 8), # Shrinks the text size
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
    
  ### return
  return(p)
}