plot_raw_data_func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### for pea-specific traits 
  if (traits %in% c("wet_pod_weight",
                     "dry_pod_weight",
                     "infected_peas_perc",
                     "pods",
                     "flowers")){
    
    df <- df %>%
      filter(species == "pea") %>%
      droplevels(.)
  }
  
  ### means
  df.means <- df %>% 
    group_by(species, contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  controls <- df %>% 
    filter(contam == "control") %>%
    group_by(species, contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  
  ### specify species names
  species_names <- c(
    lettuce = 'Lettuce',
    pea = 'Pea',
    radish = 'Radish'
  )
  
  ### plots
  p <- ggplot(df.means,
              aes(x = contam_spike, 
                  y = mean)) +
    geom_bar(aes(fill = contam_spike),
             stat = "identity") +
    geom_beeswarm(data = df,
                  aes(y = get(traits))) +
    geom_hline(data = controls,
               aes(yintercept = mean),
               linetype = 2) +
    scale_fill_manual(name = "Amendment- \n Spiking Level",
      values = c("gray",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500"),
      guide = guide_legend(nrow = 4, ncol = 2)) +
    facet_wrap(~species, scales = "free_y",
               ncol = 3,
               labeller = 
                 labeller(species = species_names)) +
    labs(y = paste0(traits), 
         x = "Amendment type + Contaminant level") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),  # removes tick marks
      axis.title.y = element_text(size = 14),
      strip.text = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
    
  ### return
  return(p)
}