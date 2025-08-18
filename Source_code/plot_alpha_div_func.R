plot_alpha_div.func <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  df_sum <- df %>%
    filter(timepoint != "other") %>%
    group_by(species, treat.rn, tp_treat, species_tp, species_treat, timepoint) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE),
              min = min(get(metrics), na.rm = TRUE),
              max = max(get(metrics), na.rm = TRUE))
  
  df_sum.c <- df %>%
    filter(treat.rn == "control-0") %>%
    droplevels(.) %>%
    group_by(species, treat.rn, tp_treat, species_tp, species_treat, timepoint) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE))
  
    ### specify species names
  species_names <- c(
    L = 'Lettuce',
    R = 'Radish',
    P = 'Pea'
  )
  
  ## custom labels
  time_labs <- c(
    TP1 = "One week post-planting",
    TP2 = "At Harvest"
  )
  
  ### plots
  p <- ggplot(df_sum,
              aes(x = treat.rn, 
                  y = mean)) +
    geom_bar(aes(fill = treat.rn),
             stat = "identity") +
    geom_beeswarm(data = df,
                  aes(y = get(metrics))) +
    geom_hline(data = df_sum.c,
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
    facet_wrap(timepoint~species, scales = "free_y",
               ncol = 3,
               labeller = 
                 labeller(species = species_names,
                          timepoint = time_labs)) +
    labs(y = paste0(metrics), 
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