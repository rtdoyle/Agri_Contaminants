plot_alpha_div <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  df_sum <- df %>%
    group_by(species, treat.rn) %>%
    summarize(count = n(),
      mean = mean(get(metrics), na.rm = TRUE),
      min = min(get(metrics), na.rm = TRUE),
      max = max(get(metrics), na.rm = TRUE))
  
  df_sum.C <- df %>%
    filter(treat.rn == "control-0") %>%
    droplevels(.) %>%
    group_by(species, treat.rn) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE))
  
  ### specify species names
  species_names <- c(
    L = 'Lettuce',
    R = 'Radish',
    P = 'Pea'
  )
  
  ### plots
  p <- ggplot(df_sum,
              aes(x = treat.rn, 
                  y = mean)) +
    geom_pointrange(aes(
      ymin = min,
      ymax = max,
      colour = treat.rn),
      size = 0.75) +
    geom_hline(aes(yintercept = 0),
               linetype = 1) +
    geom_hline(data = df_sum.C,
               aes(yintercept = mean),
               linetype = 2) +
    scale_colour_manual(name = "Amendment- \n Contaminant level",
                      values = c("gray",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500"),
                      guide = guide_legend(nrow = 4, ncol = 2)) +
    facet_wrap(~species, scales = "fixed",
               labeller = 
                 labeller(species = species_names)) +
    labs(y = paste0(metrics), 
         x = "Amendment type + Contaminant level") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),  # removes tick marks
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ### return
  return(p)
}