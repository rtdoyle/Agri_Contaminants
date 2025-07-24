plot_raw_data_func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### means
  df.means <- df %>% 
    group_by(species, contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  controls <- df %>% 
    filter(contam == "control") %>%
    group_by(species, contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  
  ### plots
  p <- ggplot(df.means,
              aes(x = contam_spike, 
                  y = mean)) +
    geom_bar(aes(fill = contam_spike),
             stat = "identity") +
    geom_beeswarm(data = df,
                  aes(y = get(traits))
    ) +
    geom_hline(data = controls,
               aes(yintercept = mean),
               linetype = 2) +
    facet_wrap(~species, scales = "free_y") +
    scale_fill_manual(values = c("gray",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue")) +
    labs(y = paste0(traits), 
         x = "Contaminant source") +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16),
      strip.text = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  
  ### return
  return(p)
}