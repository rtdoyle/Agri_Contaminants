plot_raw_data.peas_func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  df.f <- df %>%
    filter(species == "pea") %>%
    droplevels(.)
  
  ### means
  df.means <- df.f %>% 
    group_by(species, contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  controls <- df.f %>% 
    filter(contam == "control") %>%
    group_by(species, contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  
  ### plots
  p <- ggplot(df.means,
              aes(x = contam_spike, 
                  y = mean)) +
    geom_bar(aes(fill = contam_spike),
             stat = "identity") +
    geom_beeswarm(data = df.f,
                  aes(y = get(traits))
    ) +
    geom_hline(data = controls,
               aes(yintercept = mean),
               linetype = 2) +
    # facet_wrap(~species, scales = "free_y") +
    scale_fill_manual(values = c("gray",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue")) +
    labs(y = paste0(traits), 
         x = "Treatment") +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  
  ### return
  return(p)
}