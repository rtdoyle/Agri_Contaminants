plot_emms <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### means
  df.means <- df %>% 
    group_by(amend_spike, amend, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  controls <- df %>% 
    filter(amend == "control") %>%
    group_by(amend_spike, amend, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  
  ## filtered for traits
  emm.f <- emm.s %>%
    filter(!amend %in% c("saline","no_treat") &
             trait == traits) %>%
    droplevels(.)
  
  ### plots
  p <- ggplot(df.means,
              aes(x = amend_spike, 
                  y = mean)) +
    geom_bar(aes(fill = amend_spike),
             stat = "identity") +
    geom_pointrange(data = emm.f,
                  aes(y = emmean,
                      ymin = LCL,
                      ymax = UCL,
                      group = amend_spike),
                  position = position_jitter(width = 0.3,
                                             height = 0)) +
    geom_hline(data = controls,
               aes(yintercept = mean),
               linetype = 2) +
    scale_fill_manual(values = c("gray",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500")) +
    labs(y = paste0(traits), 
         x = "Amendment - Contaminant level") +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16),
      strip.text = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  
  ### return
  return(p)
}