plot_emms.func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### means
  df.means <- df %>% 
    group_by(contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  controls <- df %>% 
    filter(contam == "control") %>%
    group_by(contam_spike, contam, spikeFac) %>%
    summarize(mean = mean(get(traits), na.rm = TRUE))
  
  ## filtered for traits
  emm.f <- emm %>%
    filter(!contam %in% c("saline","no_treat") &
             trait == traits) %>%
    droplevels(.)
  
  ### plots
  p <- ggplot(df.means,
              aes(x = contam_spike, 
                  y = mean)) +
    geom_bar(aes(fill = contam_spike),
             stat = "identity") +
    geom_pointrange(data = emm.f,
                  aes(y = response,
                      ymin = lower.CL,
                      ymax = upper.CL),
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
         x = "Amendment type + Contaminant level") +
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