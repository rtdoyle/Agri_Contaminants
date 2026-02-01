plot_div_species.vs.amend <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  ## filter to metric
  df.f <- df %>%
    filter(metric == metrics) %>%
    droplevels(.)
  
  ### rename labels
  species_labs <- c(
    "L" = 'Lettuce',
    "R" = 'Radish',
    "P" = 'Pea')
  
  ##### contrasts (amendment type) #####
  p <- ggplot(data = df.f,
              aes(x = amend, 
                  y = eff, 
                  colour = amend)) +
    geom_pointrange(aes(ymin = LCL,
                        ymax = UCL),
                    position = position_dodge(0.5)) +
    ### add in significance (0.05 <= p < 0.1) 
    geom_text(
      data = df.f %>% filter(pval >= 0.05 & pval < 0.1),
      aes(x = amend, y = eff),
      label = "+",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (p < 0.05) 
    geom_text(
      data = df.f %>% filter(pval < 0.05),
      aes(x = amend, y = eff),
      label = "*",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.3
    ) +
    ## additional formatting
    geom_hline(aes(yintercept = 0),
               linetype = 2) +
    coord_flip() +
    scale_colour_manual(values = c(
      "blue",
      "#B06500")) +
    facet_grid(timepoint~species, scales = "fixed",
               labeller = labeller(species = species_labs)) +
    labs(x = NULL, y = NULL) +
    guides(colour = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      strip.text.y = element_blank(),
      strip.text.x = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  
  ## save
  return(p)
  
}