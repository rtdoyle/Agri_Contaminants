plot_div.RR_species.vs.contam <- function(metrics, df){
  
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
  
  ### rename labels
  time_labs <- c(
    "TP1" = '1 wk PP',
    "TP2" = 'Harvest')
  
  ##### Contrasts (contamination level) #####
  p <- ggplot(data = df.f,
              aes(x = factor(spikeFac), 
                  y = emmean, 
                  colour = amend)) +
    geom_line(aes(group = amend),
              position = position_dodge(0.5)) +
    geom_pointrange(aes(ymin = lower.CL,
                        ymax = upper.CL),
                    position = position_dodge(0.5)) +
    ### add in significance (0.05 <= p < 0.1) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-2" &
                               pval_linear >= 0.05 & pval_linear < 0.1),
      aes(x = factor(spikeFac), y = emmean),
      label = "+",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (p < 0.05) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-2" &
                               pval_linear < 0.05),
      aes(x = factor(spikeFac), y = emmean),
      label = "*",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (0.05 <= p < 0.1) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-0" &
                               pval_quadratic >= 0.05 & pval_quadratic < 0.1),
      aes(x = factor(spikeFac), y = emmean),
      label = "+",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (p < 0.05) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-0" &
                               pval_quadratic < 0.05),
      aes(x = factor(spikeFac), y = emmean),
      label = "*",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ## additional formatting
    scale_colour_manual(values = c(
      "blue",
      "#B06500")) +
    scale_x_discrete(labels = c("0","1","2")) +
    facet_wrap(.~species, scales = "fixed",
               labeller = labeller(species = species_labs,
                                   time = time_labs)) +
    labs(x = NULL, y = NULL) +
    guides(colour = "none") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
  
  ## save
  return(p)
  
}