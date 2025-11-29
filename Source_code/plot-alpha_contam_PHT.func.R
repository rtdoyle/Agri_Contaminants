plot_div_species.vs.contam <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  
  ##Separate into two columns
  df <- df %>%
    separate(metric, into = c("metric", "amend2"), sep = "_")
  
  
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
                        y = log(response), 
                        colour = amend)) +
    geom_line(aes(group = amend),
              position = position_dodge(0.5)) +
    geom_pointrange(aes(ymin = log(lower.CL),
                        ymax = log(upper.CL)),
                    position = position_dodge(0.5)) +
    ### add in significance (0.05 <= p < 0.1) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-2" &
                                 pval_linear >= 0.05 & pval_linear < 0.1),
      aes(x = factor(spikeFac), y = log(response)),
      label = "+",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (p < 0.05) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-2" &
                               pval_linear < 0.05),
      aes(x = factor(spikeFac), y = log(response)),
      label = "*",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (0.05 <= p < 0.1) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-0" &
                               pval_quadratic >= 0.05 & pval_quadratic < 0.1),
      aes(x = factor(spikeFac), y = log(response)),
      label = "+",
      size = 6,
      position = position_dodge(0.5),
      vjust = -0.5
    ) +
    ### add in significance (p < 0.05) 
    geom_text(
      data = df.f %>% filter(spikeFac == "SL-0" &
                               pval_quadratic < 0.05),
      aes(x = factor(spikeFac), y = log(response)),
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
    facet_grid(timepoint~species, scales = "fixed",
               labeller = labeller(species = species_labs,
                                   timepoint = time_labs)) +
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