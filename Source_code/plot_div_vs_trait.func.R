plot_div_vs_trait.func <- function(metrics, df){
  
  ## print metric
  print(metrics)
  
  # Convert metric to symbol for tidy evaluation
  metric <- df[[metrics]]
  
  
  ## trait labs
  ### rename traits
  trait_names <- c(
    resp_norm_per_plant_weight = 'Per plant weight (g)',
    resp_norm_dry_shoot_weight = 'Dry shoot weight (g)',
    resp_norm_dry_total_weight = 'Dry total weight (g)',
    resp_norm_dry_root_weight = 'Dry root weight (g)',
    resp_norm_wet_shoot_weight = 'Wet shoot weight (g)',
    resp_norm_wet_root_weight = 'Wet root weight (g)',
    resp_norm_shoot_root_ratio = 'Shoot to root ratio',
    resp_norm_shoot_moisture = 'Shoot moisture (%)',
    resp_norm_germination_perc = 'Germination (%)',
    resp_norm_survival_perc.corr = 'Survival (%)',
    resp_norm_pods = "Pods (no.)",
    resp_norm_plant_height4 = "Plant height (cm)",
    resp_norm_chlorophyll_content4_mean = "NDVI",
    resp_norm_leaf_number2 = "Leaves (no.)"
  )
  
  ## custom labels
  time_names <- c(
    TP1 = "1 wk PP",
    TP2 = "@ harvest"
  )
  
  ## plot
  p <-
    ggplot(data = df,
           aes(y = metric, x = units)) +
    geom_point(size = 3, alpha = 0.8,
               aes(colour = treat.rn)) + 
    geom_vline(aes(xintercept = 0), linetype = 2) +
    geom_smooth(method = "lm", colour = "black",
                linetype = 2, se = FALSE, formula = 'y ~ x') +
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    guides(color = "none", shape = "none") +
    # ggtitle(paste0(metrics)) + 
    labs(x = NULL, y = NULL) +
    facet_wrap(~trait, scales = "free", ncol = 2,
               labeller = 
                 labeller(trait = trait_names)) +  
    theme_bw() + 
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  ### outputs
  return(p)
  
} 