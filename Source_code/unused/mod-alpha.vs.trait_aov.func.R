mod_alpha.vs.trait_aov <- function(combs, metrics, traits, df){
  
  ## print metric
  print(combs)
  
## filter for each trait
  df.f <- df %>%
    filter(trait == traits) %>%
    droplevels(.)
  
  # dynamically assign the column to a new variable
  df.f$metric_value <- df.f[[metrics]]  
  
  ## specify non pea-specific traits
  if (!traits %in% c("direct.dry_pod_weight.corr",
                     "direct.per_pod_weight",
                     "direct.infected_peas_perc.corr",
                     "direct.NDVI_mean",
                     "direct.pods",
                     "direct.flowers")){
  
  ### linear model
  lm <- lm(units ~ log(metric_value)*species*timepoint, data = df.f)
  }
  
  else if (traits %in% c("direct.dry_pod_weight.corr",
                         "direct.per_pod_weight",
                          "direct.infected_peas_perc.corr",
                          "direct.NDVI_mean",
                          "direct.pods",
                          "direct.flowers")){
    
    ### linear model
    lm <- lm(units ~ log(metric_value)*timepoint, data = df.f)
  }
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, 
                        "resfits_plots/mod-alpha.vs.trait_aov/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  
  
  ### res vs fit plot
  png(paste0(output_dir2, combs, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lm)
  layout(1)
  dev.off()
  
  ###anova
  aov <- as.data.frame(Anova(lm, type = 3))
  
  ### processing
  aov$metric <- paste0(metrics)
  aov$term <- rownames(aov)
  aov$trait <- paste0(traits)
  
  ### model coefficients
  lm_sum <- summary(lm)
  lm_coeff <- as.data.frame(lm_sum$coefficients)
  lm_coeff$metric <- paste0(metrics)
  lm_coeff$trait <- paste0(traits)
  lm_coeff$term <- rownames(lm_coeff)
  
  return(list(aov, lm_coeff))
}

