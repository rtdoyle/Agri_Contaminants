lin_reg.func <- function(traits, df){
  
  ## print traits
  print(paste0(traits))
  
  ## filter for each trait
  df.f <- df %>%
    filter(trait == traits) %>%
    droplevels(.)
  
  ## run model
  lm <- lm(resp_norm ~ raw_norm * contam, data = df.f)
  ### res vs fit plot
  png(paste0("./indirect/model_outputs/resfits_plots/resfits4_", 
             combs, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lmm)
  layout(1)
  dev.off()
  
  ## model coefficients
  lm_sum <- summary(lm)
  lm_coeff <- as.data.frame(lm_sum$coefficients)
  lm_coeff$trait <- paste0(traits)
  
  ## save
  return(list(lm_coeff))
}