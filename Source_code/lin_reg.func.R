lin_reg.func <- function(traits, df){
  
  ## print traits
  print(paste0(traits))
  
  ## filter for each trait
  df.f <- df %>%
    filter(trait == traits) %>%
    droplevels(.)
  
  ## run model
  lm <- lm(resp_norm ~ raw_norm * contam, data = df.f)
  lm_sum <- summary(lm)
  lm_coeff <- as.data.frame(lm_sum$coefficients)
  lm_coeff$trait <- paste0(traits)
  
  ## save
  return(list(lm_coeff))
}