mod_contam <- function(combs, traits, amends, df){
  
  ### print traits
  print(combs)
  
  ### subset data to amendment types
  df.f <- df %>%
    filter(trait == traits &
      amend == amends) %>%
    droplevels(.)
  
  ### specify non-count data
  if (!traits %in% c("leaf_num2")){
    
    ### model
    lmm <- lm(log(emmean) ~ spikeFac, data = df.f)
    
  }
  
  ### specify count data
  else if (traits %in% c("leaf_num2")){
    
    ### model
    lmm <- glm(emmean ~ spikeFac,
               family = poisson(link = "log"),
               data = df.f)
    
  }
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m4-contam/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  

  ### res vs fit plot
  png(paste0(output_dir2, combs, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lmm)
  layout(1)
  dev.off()
  ### model output
  print(summary(lmm))
  aov <- Anova(lmm, type = 2)
  ### processing
  aov$trait <- paste0(traits)
  aov$term <- rownames(aov)
  aov$amend <- paste0(amends)
  
  ### compare treatments within species
  lmm.bt <- update(ref_grid(lmm), tran = "log")
    
  lmm.emm <- emmeans(lmm.bt, poly ~ 
                         spikeFac,
                       infer = TRUE)
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$trait <- paste0(traits)
  emm$amend <- paste0(amends)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$trait <- paste0(traits)
  cont$amend <- paste0(amends)
  
  ### return dfs
  return(list(aov, ## [[1]]
              emm, ## [[2]]
              cont ## [[3]]
  )) 
  
}