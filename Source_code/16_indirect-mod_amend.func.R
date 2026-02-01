mod_amend <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### filter to trait
  df.f <- df %>%
    filter(trait == traits) %>%
    droplevels(.)
  
  ### specify non-count data
  if (!traits %in% c("leaf_num2")){
    
    ### model
    lmm <- lm(log(emmean) ~ amend, data = df.f)
    
  }
  
  ### specify count data
  else if (traits %in% c("leaf_num2")){
    
    ### model
    lmm <- glm(emmean ~ amend,
               family = poisson(link = "log"),
               data = df.f)
    
  }
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m3-amend/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  
  
  ### res vs fit plot
  png(paste0(output_dir2, traits, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lmm)
  layout(1)
  dev.off()
  
  ## print off summary  
  print(summary(lmm))
  aov <- Anova(lmm, type = 2)
  
  ### processing
  aov$trait <- paste0(traits)
  aov$term <- rownames(aov)
  
  ### compare spiking levels within source
  lmm.bt <- update(ref_grid(lmm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ amend,
                     infer = TRUE,
                     type = "response")
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$trait <- paste0(traits)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$trait <- paste0(traits)
  
  ### return dfs
  return(list(aov, ## [[1]]
              emm, ## [[2]]
              cont ## [[3]]
  )) 
  
}