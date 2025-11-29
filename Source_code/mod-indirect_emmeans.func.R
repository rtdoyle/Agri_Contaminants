mod_indirect_emms <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ## split out count data from non-count data
  ### specify non-count data
  if (!traits %in% c("leaf_num","leaf_num2")){
    
    ### model
    lmm <- lm(log(get(traits)) ~ treatment, data = df)
    
  }
  ### specify count data
  else if (traits %in% c("leaf_num","leaf_num2")){
    
    ### model
    form <- as.formula(paste(traits, "~ treatment"))
    lmm <- glm(form,
               family = poisson(link = "log"),
               data = df)
    
  }
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m1-emmeans/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }
  
  ### res vs fit plots
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
  
  ## emmeans
  lmm.bt <- update(ref_grid(lmm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ treatment,
                     infer = TRUE,
                     type = "response")
  
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$trait <- paste0(traits)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$trait <- paste0(traits)
  
  ## add in contam
  emm$amend <- df$amend[match(emm$treatment, df$treatment)]
  
  ## add in spike
  emm$spikeFac <- df$spikeFac[match(emm$treatment, df$treatment)]
  
  ## add in contam_spike
  emm$amend_spike <- df$amend_spike[match(emm$treatment, df$treatment)]
  
  ### return dfs
  return(list(aov, ## [[1]]
              emm,   ## [[2]]
              cont)) ## [[3]] 
  
}