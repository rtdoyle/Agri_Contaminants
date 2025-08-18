emms.func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### model to calc geometric means for each slurry
  lmm <- lm(log(get(traits)) ~ treatment, data = df)
  
  ### res vs fit plots
  png(paste0("./indirect/model_outputs/resfits_plots/resfits1_", 
             traits, ".png"), 
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
  emm$contam <- df$contam[match(emm$treatment, df$treatment)]
  
  ## add in spike
  emm$spikeFac <- df$spikeFac[match(emm$treatment, df$treatment)]
  
  ## add in contam_spike
  emm$contam_spike <- df$contam_spike[match(emm$treatment, df$treatment)]
  
  ### return dfs
  return(list(aov, ## [[1]]
              emm,   ## [[2]]
              cont)) ## [[3]] 
  
}