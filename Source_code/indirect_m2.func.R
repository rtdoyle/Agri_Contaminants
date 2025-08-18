indirect_m2.func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### filter to trait
  df.f <- df %>%
    filter(trait == traits) %>%
    droplevels(.)
  
  ### model to compare treatments
  lmm <- lm(log(response) ~ contam, data = df.f)
  
  ### res vs fit plot
  png(paste0("./indirect/model_outputs/resfits_plots/resfits2_", 
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
  
  ### compare spiking levels within source
  lmm.bt <- update(ref_grid(lmm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ contam,
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