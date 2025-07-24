lmm1_func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### subset data to non-spiked samples
  df.f <- df %>%
    filter(spiking_level == 0) %>%
    droplevels(.)
  
  ### model
  lmm <- lm(log(get(traits)) ~ species*contam,
            data = df.f)
  ### res vs fit plot
  png(paste0("./model_outputs/direct/resfits1_", 
             traits, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lmm)
  layout(1)
  dev.off()
  ### model output
  print(summary(lmm))
  aov <- Anova(lmm, type = 3)
  ### processing
  aov$trait <- paste0(traits)
  aov$term <- rownames(aov)
  
  ### compare treatments within species
  lmm.bt <- update(ref_grid(lmm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                       contam | species,
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