div_amend.func <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  ### model to compare treatments
  lm <- lm(log(metric_value) ~ species*treat.rn*timepoint, 
           data = df)
  
  ### res vs fit plot
  png(paste0("./16S_outputs/models/div_treat.resfits_", 
             metrics, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lm)
  layout(1)
  dev.off()
  
  ## print off summary  
  print(summary(lm))
  aov <- Anova(lm, type = 3)
  
  ### processing
  aov$metric <- paste0(metrics)
  aov$term <- rownames(aov)
  
  ### compare spiking levels within source
  lmm.bt <- update(ref_grid(lm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                       timepoint | treat.rn + species,
                     infer = TRUE,
                     type = "response")
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$metric <- paste0(metrics)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$metric <- paste0(metrics)
  
  ### return dfs
  return(list(
              aov, ## [[1]]
              emm, ## [[2]]
              cont ## [[3]]
  )) 
}