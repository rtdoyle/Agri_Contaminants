div_amend.func <- function(combs, times, metrics, df){
  
  ### print traits
  print(combs)
  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  df.f <- df %>%
    filter(timepoint == times) %>%
    droplevels(.)
  
  ### model to compare treatments within each timepoint
  lm <- lm(log(metric_value) ~ species*treat.rn, 
           data = df.f)
  
  ### res vs fit plot
  png(paste0("./16S_outputs/models/resfits_plots/div_amend", 
             combs, ".png"), 
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
  aov$time <- paste0(times)
  
  ### compare spiking levels within source
  lmm.bt <- update(ref_grid(lm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                       treat.rn | species,
                     infer = TRUE,
                     type = "response")
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$metric <- paste0(metrics)
  emm$time <- paste0(times)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$metric <- paste0(metrics)
  cont$time <- paste0(times)
  
  ### return dfs
  return(list(
              aov, ## [[1]]
              emm, ## [[2]]
              cont ## [[3]]
  )) 
}