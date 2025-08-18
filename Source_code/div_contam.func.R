div_contam.func <- function(combs, times, metrics, amends, df){
  
  ### print traits
  print(combs)
  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  df.f <- df %>%
    filter(timepoint == times &
             Source == amends) %>%
    droplevels(.)
  
  ### model to compare contamination lvl within each timepoint + amendment type
  lm <- lm(log(metric_value) ~ species*spikeFac, 
           data = df.f)
  
  ### res vs fit plot
  png(paste0("./16S_outputs/models/resfits_plots/div_contam-", 
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
  aov$amend <- paste0(amends)
  
  ### compare spiking levels within source
  #lmm.bt <- update(ref_grid(lm), tran = "log")
  lmm.emm <- emmeans(lm, poly ~ 
                       spikeFac | species,
                     infer = TRUE)
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$metric <- paste0(metrics)
  emm$time <- paste0(times)
  emm$amend <- paste0(amends)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$metric <- paste0(metrics)
  cont$time <- paste0(times)
  cont$amend <- paste0(amends)
  
  ### return dfs
  return(list(
              aov, ## [[1]]
              emm, ## [[2]]
              cont ## [[3]]
  )) 
}