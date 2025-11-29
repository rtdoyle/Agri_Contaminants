div.RR_species.vs.contam <- function(combs, metrics, amends, df){
  
  ### print traits
  print(combs)
  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  df.f <- df %>%
    filter(Amendment == amends) %>%
    droplevels(.)
  
  ### model to compare contamination lvl within each timepoint + amendment type
  lm <- lm(metric_value ~ species*spikeFac, 
           data = df.f)
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m2b-species.vs.contam_RR/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  
  
  ### res vs fit plot
  png(paste0(output_dir2, combs, ".png"), 
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
  aov$amend <- paste0(amends)
  
  ### compare spiking levels within source
  #lmm.bt <- update(ref_grid(lm), tran = "log")
  lmm.emm <- emmeans(lm, poly ~ 
                       spikeFac | species,
                     infer = TRUE)
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$metric <- paste0(metrics)
  emm$amend <- paste0(amends)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$metric <- paste0(metrics)
  cont$amend <- paste0(amends)
  
  ### return dfs
  return(list(
    aov, ## [[1]]
    emm, ## [[2]]
    cont ## [[3]]
  )) 
}