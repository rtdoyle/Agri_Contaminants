mod_treat.v.time <- function(traits, df){
  
  ### print traits
  print(traits)
  
  df.f <- df %>%
    filter(trait == traits) %>%
    droplevels(.)
  
  ### model to compare treatments within each timepoint
  lm <- lm(log(emmean) ~ amend_spike * weekFac,  
           data = df.f)
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m2-treat.vs.time/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }
  
  ### res vs fit plot
  png(paste0(output_dir2, traits, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lm)
  layout(1)
  dev.off()
  
  ## print off summary  
  print(summary(lm))
  aov <- Anova(lm, type = 3)
  
  ### processing
  aov$trait <- paste0(traits)
  aov$term <- rownames(aov)
  
  ### compare spiking levels within source
  lmm.bt <- update(ref_grid(lm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, pairwise ~ 
                       amend_spike | weekFac,
                     infer = TRUE,
                     type = "response")
  
  lmm.emm.poly <- emmeans(lm, poly ~
                       weekFac | amend_spike,
                     infer = TRUE)
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$trait <- paste0(traits)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$trait <- paste0(traits)
  cont.p <- as.data.frame(lmm.emm.poly$contrasts)
  cont.p$trait <- paste0(traits)
  
  ### return dfs
  return(list(
    aov, ## [[1]]
    emm, ## [[2]]
    cont, ## [[3]]
    cont.p ## [[4]]
  )) 
}