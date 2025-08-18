indirect_m3.func <- function(combs, traits, amends, df){
  
  ### print traits
  print(combs)
  
  ### subset data to amendment types
  df.f <- df %>%
    filter(trait == traits &
      contam == amends) %>%
    droplevels(.)
  
  ### model
  lmm <- lm(log(response) ~ spikeFac,
            data = df.f)
    

  ### res vs fit plot
  png(paste0("./indirect/model_outputs/resfits_plots/resfits3_", 
             combs, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lmm)
  layout(1)
  dev.off()
  ### model output
  print(summary(lmm))
  aov <- Anova(lmm, type = 2)
  ### processing
  aov$trait <- paste0(traits)
  aov$term <- rownames(aov)
  aov$amend <- paste0(amends)
  
  ### compare treatments within species
  lmm.bt <- update(ref_grid(lmm), tran = "log")
    
  lmm.emm <- emmeans(lmm.bt, poly ~ 
                         spikeFac,
                       infer = TRUE)
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$trait <- paste0(traits)
  emm$amend <- paste0(amends)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$trait <- paste0(traits)
  cont$amend <- paste0(amends)
  
  ### return dfs
  return(list(aov, ## [[1]]
              emm, ## [[2]]
              cont ## [[3]]
  )) 
  
}