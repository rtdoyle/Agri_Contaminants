emms_func <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ### model to calc geometric means for each slurry
  lmm <- lm(log(get(traits)) ~ treatment, data = df)
  
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
  return(list(emm,   ## [[1]]
              cont)) ## [[2]] 
  
}