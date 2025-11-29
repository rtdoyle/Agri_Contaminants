mod_species.vs.amend <- function(traits, df){
  
  ### print traits
  print(traits)
  
  ## specify non pea-specific traits
  if (!traits %in% c("wet_pod_weight.corr",
                     "dry_pod_weight.corr",
                     "per_pod_weight",
                     "infected_peas_perc.corr",
                     "NDVI_mean",
                     "pods",
                     "flowers")){
  
  ### subset data to non-spiked samples
  df.f <- df %>%
    filter(spiking_level == 0) %>%
    droplevels(.)
  
  ### model
  lmm <- lm(log(get(traits)) ~ species*amend,
            data = df.f)
  
  }
  
  ## specify pea-specific traits
  else if (traits %in% c("wet_pod_weight.corr",
                         "dry_pod_weight.corr",
                         "per_pod_weight",
                         "infected_peas_perc.corr",
                         "NDVI_mean",
                          "pods",
                          "flowers")){
    
    ### subset data to non-spiked samples
    df.f <- df %>%
      filter(species == "pea" &
               spiking_level == 0) %>%
      droplevels(.)
      ### specify non-count data
      if (!traits %in% c("pods","flowers")){
        
        ### model
        lmm <- lm(log(get(traits)) ~ amend,
                  data = df.f)
        
      }
      ### specify count data
      else if (traits %in% c("pods","flowers")){
        
        ### model
        form <- as.formula(paste(traits, "~ amend"))
        lmm <- glm(form,
                   family = poisson(link = "log"),
                   data = df.f)
        
      }
  }
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m1-species.vs.amend/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  
  
  ### res vs fit plot
  png(paste0(output_dir2, traits, ".png"), 
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
  
  ## specify non pea-specific traits
  if (!traits %in% c("wet_pod_weight.corr",
                     "dry_pod_weight.corr",
                     "per_pod_weight",
                     "NDVI_mean",
                     "infected_peas_perc.corr",
                     "pods",
                     "flowers")){
    
  lmm.emm <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                       amend | species,
                     infer = TRUE,
                     type = "response")
  
  emm.amend <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                        amend,
                      infer = TRUE,
                      type = "response")
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$trait <- paste0(traits)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$trait <- paste0(traits)
  cont.amend <- as.data.frame(emm.amend$contrasts)
  cont.amend$trait <- paste0(traits)
  
  ### return dfs
  return(list(aov, ## [[1]]
              emm, ## [[2]]
              cont, ## [[3]]
              cont.amend ## [[4]]
  )) 
  
  }
  ## specify non pea-specific traits
  else if (traits %in% c("wet_pod_weight.corr",
                     "dry_pod_weight.corr",
                     "per_pod_weight",
                     "NDVI_mean",
                     "infected_peas_perc.corr",
                     "pods",
                     "flowers")){
    
    emm.amend <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                         amend,
                       infer = TRUE,
                       type = "response")
    
    ### processing
    emm <- as.data.frame(emm.amend$emmeans)
    emm$trait <- paste0(traits)
    cont.amend <- as.data.frame(emm.amend$contrasts)
    cont.amend$trait <- paste0(traits)
    
    ### return dfs
    return(list(aov, ## [[1]]
                emm, ## [[2]]
                cont.amend ## [[3]]
    )) 
  }
}