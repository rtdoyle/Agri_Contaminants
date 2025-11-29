mod_species.vs.contam <- function(combs, traits, amends, df){
  
  ### print traits
  print(combs)
  
  ### non-pea specific traits
  if (!traits %in% c("wet_pod_weight.corr",
                     "dry_pod_weight.corr",
                     "NDVI_mean",
                     "infected_peas_perc.corr",
                     "pods",
                     "per_pod_weight",
                     "flowers")){
    
    ### subset data to amendment types
    df.f <- df %>%
      filter(amend == amends) %>%
      droplevels(.)
    
    ### model
    lmm <- lm(log(get(traits)) ~ species*spikeFac,
              data = df.f)
    
  }
  
  ### pea-specific traits
  else if (traits %in% c("wet_pod_weight.corr",
                         "dry_pod_weight.corr",
                         "NDVI_mean",
                         "infected_peas_perc.corr",
                         "pods",
                         "per_pod_weight",
                         "flowers")){
    
    ### subset data to peas
    df.f <- df %>%
      filter(species == "pea" &
               amend == amends) %>%
      droplevels(.)
    
    if (!traits %in% c("pods","flowers")){
      
      ### model
      lmm <- lm(log(get(traits)) ~ spikeFac,
                data = df.f)
      
    }
    
    else if (traits %in% c("pods","flowers")){
      
      ### model
      form <- as.formula(paste(traits, "~ spikeFac"))
      lmm <- glm(form,
                 family = poisson(link = "log"),
                 data = df.f)
      
    }
  }
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/m2-species.vs.contam/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  
  
  ### res vs fit plot
  png(paste0(output_dir2, combs, ".png"), 
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
  aov$amend <- paste0(amends)
  
  ### compare treatments within species
  
  ### specify for pea vs other species
  if (!traits %in% c("wet_pod_weight.corr",
                     "dry_pod_weight.corr",
                     "NDVI_mean",
                     "infected_peas_perc.corr",
                     "pods",
                     "per_pod_weight",
                     "flowers")){
    
    lmm.emm <- emmeans(lmm, poly ~ 
                         spikeFac | species,
                       infer = TRUE)
    
    emm.SL <- emmeans(lmm, poly ~ 
                        spikeFac,
                         infer = TRUE)
    
    ### processing
    emm <- as.data.frame(lmm.emm$emmeans)
    emm$trait <- paste0(traits)
    emm$amend <- paste0(amends)
    cont <- as.data.frame(lmm.emm$contrasts)
    cont$trait <- paste0(traits)
    cont$amend <- paste0(amends)
    ## across species 
    cont.SL <- as.data.frame(emm.SL$contrasts)
    cont.SL$trait <- paste0(traits)
    cont.SL$amend <- paste0(amends)
    
    ### return dfs
    return(list(aov, ## [[1]]
                emm, ## [[2]]
                cont, ## [[3]]
                cont.SL ## [[4]]
    )) 
  }
  
  else if (traits %in% c("wet_pod_weight.corr",
                         "dry_pod_weight.corr",
                         "NDVI_mean",
                         "infected_peas_perc.corr",
                         "pods",
                         "per_pod_weight",
                         "flowers")){
    
    emm.SL <- emmeans(lmm, poly ~ 
                        spikeFac,
                      infer = TRUE)
    
    ### processing
    emm <- as.data.frame(emm.SL$emmeans)
    emm$trait <- paste0(traits)
    emm$amend <- paste0(amends)
    ## across species 
    cont.SL <- as.data.frame(emm.SL$contrasts)
    cont.SL$trait <- paste0(traits)
    cont.SL$amend <- paste0(amends)
    
    ### return dfs
    return(list(aov, ## [[1]]
                emm, ## [[2]]
                cont.SL ## [[3]]
    )) 
  }
  
  
  
}