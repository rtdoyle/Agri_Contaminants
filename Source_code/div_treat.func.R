div_treat.func <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  ### plot
  richness_sum <- df %>%
    group_by(timepoint, species, treat.rn, tp_treat, species_tp) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE),
              min = min(get(metrics), na.rm = TRUE),
              max = max(get(metrics), na.rm = TRUE))
  
  # Custom labels for the 'species_tp' facet
  labs <- c(
    no_plant = "No Plant",
    L_TP1 = "Lettuce @ Planting",
    L_TP2 = "Lettuce @ Harvest",
    P_TP1 = "Pea @ Planting",
    P_TP2 = "Pea @ Harvest",
    R_TP1 = "Radish @ Planting",
    R_TP2 = "Radish @ Harvest"
  )

  ## change order of species_tp
  richness_sum$species_tp <- factor(richness_sum$species_tp,
                                    levels = c(
                                      "L_TP1", "P_TP1", "R_TP1",
                                      "L_TP2", "P_TP2", "R_TP2"
                                     ))
  
  p <- richness_sum %>%
    ggplot(aes(x = treat.rn, y = mean, colour = treat.rn, shape = species)) +
    geom_pointrange(aes(ymin = min,
                        ymax = max)) + 
    geom_hline(data = richness_sum %>%
                 filter(treat.rn == "control-0"),
               aes(yintercept = mean),
               linetype = 2) +
    scale_color_manual(values = c(
      "gray",
      "lightblue",
      "cornflowerblue",
      "blue",
      "#EABD8C",
      "#FFAD00",
      "#B06500")) +
    guides(color = "none", shape = "none") +
    scale_shape_manual(values = c(16,17,15),
                       labels = c("Lettuce", "Pea", "Radish")) + 
    ggtitle(NULL) + 
    labs(x = NULL, y = NULL) +
    facet_wrap(~species_tp, ncol = 6,
               scales = "fixed",
               labeller = labeller("species_tp" = labs)) +  
    #coord_flip() +
    theme_bw() + 
    theme(
      strip.text = element_text(size = 16, face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  ### model to compare treatments
  lm <- lm(log(metric_value) ~ species*timepoint*treat.rn, data = df)
  
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
                       treat.rn | species + timepoint,
                     infer = TRUE,
                     type = "response")
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$metric <- paste0(metrics)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$metric <- paste0(metrics)
  
  ### return dfs
  return(list(p, 
              aov, ## [[2]]
              emm, ## [[3]]
              cont ## [[4]]
  )) 
}