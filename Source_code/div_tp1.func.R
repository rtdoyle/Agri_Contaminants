div_tp1.func <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  richness_sum <- df %>%
    group_by(species, treat.rn, tp_treat, species_tp, species_treat) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE),
              min = min(get(metrics), na.rm = TRUE),
              max = max(get(metrics), na.rm = TRUE))
  
  p <- richness_sum %>%
    ggplot(aes(x = species, y = mean, colour = treat.rn, shape = species)) +
    geom_pointrange(
      size = 1, aes(ymin = min,
                    ymax = max),
      position = position_dodge2(0.5)) + 
    geom_hline(data = richness_sum %>%
                 filter(species_treat == "no_plant-Garden Soil"),
               aes(yintercept = mean),
               linetype = 2) +
    scale_color_manual(values = c(
      "gray",
      "black",
      "brown",
      "green3",
      "green4",
      "lightblue",
      "cornflowerblue",
      "blue",
      "#EABD8C",
      "#FFAD00",
      "#B06500")) +
    guides(color = "none", shape = "none") +
    scale_shape_manual(values = c(1,16,17,15),
                       labels = c("No Plant", "Lettuce", 
                                  "Pea", "Radish")) +
    scale_x_discrete(labels = c(
      "no_plant" = "No Plant",
      "L" = "Lettuce",
      "P" = "Pea",
      "R" = "Radish")) +
    ggtitle(NULL) + 
    labs(x = NULL, y = NULL) +
    theme_bw() + 
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  ### model to compare treatments
  lm <- lm(log(metric_value) ~ species_treat, data = df)
  
  ### res vs fit plot
  png(paste0("./16S_outputs/models/div_tp1.resfits_", 
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
                       species_treat,
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