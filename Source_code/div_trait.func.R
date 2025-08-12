div_trait.func <- function(combs, traits, metrics, df){
  
  ## print trait-metric comb
  print(combs)
  
  ## plot
  p <-
    ggplot(data = df,
           aes(y = get(metrics), x = get(traits))) +
    geom_point(size = 3, alpha = 0.8,
               aes(colour = treat.rn)) + 
    geom_vline(aes(xintercept = 0), linetype = 2) +
    geom_smooth(method = "lm", colour = "black",
                linetype = 2, se = FALSE, formula = 'y ~ x') +
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    guides(color = "none", shape = "none") +
    ggtitle(NULL) + 
    labs(x = paste0(traits), y = NULL) +
    # facet_wrap(~Trait, scales = "free", ncol = 2) +  
    theme_bw() + 
    theme(axis.title.x = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  ### linear model
  lm <- lm(log(get(metrics)) ~ get(traits), data = df)
  
  ### res vs fit plot
  png(paste0("./16S_outputs/models/div_trait.resfits_", 
             combs, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lm)
  layout(1)
  dev.off()
  ### model coefficients
  lm_sum <- summary(lm)
  lm_coeff <- as.data.frame(lm_sum$coefficients)
  lm_coeff$trait <- paste0(traits)
  
  ### outputs
  return(list(p, lm_coeff))
  
} 