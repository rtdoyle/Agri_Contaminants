div_vs_trait.func <- function(combs, metrics, traits, df){
  
  ## print metric
  print(combs)
  
  # Convert metric to symbol for tidy evaluation
  metric <- df[[metrics]]
  trait <- df[[traits]]

### linear model
lm <- lm(log(metric) ~ trait, data = df)

### res vs fit plot
png(paste0("./16S_outputs/models/resfits_plots/div_vs_trait.resfits_", 
           combs, ".png"), 
    width=6, height=6, units='in', res=300)
layout(matrix(1:4, ncol = 2))
plot(lm)
layout(1)
dev.off()
### model coefficients
lm_sum <- summary(lm)
lm_coeff <- as.data.frame(lm_sum$coefficients)
lm_coeff$metric <- paste0(metrics)
lm_coeff$trait <- paste0(traits)

return(list(lm_coeff))
}

