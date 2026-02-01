library(lmerTest)
library(car)
library(dplyr)
library(tibble)
library(stringr)
library(MuMIn)  # For pseudo-R²
library(emmeans)

div_species.vs.time <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  df.f <- df %>%
    ## filter(timepoint == times) %>% ## no filtering
    droplevels(.)
  
  ### model to compare treatments within each timepoint
  lm <- lmer(log(metric_value) ~ species*timepoint + (1 | pot_ID), 
           data = df.f)
  
  # Create output directory if it doesn't exist
  ## set output directory
  output_dir2 <- paste0(output_dir, "resfits_plots/mod-species.vs.time/")
  if (!dir.exists(output_dir2)) {
    dir.create(output_dir2, recursive = TRUE)
  }  
  
  ### res vs fit plot
  png(paste0(output_dir2, metrics, ".png"), 
      width=6, height=6, units='in', res=300)
  layout(matrix(1:4, ncol = 2))
  plot(lm)
  layout(1)
  dev.off()
  
  ## print off summary  
  aov <- Anova(lm, type = 3)
  
  ### processing
  aov$metric <- paste0(metrics)
  aov$term <- rownames(aov)
  #aov$time <- paste0(times)
  
  ### compare spiking levels within source
  lmm.bt <- update(ref_grid(lm), tran = "log")
  lmm.emm <- emmeans(lmm.bt, pairwise ~ 
                       species | timepoint,
                     infer = TRUE,
                     type = "response")
  
  ### processing
  emm <- as.data.frame(lmm.emm$emmeans)
  emm$metric <- paste0(metrics)
  #emm$time <- paste0(times)
  cont <- as.data.frame(lmm.emm$contrasts)
  cont$metric <- paste0(metrics)
  #cont$time <- paste0(times)
  
  ### return dfs
  return(list(
    aov, ## [[1]]
    emm, ## [[2]]
    cont ## [[3]]
  )) 
}