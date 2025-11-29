library(lme4)
library(lmerTest)

div_species.vs.amend <- function(metrics, df, output_dir){
  
  ### print traits
  print(metrics)
  
  ## model
  df$metric_value <- df[[metrics]]  # dynamically assign the column to a new variable
  
  df.f <- df %>%
    filter(spikeFac == "SL-0") %>%
    droplevels(.)
  
  ### model to compare treatments within each timepoint
  lm <- lmer(log(metric_value) ~ species*Amendment*timepoint + (1 | pot_ID), 
           data = df.f)
  
  # Create output directory if missing
  output_dir2 <- file.path(output_dir, "resfits_plots/mod-amend.sp.time/")
  if (!dir.exists(output_dir2)) dir.create(output_dir2, recursive = TRUE)
  
  # Open PNG device
  png(file.path(output_dir2, paste0(metrics, "_diagnostics.png")),
      width = 8, height = 8, units = "in", res = 300)
  
  par(mfrow = c(2, 2))  # 4 plots in one page
  
  # 1. Residuals vs Fitted
  plot(fitted(lm), resid(lm), main = "Residuals vs Fitted",
       xlab = "Fitted values", ylab = "Residuals")
  abline(h = 0, col = "red")
  
  # 2. Normal Q-Q
  qqnorm(resid(lm), main = "Normal Q-Q")
  qqline(resid(lm), col = "red")
  
  # 3. Scale-Location (sqrt standardized residuals vs fitted)
  std_resid <- resid(lm) / sd(resid(lm))
  plot(fitted(lm), sqrt(abs(std_resid)), main = "Scale-Location",
       xlab = "Fitted values", ylab = "√|Standardized residuals|")
  
  # 4. Histogram of residuals (since leverage isn't standard for lmer)
  hist(resid(lm), main = "Histogram of Residuals", xlab = "Residuals")
  
  dev.off()
  
  ## print off summary  
  aov <- Anova(lm, type = 3)
  print(aov)
  
  ### processing
  aov$metric <- paste0(metrics)
  aov$term <- rownames(aov)
  
  ## emmeans
  lmm.bt <- update(ref_grid(lm), tran = "log")
  
  ### compare amendment diffs within species and timepoint
  lmm.emm.amend <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                       Amendment | species + timepoint,
                     infer = TRUE,
                     type = "response")
  lmm.emm.sp <- emmeans(lmm.bt, pairwise ~ 
                             species | Amendment + timepoint,
                           infer = TRUE,
                           type = "response")
  lmm.emm.time <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                            timepoint | Amendment + species,
                          infer = TRUE,
                          type = "response")
  
  ### amend
  emm.amend <- as.data.frame(lmm.emm.amend$emmeans)
  emm.amend$metric <- paste0(metrics)
  cont.amend <- as.data.frame(lmm.emm.amend$contrasts)
  cont.amend$metric <- paste0(metrics)
  cont.amend$treat <- "SL-0"
  cont.amend$eff <- log2(cont.amend$ratio)
  cont.amend$LCL <- log2(cont.amend$lower.CL)
  cont.amend$UCL <- log2(cont.amend$upper.CL)
  cont.amend <- cont.amend %>%
    select(-ratio)
  
  ## species
  emm.sp <- as.data.frame(lmm.emm.sp$emmeans)
  emm.sp$metric <- paste0(metrics)
  cont.sp <- as.data.frame(lmm.emm.sp$contrasts)
  cont.sp$metric <- paste0(metrics)
  cont.sp$treat <- cont.sp$Amendment
  cont.sp$eff <- log2(cont.sp$ratio)
  cont.sp$LCL <- log2(cont.sp$lower.CL)
  cont.sp$UCL <- log2(cont.sp$upper.CL)
  cont.sp <- cont.sp %>%
    select(-ratio)
    
  ## time
  emm.time <- as.data.frame(lmm.emm.time$emmeans)
  emm.time$metric <- paste0(metrics)
  cont.time <- as.data.frame(lmm.emm.time$contrasts)
  cont.time$metric <- paste0(metrics)
  cont.time$treat <- cont.time$Amendment
  cont.time$eff <- log2(cont.time$ratio)
  cont.time$LCL <- log2(cont.time$lower.CL)
  cont.time$UCL <- log2(cont.time$upper.CL)
  cont.time <- cont.time %>%
    select(-ratio)
  
  ### return dfs
  return(list(
              aov, ## [[1]]
              emm.amend, ## [[2]]
              cont.amend, ## [[3]]
              cont.sp, ## [[4]]
              cont.time ## [[5]]
  )) 
}