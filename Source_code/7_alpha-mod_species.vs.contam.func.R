library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(dplyr)

div_species.vs.contam <- function(combs, metrics, amends, df, output_dir) {
  
  # Print combination info
  print(combs)
  
  # Dynamically assign metric column
  df$metric_value <- df[[metrics]]
  
  # Subset data for given amendment and exclude controls
  df.f <- df %>%
    filter(Amendment == amends & treat.rn != "control-0") %>%
    droplevels()
  
  # Fit mixed model: contamination level (spikeFac) within species and timepoint
  lmm <- lmer(log(metric_value) ~ species * spikeFac * timepoint + (1 | pot_ID),
              data = df.f)
  
  # Create output directory if missing
  output_dir2 <- file.path(output_dir, "resfits_plots/mod-contam.sp.time/")
  if (!dir.exists(output_dir2)) dir.create(output_dir2, recursive = TRUE)
  
  # Save diagnostic plots (4-panel)
  png(file.path(output_dir2, paste0(combs, "_diagnostics.png")),
      width = 8, height = 8, units = "in", res = 300)
  par(mfrow = c(2, 2))
  plot(fitted(lmm), resid(lmm), main = "Residuals vs Fitted",
       xlab = "Fitted values", ylab = "Residuals")
  abline(h = 0, col = "red")
  qqnorm(resid(lmm), main = "Normal Q-Q")
  qqline(resid(lmm), col = "red")
  std_resid <- resid(lmm) / sd(resid(lmm))
  plot(fitted(lmm), sqrt(abs(std_resid)), main = "Scale-Location",
       xlab = "Fitted values", ylab = "√|Standardized residuals|")
  hist(resid(lmm), main = "Histogram of Residuals", xlab = "Residuals")
  dev.off()
  
  # generate summary
  mod.sum <- summary(lmm)
  
  # Type III ANOVA for fixed effects
  aov <- Anova(lmm, type = 3)
  aov$metric <- metrics
  aov$term <- rownames(aov)
  aov$amend <- amends
  print(aov)
  
  # Estimated marginal means 
  lmm.bt <- update(ref_grid(lmm), tran = "log")
  
  ### for contamination level (don't use response scale)
  lmm.contam <- emmeans(lmm, poly ~ 
                          spikeFac | species + timepoint, 
                        infer = TRUE,
                        type = "response")
  emm.contam <- as.data.frame(lmm.contam$emmeans)
  emm.contam$metric <- metrics
  emm.contam$amend <- amends
  
  # Contrasts for contamination level
  cont.contam <- as.data.frame(lmm.contam$contrasts)
  cont.contam$metric <- metrics
  cont.contam$amend <- amends
  cont.contam$eff <- log2(cont.contam$ratio)
  cont.contam$LCL <- log2(cont.contam$lower.CL)
  cont.contam$UCL <- log2(cont.contam$upper.CL)
  cont.contam <- cont.contam %>%
    select(-ratio)
  
  # Estimated marginal means for species level
  lmm.sp <- emmeans(lmm.bt, pairwise ~ 
                      species | spikeFac + timepoint, 
                    infer = TRUE,
                    type = "response")
  emm.sp <- as.data.frame(lmm.sp$emmeans)
  emm.sp$metric <- metrics
  emm.sp$amend <- amends
  
  # Contrasts for species level
  cont.sp <- as.data.frame(lmm.sp$contrasts)
  cont.sp$metric <- metrics
  cont.sp$amend <- amends
  cont.sp$treat <- cont.sp$spikeFac
  cont.sp$eff <- log2(cont.sp$ratio)
  cont.sp$LCL <- log2(cont.sp$lower.CL)
  cont.sp$UCL <- log2(cont.sp$upper.CL)
  cont.sp <- cont.sp %>%
    select(-spikeFac, -ratio)
  
  # Estimated marginal means for time level
  lmm.time <- emmeans(lmm.bt, trt.vs.ctrl ~ 
                        timepoint | spikeFac + species, 
                      infer = TRUE,
                      type = "response")
  emm.time <- as.data.frame(lmm.time$emmeans)
  emm.time$metric <- metrics
  emm.time$amend <- amends
  
  # Contrasts for time level
  cont.time <- as.data.frame(lmm.time$contrasts)
  cont.time$metric <- metrics
  cont.time$amend <- amends
  cont.time$treat <- cont.time$spikeFac
  cont.time$eff <- log2(cont.time$ratio)
  cont.time$LCL <- log2(cont.time$lower.CL)
  cont.time$UCL <- log2(cont.time$upper.CL)
  cont.time <- cont.time %>%
    select(-spikeFac, -ratio)
  
  # Return results
  return(list(
              aov,       # [[1]] ANOVA table
              emm.contam,    # [[2]] Estimated marginal means
              cont.contam,   # [[3]]  Contrasts (poly)
              cont.sp,   # [[4]] pw  Contrasts (pw)
              cont.time # [[5]] time  (trt.crtl)
  ))
}