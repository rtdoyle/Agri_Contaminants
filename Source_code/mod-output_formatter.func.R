format_anova_and_contrasts_flex <- function(results_list, output_dir, 
                                            analysis_type = c("amendment", 
                                                              "contamination")) {
  library(dplyr)
  library(tidyr)
  library(readr)
  library(kableExtra)
  
  analysis_type <- match.arg(analysis_type)
  
  # 1. Combine ANOVA tables across metrics
  aov_combined <- bind_rows(
    lapply(names(results_list), function(m) {
      df <- results_list[[m]][[1]]  # ANOVA table
      df$metric <- m
      df
    })
  )
  
  # 2. Format ANOVA for manuscript (Chi-sq)
  aov_combined <- aov_combined %>%
    mutate(
      pval = signif(`Pr(>Chisq)`, 3),
      sig = ifelse(pval < 0.001, ", p < 0.001",
                   ifelse(pval < 0.01, paste0(", p = ", pval, "**"),
                          ifelse(pval < 0.05, paste0(", p = ", pval, "*"),
                                 ifelse(pval < 0.1, paste0(", p = ", pval, "."),
                                        ", p > 0.1")))),
      stat = signif(Chisq, 3),
      formatted = paste0("χ² [", Df, "] = ", stat, sig)
    )
  
  # 3. Pivot wider so metrics are columns, terms are rows
  aov_wide <- aov_combined %>%
    select(term, metric, formatted) %>%
    pivot_wider(names_from = metric, values_from = formatted)
  
  # 4. Arrange rows in logical order
  if (analysis_type == "amendment") {
    aov_wide <- aov_wide %>%
      arrange(match(term, c("species", "Amendment", "timepoint",
                            "species:Amendment", "species:timepoint",
                            "Amendment:timepoint", "species:Amendment:timepoint")))
  } else {
    aov_wide <- aov_wide %>%
      arrange(match(term, c("species", "spikeFac", "timepoint",
                            "species:spikeFac", "species:timepoint",
                            "spikeFac:timepoint", "species:spikeFac:timepoint")))
  }
  
  # Combine contrasts across metrics (indices 4–6)
  cont_combined <- bind_rows(
    lapply(names(results_list), function(m) {
      df <- bind_rows(results_list[[m]][3:5])  # combine contrasts tables
      df <- df %>% mutate(metric_full = m) ## store original name
      df
    })
  )
  
  
  # Split metric_full into metric and amendment
  cont_combined <- cont_combined %>%
    separate(metric_full, into = c("metric", "amendment"), sep = "_")
  
  # Format contrasts
  cont_combined <- cont_combined %>%
    mutate(
      pval = signif(p.value, 3),
      sig = ifelse(p.value < 0.001, ", p < 0.001",
                   ifelse(p.value < 0.01, paste0(", p = ", pval, "**"),
                          ifelse(p.value < 0.05, paste0(", p = ", pval, "*"),
                                 ifelse(p.value < 0.1, paste0(", p = ", pval, "."),
                                        ", p > 0.1")))),
      stat = signif(t.ratio, 3),
      df.r = signif(df, 1),
      stat.f = paste0("t = ", stat, ", df = ", df.r, sig),
      eff_fmt = paste0(signif(eff, 3), " [",
                         signif(LCL, 3), " to ",
                         signif(UCL, 3), "]")
    )
  
  # contrasts: type = "response" shows ratio
  cont_summary <- cont_combined %>%
    filter(p.value < 0.1) %>%
    mutate(
      contrast_id = paste0(contrast, ": ", species, "_", treat, "_", timepoint),
      contrast_info = paste0(stat.f, "; Est. [95% CL] = ",
                             eff_fmt)
    ) %>%
    select(metric, amendment, contrast_id, contrast_info) %>%
    pivot_wider(
      names_from = c("metric"),
      values_from = contrast_info)
  
  # Contrast summary table
  contrast_table <- cont_summary %>%
    kable(format = "html", escape = FALSE,
          caption = "Significant contrasts (p < 0.1)") %>%
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
    column_spec(2, width = "12cm")
  print(contrast_table)
  
  # Combine emmeans across metrics
  emm_combined <- bind_rows(
    lapply(names(results_list), function(m) {
      df <- results_list[[m]][[2]]  # emmeans
      df$metric <- m
      df
    })
  )
  
  # Save outputs
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write_excel_csv(aov_wide, file.path(output_dir, paste0("anova.alpha_", analysis_type, ".csv")))
  write_excel_csv(cont_combined, file.path(output_dir, paste0("contrasts.alpha_", analysis_type, ".csv")))
  write_excel_csv(cont_summary, file.path(output_dir, paste0("contrasts.sum.alpha_", analysis_type, ".csv")))
  write_excel_csv(emm_combined, file.path(output_dir, paste0("emmeans.alpha_", analysis_type, ".csv")))
  
  message("Saved ANOVA + contrasts + emmeans tables for ", analysis_type, " analysis.")
  
  return(list(aov_table = aov_wide, contrasts = contrast_table, emmeans = emm_combined))
}
  