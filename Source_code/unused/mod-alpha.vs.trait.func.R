library(dplyr)
library(tibble)
library(purrr)

univariate_diversity <- function(data, metrics, traits, 
                                 species_var = "species", 
                                 time_var = "timepoint") {
  
  all_results <- list()
  
  for (metric in metrics) {
    metric_results <- list()
    
    for (sp in unique(data[[species_var]])) {
      for (tp in unique(data[[time_var]])) {
        
        # Subset data for species × timepoint
        subset_df <- data %>%
          filter(.data[[species_var]] == sp, .data[[time_var]] == tp)
        
        if (nrow(subset_df) == 0) next
        
        # Fit univariate models for each trait
        trait_results <- map_df(traits, function(trait) {
          mod <- lm(as.formula(paste("log(", metric, ") ~", trait)), 
                    data = subset_df)
          
          # Extract ANOVA info
          aov_tab <- anova(mod) %>%
            as.data.frame() %>%
            slice(1) %>%  # Only trait row
            mutate(term = trait)
          
          # Extract coefficients
          coef_tab <- summary(mod)$coefficients
          intercept <- coef_tab[1, "Estimate"]
          intercept_se <- coef_tab[1, "Std. Error"]
          slope <- coef_tab[2, "Estimate"]
          slope_se <- coef_tab[2, "Std. Error"]
          
          # Format slope/intercept summary
          slope_summary <- sprintf("%.2g (± %.2g)", slope, slope_se)
          intercept_summary <- sprintf("%.2g (± %.2g)", intercept, intercept_se)
          equation <- sprintf("y = %s + %s", slope_summary, intercept_summary)
          
          aov_tab <- aov_tab %>%
            mutate(intercept = intercept,
                   intercept_se = intercept_se,
                   slope = slope,
                   slope_se = slope_se,
                   equation = equation)
          
          return(aov_tab)
        })
        
        trait_results <- trait_results %>%
          mutate(metric = metric, species = sp, timepoint = tp)
        
        metric_results[[paste(sp, tp, sep = "_")]] <- trait_results
      }
    }
    
    metric_combined <- bind_rows(metric_results)
    all_results[[metric]] <- metric_combined
  }
  
  combined <- bind_rows(all_results)
  
  # Dynamically find p-value column
  p_col <- grep("^Pr", names(combined), value = TRUE)
  
  if (length(p_col) == 0) {
    warning("No p-value column found in ANOVA output.")
    combined$p_value <- NA
  } else {
    combined <- combined %>% rename(p_value = all_of(p_col[1]))
  }
  
  # Apply BH correction and format
  combined <- combined %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      p_formatted = case_when(
        p_value < 0.001 ~ "<0.001",
        p_value < 0.01  ~ "<0.01",
        p_value < 0.05  ~ "<0.05",
        p_value < 0.1   ~ "<0.1",
        p_value > 0.1   ~ "ns",
        TRUE ~ sprintf("%.3g", p_value)
      ),
      sig = ifelse(p_adj > 0.05, "", "*"),
      p_display = paste0(p_formatted, sig), # Combine p-value and asterisk
      stat = sprintf("F = %.2g, p = %s, df = %.2g", `F value`, p_display, Df)
    ) %>%
    select(metric, species, timepoint, term, stat, p_adj, equation)
  
  return(combined)
}