library(lmerTest)   # Adds p-values to lmer models
library(car)
library(dplyr)
library(tibble)
library(stringr)
library(MuMIn)      # For pseudo-R²

forward_select_lmer <- function(response, fixed_terms, random_term, data, metric_name) {
  
  # Start with random intercept only
  null_formula <- as.formula(paste(response, "~ 1 + (1 |", random_term, ")"))
  current_model <- lmer(null_formula, data = data, REML = FALSE)
  
  remaining_terms <- fixed_terms
  selected_terms <- c()
  
  # Forward selection loop
  while(length(remaining_terms) > 0) {
    candidate_models <- list()
    
    for(term in remaining_terms) {
      new_formula <- as.formula(
        paste(response, "~", paste(c(selected_terms, term), collapse = " +"),
              "+ (1 |", random_term, ")")
      )
      candidate_models[[term]] <- lmer(new_formula, data = data, REML = FALSE)
    }
    
    # Compare AIC
    aic_values <- sapply(candidate_models, AIC)
    best_term <- names(which.min(aic_values))
    
    # Compare with current model using LRT
    best_model <- candidate_models[[best_term]]
    lrt <- anova(current_model, best_model)
    
    if(lrt$`Pr(>Chisq)`[2] < 0.05) {
      # Accept term
      selected_terms <- c(selected_terms, best_term)
      current_model <- best_model
      remaining_terms <- setdiff(remaining_terms, best_term)
    } else {
      break
    }
  }
  
  # Final model summary
  mod_sum <- summary(current_model)
  mod_coeff <- as.data.frame(mod_sum$coefficients) %>%
    rownames_to_column("coeff_term") %>%
    mutate(coeff_term = str_replace(coeff_term, "\\.(L|Q|C)$", "")) %>%
    rename(estimate = Estimate,
           std_error = `Std. Error`,
           t_value = `t value`,
           p_value = `Pr(>|t|)`)  # Now works because lmerTest adds p-values
  
  # ANOVA table
  aov <- Anova(current_model, type = 3)
  aov_df <- as.data.frame(aov) %>%
    rownames_to_column("aov_term") %>%
    mutate(metric = metric_name)
  
  # Merge ANOVA and coefficients
  merged_df <- full_join(aov_df, mod_coeff, by = c("aov_term" = "coeff_term"))
  
  # Compute pseudo-R²
  r2_vals <- r.squaredGLMM(current_model)
  merged_df <- merged_df %>%
    mutate(marginal_R2 = r2_vals[1],
           conditional_R2 = r2_vals[2])
  
  return(list(
    model = current_model,
    merged = merged_df,
    selected_terms = selected_terms,
    pseudo_R2 = r2_vals
  ))
}

# Example usage:
# response_vars <- c("log(Shannon)", "log(InvSimpson)")
# fixed_terms <- c("timepoint", "species", "direct.dry_pod_weight.corr")
# random_term <- "pot_ID"
# results <- lapply(response_vars, function(resp) {
#   forward_select_lmer(resp, fixed_terms, random_term, meta.p, resp)
# })