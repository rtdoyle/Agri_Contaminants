capscale.func <- function(times, ps, df){
  
  ## print time-measure comb
  print(times)
  
  ## set to characters
  time <- as.character(times)
  
  ## subset to specific timepoints
  ps.f <- prune_samples(sample_data(ps)$timepoint %in% time, ps)
  
  ## calculate Bray Curtis distnaces
  dist_bc <- phyloseq::distance(ps.f, method = "bray")
  
  ## filter meta to times
  df.f <- df %>%
    filter(timepoint %in% times) %>%
    droplevels(.)
  
  # Run capscale
  cap_result <- capscale(dist_bc ~ 
                           treat.rn + species, 
                         data = df.f)
  
  # Total constrained inertia
  (total_constrained <- cap_result$CCA$tot.chi)
  
  # Run ANOVA by term
  anova_terms <- anova(cap_result, by = "term")
  
  # Extract Sum of Squares
  sum_sqs <- anova_terms$"SumOfSqs"
  
  # Calculate total variation (including residuals)
  total_variation <- sum(sum_sqs)
  
  # Calculate proportion explained
  prop_explained <- sum_sqs / total_variation
  
  # Combine into a data frame
  explained_df <- data.frame(
    Term = rownames(anova_terms),
    SumOfSqs = sum_sqs,
    Proportion = prop_explained
  )
  
  print(explained_df)
  
  PVE <- as.data.frame(explained_df)
  PVE$timepoint <- paste0(time)
  
  return(list(PVE))
}