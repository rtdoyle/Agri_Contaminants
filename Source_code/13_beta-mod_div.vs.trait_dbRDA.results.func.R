flatten_dbRDA_results <- function(results_list) {
  
  # Initialize empty lists
  all_aov <- list()
  all_vectors <- list()
  
  # Loop through species
  for (species_name in names(results_list)) {
    species_data <- results_list[[species_name]]
    
    # Combined model
    if (!is.null(species_data$combined)) {
      all_aov[[paste0(species_name, "_combined")]] <- species_data$combined$aov_df
      all_vectors[[paste0(species_name, "_combined")]] <- species_data$combined$vector_df
    }
    
    # Per-time models
    if (!is.null(species_data$per_time)) {
      for (tp in names(species_data$per_time)) {
        all_aov[[paste0(species_name, "_", tp)]] <- species_data$per_time[[tp]]$aov_df
        all_vectors[[paste0(species_name, "_", tp)]] <- species_data$per_time[[tp]]$vector_df
      }
    }
  }
  
  # Combine into data frames
  aov_df_all <- do.call(rbind, all_aov)
  vector_df_all <- do.call(rbind, all_vectors)
  
  return(list(aov_df_all = aov_df_all, vector_df_all = vector_df_all))
}