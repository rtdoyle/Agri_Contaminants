mod_beta.vs.trait_dbRDA_species <- function(df, return_mode = c("both", "combined", "per_time")) {
  
  return_mode <- match.arg(return_mode)
  species_list <- unique(sample_data(df)$species)
  
  # Inner function: dynamic dbRDA with combined + per-time logic
  run_dbRDA_dynamic <- function(df_sub, species_label, return_mode) {
    all_times <- unique(sample_data(df_sub)$timepoint)
    
    # Helper: build dynamic formula and run dbRDA
    run_dbRDA <- function(df_inner, time_label, species_label) {
      dist_bc <- phyloseq::distance(df_inner, method = "bray")
      sample_data_df <- as.data.frame(sample_data(df_inner))
      
      # Select traits
      selected_traits <- c("direct.dry_shoot_weight", "direct.shoot_root_ratio",
                           "direct.per_plant_weight", "direct.shoot_moisture",
                           "direct.survival_perc.corr")
      
      # Add pea-specific traits if only peas present
      species_in_df <- unique(sample_data_df$species)
      if (length(species_in_df) == 1 && species_in_df == "P") {
        selected_traits <- c(selected_traits,
                             "direct.dry_pod_weight.corr",
                             "direct.NDVI_mean",
                             "direct.pods",
                             "direct.flowers")
      }
      
      # Subset trait data
      trait_data_subset <- sample_data_df[, selected_traits]
      
      # Impute missing values
      trait_data_subset[] <- lapply(trait_data_subset, function(x) {
        if (is.numeric(x)) ifelse(is.na(x), mean(x, na.rm = TRUE), x) else x
      })
      
      # Add categorical predictors
      trait_data_final <- cbind(trait_data_subset,
                                sample_data_df[, c("treat.rn", "species", "timepoint")])
      
      # Dynamic formula
      numeric_predictors <- paste(selected_traits, collapse = " + ")
      terms <- c("treat.rn")
      if (length(unique(sample_data_df$species)) > 1) terms <- c(terms, "species")
      if (length(unique(sample_data_df$timepoint)) > 1) terms <- c(terms, "timepoint")
      interaction_terms <- paste(terms, collapse = " * ")
      full_formula <- as.formula(paste("dist_bc ~", interaction_terms, "+", numeric_predictors))
      
      # Run dbRDA
      dbRDA_model <- capscale(full_formula, data = trait_data_final)
      
      # Axis labels
      eig_vals <- dbRDA_model$CCA$eig
      var_explained <- eig_vals / sum(eig_vals) * 100
      x_lab <- paste0("CAP1 (", round(var_explained[1], 1), "%)")
      y_lab <- paste0("CAP2 (", round(var_explained[2], 1), "%)")
      
      # Fit trait vectors
      fit <- envfit(dbRDA_model, trait_data_subset, permutations = 999)
      trait_vectors <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
      trait_vectors$Trait <- rownames(trait_vectors)
      trait_vectors$pval <- fit$vectors$pvals
      sig_traits <- trait_vectors$pval < 0.05
      trait_vectors_sig <- trait_vectors[sig_traits, ]
      
      # Shortened labels
      trait_labels <- c(
        "direct.dry_shoot_weight" = "Shoot dry wt",
        "direct.shoot_root_ratio" = "Shoot:root",
        "direct.per_plant_weight" = "Plant wt",
        "direct.shoot_moisture" = "Shoot moist",
        "direct.survival_perc.corr" = "Survival %",
        "direct.dry_pod_weight.corr" = "Pod dry wt",
        "direct.flowers" = "Flowers",
        "direct.pods" = "Pods"
      )
      trait_vectors_sig$Trait <- trait_labels[trait_vectors_sig$Trait]
      
      # Plot
      plot_data <- cbind(scores(dbRDA_model, display = "sites"),
                         trait_data_final[rownames(scores(dbRDA_model, display = "sites")), ])
      plot_data$species <- factor(plot_data$species, levels = c("L", "R", "P"))
      
      p <- ggplot(plot_data, aes(x = CAP1, y = CAP2)) +
        geom_point(aes(colour = treat.rn, shape = species),
                   size = 4, stroke = 2, alpha = 0.7) +
        geom_segment(data = trait_vectors_sig,
                     aes(x = 0, y = 0, xend = CAP1 * 2.5, yend = CAP2 * 2.5),
                     arrow = arrow(length = unit(0.2, "cm")), color = "black") +
        geom_text_repel(data = trait_vectors_sig,
                        aes(x = CAP1 * 2.5, y = CAP2 * 2.5, label = Trait),
                        color = "black", size = 4, max.overlaps = 100) +
        scale_color_manual(values = c("gray", "lightblue", "cornflowerblue", "blue",
                                      "#EABD8C", "#FFAD00", "#B06500")) +
        scale_shape_manual(values = c("L" = 16, "R" = 15, "P" = 17),
                           labels = c("Lettuce", "Radish", "Pea")) +
        labs(color = "Amendment-\nSpiking level", shape = "Crop Species",
             x = x_lab, y = y_lab) +
        theme_bw() +
        theme(axis.title.x = element_text(size=14), 
      axis.title.y = element_text(size=14), 
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size= 12), 
      legend.text = element_text(size = 10), 
      legend.title = element_text(size =10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
      
      ## ---- Results Formatting ----
      # Continuous trait results
      vector_df <- data.frame(
        Term = rownames(fit$vectors$arrows),
        R_squared = fit$vectors$r,
        P_value = fit$vectors$pvals,
        Time = time_label,
        Species = species_label
      )
      
      # ANOVA tables
      aov_all <- anova(dbRDA_model)
      aov_terms <- anova(dbRDA_model, by = "terms")
      
      sum_sqs <- aov_terms$"SumOfSqs"
      total_variation <- sum(sum_sqs)
      prop_explained <- sum_sqs / total_variation
      
      aov_df <- data.frame(
        Term = rownames(aov_terms),
        SumOfSqs = sum_sqs,
        Proportion = prop_explained,
        F_value = aov_terms$`F`,
        Sig = aov_terms$`Pr(>F)`,
        Df = aov_terms$Df,
        Time = time_label,
        Species = species_label
      )
      
      return(list(plot = p,
                  model = dbRDA_model,
                  anova_terms = aov_terms,
                  vector_df = vector_df,
                  aov_df = aov_df))
    }
    
    # Results container
    results <- list()
    
    if (return_mode %in% c("combined", "both")) {
      results$combined <- run_dbRDA(df_sub, time_label = "all_timepoints", species_label = species_label)
    }
    
    if (return_mode %in% c("per_time", "both")) {
      per_time_results <- lapply(all_times, function(tp) {
        df_tp <- prune_samples(sample_data(df_sub)$timepoint == tp, df_sub)
        run_dbRDA(df_tp, time_label = tp, species_label = species_label)
      })
      names(per_time_results) <- all_times
      results$per_time <- per_time_results
    }
    
    return(results)
  }
  
  # Species-level iteration
  species_results <- lapply(species_list, function(sp) {
    df_sp <- prune_samples(sample_data(df)$species == sp, df)
    run_dbRDA_dynamic(df_sp, species_label = sp, return_mode)
  })
  names(species_results) <- species_list
  
  # Add all species combined
  species_results[["all_species"]] <- run_dbRDA_dynamic(df, species_label = "all_species", return_mode)
  
  return(species_results)
}