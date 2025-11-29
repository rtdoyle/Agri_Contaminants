
mod_beta.vs.trait_dbRDA_species <- function(df, return_mode = c("both", "combined", "per_time")) {
  
  return_mode <- match.arg(return_mode)
  species_list <- unique(sample_data(df)$species)
  
  # Inner function: dynamic dbRDA with combined + per-time logic
  run_dbRDA_dynamic <- function(df_sub, species_label, return_mode) {
    all_times <- unique(sample_data(df_sub)$timepoint)
    
    # Helper: build dynamic formula and run dbRDA with forward selection
    run_dbRDA <- function(df_inner, time_label, species_label) {
      dist_bc <- phyloseq::distance(df_inner, method = "bray")
      sample_data_df <- as.data.frame(sample_data(df_inner))
      
      # Select traits
      selected_traits <- c(
        "direct.dry_total_weight", "direct.shoot_root_ratio",
        "direct.per_plant_weight", "direct.shoot_moisture",
        "direct.survival_perc.corr"
      )
      
      # Add pea-specific traits if only peas present
      species_in_df <- unique(sample_data_df$species)
      if (length(species_in_df) == 1 && species_in_df == "P") {
        selected_traits <- c(
          selected_traits,
          "direct.dry_pod_weight.corr",
          "direct.pods",
          "direct.flowers"
        )
      }
      
      # Subset trait data (numeric predictors)
      trait_data_subset <- sample_data_df[, selected_traits, drop = FALSE]
      
      # Impute missing values
      trait_data_subset[] <- lapply(trait_data_subset, function(x) {
        if (is.numeric(x)) ifelse(is.na(x), mean(x, na.rm = TRUE), x) else x
      })
      
      # Add categorical predictors (design factors)
      trait_data_final <- cbind(
        trait_data_subset,
        sample_data_df[, c("treat.rn", "species", "timepoint")]
      )
      
      # Dynamic formula components
      numeric_predictors <- paste(selected_traits, collapse = " + ")
      terms <- c("treat.rn")
      if (length(unique(sample_data_df$species)) > 1) terms <- c(terms, "species")
      if (length(unique(sample_data_df$timepoint)) > 1) terms <- c(terms, "timepoint")
      interaction_terms <- paste(terms, collapse = " *")
      interaction_terms <- gsub("\\*\\*", "*", interaction_terms)  # just in case
      full_formula <- as.formula(paste("dist_bc ~", interaction_terms, "+", numeric_predictors))
      
      # ---- Forward Selection ----
      null_model <- vegan::capscale(dist_bc ~ 1, data = trait_data_final)
      full_model <- vegan::capscale(full_formula, data = trait_data_final)
      
      forward_model <- vegan::ordistep(
        null_model, scope = formula(full_model),
        direction = "forward", permutations = 999
      )
      
      dbRDA_model <- forward_model
      
      # >>> Identify which trait terms were retained by forward selection
      #     We read the final formula and get its term labels; keep only the numeric trait names.
      final_formula <- formula(dbRDA_model)
      final_terms <- attr(terms(final_formula), "term.labels")
      selected_trait_terms <- intersect(selected_traits, final_terms)
      
      # >>> Use only forward-selected traits for envfit
      trait_data_envfit <- if (length(selected_trait_terms) > 0) {
        trait_data_subset[, selected_trait_terms, drop = FALSE]
      } else {
        # If no traits selected, make a 0-row/0-col data.frame so envfit can be skipped cleanly
        data.frame()
      }
      
      # Axis labels
      eig_vals <- dbRDA_model$CCA$eig
      var_explained <- eig_vals / sum(eig_vals) * 100
      x_lab <- paste0("CAP1 (", round(var_explained[1], 1), "%)")
      y_lab <- paste0("CAP2 (", round(var_explained[2], 1), "%)")
      
      # Fit trait vectors (only on significant forward-selected traits)
      if (ncol(trait_data_envfit) > 0) {
        fit <- vegan::envfit(dbRDA_model, trait_data_envfit, permutations = 999)
        trait_vectors <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
        trait_vectors$Trait <- rownames(trait_vectors)
        trait_vectors$pval <- fit$vectors$pvals
        sig_traits <- trait_vectors$pval < 0.05
        trait_vectors_sig <- trait_vectors[sig_traits, ]
      } else {
        # No traits retained by forward step -> no arrows
        fit <- NULL
        trait_vectors_sig <- data.frame(CAP1 = numeric(0), CAP2 = numeric(0), Trait = character(0), pval = numeric(0))
      }
      
      # Shortened labels
      trait_labels <- c(
        "direct.dry_shoot_weight"   = "Shoot dry wt",
        "direct.shoot_root_ratio"   = "Shoot:root",
        "direct.per_plant_weight"   = "Plant wt",
        "direct.shoot_moisture"     = "Shoot moist",
        "direct.survival_perc.corr" = "Survival %",
        "direct.dry_pod_weight.corr"= "Pod dry wt",
        "direct.flowers"            = "Flowers",
        "direct.pods"               = "Pods",
        "direct.NDVI_mean"          = "NDVI"
      )
      if (nrow(trait_vectors_sig) > 0) {
        trait_vectors_sig$Trait <- trait_labels[trait_vectors_sig$Trait]
      }
      
      # Plot
      plot_scores <- vegan::scores(dbRDA_model, display = "sites")
      plot_data <- cbind(plot_scores, trait_data_final[rownames(plot_scores), ])
      plot_data$species <- factor(plot_data$species, levels = c("L", "R", "P"))
      
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = CAP1, y = CAP2)) +
        ggplot2::geom_point(ggplot2::aes(colour = treat.rn, shape = species),
                            size = 4, stroke = 2, alpha = 0.7) +
        ggplot2::geom_segment(data = trait_vectors_sig,
                              ggplot2::aes(x = 0, y = 0, xend = CAP1 * 2.5, yend = CAP2 * 2.5),
                              arrow = grid::arrow(length = grid::unit(0.2, "cm")), color = "black") +
        ggrepel::geom_text_repel(data = trait_vectors_sig,
                                 ggplot2::aes(x = CAP1 * 2.5, y = CAP2 * 2.5, label = Trait),
                                 color = "black", size = 4, max.overlaps = 100) +
        ggplot2::scale_color_manual(values = c("gray", "lightblue", "cornflowerblue", "blue",
                                               "#EABD8C", "#FFAD00", "#B06500")) +
        ggplot2::scale_shape_manual(values = c("L" = 16, "R" = 15, "P" = 17),
                                    labels = c("Lettuce", "Radish", "Pea")) +
        ggplot2::labs(color = "Amendment-\nSpiking level", shape = "Crop Species",
                      x = x_lab, y = y_lab) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size=14), 
                       axis.title.y = ggplot2::element_text(size=14), 
                       axis.text.x = ggplot2::element_text(size = 12), 
                       axis.text.y = ggplot2::element_text(size= 12), 
                       legend.text = ggplot2::element_text(size = 10), 
                       legend.title = ggplot2::element_text(size =10),
                       plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"), 
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank())
      
      ## ---- Results Formatting ----
      if (!is.null(fit) && nrow(trait_vectors_sig) + (nrow(trait_vectors_sig) == 0) >= 0) {
        vector_df <- data.frame(
          Term = rownames(fit$vectors$arrows),
          R_squared = fit$vectors$r,
          P_value = fit$vectors$pvals,
          Time = time_label,
          Species = species_label
        )
      } else {
        vector_df <- data.frame(
          Term = character(0),
          R_squared = numeric(0),
          P_value = numeric(0),
          Time = character(0),
          Species = character(0)
        )
      }
      
      aov_terms <- stats::anova(dbRDA_model, by = "terms")
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
                  aov_df = aov_df,
                  selected_trait_terms = selected_trait_terms))  # >>> expose which traits were used
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
