mod_beta.vs.trait_dbRDA <- function(times, df, 
                                    include_species_term = TRUE){
  
  ## combs
  print(times)
  
  ## set to characters
  time <- as.character(times)
  
  ## detect species (for determining which traits to include)
  species_in_df <- unique(sample_data(df)$species)
  print(species_in_df)
  
  ## subset to specific timepoints and species
  df.f <- prune_samples(sample_data(df)$timepoint %in% time, df)
  
  # Step 1: Calculate Bray-Curtis distance from phyloseq object
  dist_bc <- phyloseq::distance(df.f, method = "bray")
  
  # Step 2: Extract OTU table and sample data
  otu_table_df <- as.data.frame(phyloseq::otu_table(df.f))
  if (phyloseq::taxa_are_rows(df.f)) {
    otu_table_df <- t(otu_table_df)
  }
  sample_data_df <- as.data.frame(phyloseq::sample_data(df.f))
  
  selected_traits <- c(
    "direct.dry_shoot_weight",   
    "direct.shoot_root_ratio",    
    "direct.per_plant_weight",   
    "direct.shoot_moisture",      
    "direct.survival_perc.corr"
  )
  
  # Add pea-specific traits if only peas are present
  if (length(species_in_df) == 1 && species_in_df == "P") {
    selected_traits <- c(selected_traits, 
                         "direct.dry_pod_weight.corr",
                         "direct.NDVI_mean",
                         "direct.pods",
                         "direct.flowers")
  }
  
  # Subset trait data to only those columns
  trait_data_subset <- sample_data_df[, selected_traits]
  
  # Impute numeric columns only
  trait_data_subset[] <- lapply(trait_data_subset, function(x) {
    if (is.numeric(x)) {
      ifelse(is.na(x), mean(x, na.rm = TRUE), x)
    } else {
      x  # leave non-numeric columns unchanged
    }
  })
  
  # Add categorical variables
  trait_data_final <- cbind(trait_data_subset, 
                            sample_data_df[, c("treat.rn",
                                               "species",
                                               "timepoint",
                                               "species_tp")])
  
  # Create a formula string with numeric predictors
  numeric_predictors <- paste(selected_traits, collapse = " + ")
  
  # Combine with categorical predictors and interactions
  if (include_species_term) {
    full_formula <- as.formula(
      paste("dist_bc ~ (species * treat.rn) +", numeric_predictors)
    )
  } else {
    full_formula <- as.formula(
      paste("dist_bc ~ treat.rn +", numeric_predictors)
    )
  }
  
  # Step 4: Run dbRDA using capscale with all predictors
  dbRDA_model <- capscale(full_formula, data = trait_data_final)
  
  # Extract eigenvalues from capscale model
  eig_vals <- dbRDA_model$CCA$eig  # eigenvalues for constrained axes
  var_explained <- eig_vals / sum(eig_vals) * 100  # percent variation
  
  # Format axis labels with percent explained (rounded to 1 decimal)
  x_lab <- paste0("CAP1 (", round(var_explained[1], 1), "%)")
  y_lab <- paste0("CAP2 (", round(var_explained[2], 1), "%)")
  
  ## to get consistent results
  set.seed(123)  # Use any fixed number
  
  # Step 5: Fit trait vectors using envfit
  fit <- envfit(dbRDA_model, trait_data_subset, permutations = 999)
  
  # Step 6: Extract site scores and merge with sample data
  site_scores <- scores(dbRDA_model, display = "sites")
  site_df <- as.data.frame(site_scores)
  site_df$SampleID <- rownames(site_df)
  
  plot_data <- cbind(site_df, trait_data_final[rownames(site_df), ])
  
  # Step 7: Extract trait vectors
  trait_vectors <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
  trait_vectors$Trait <- rownames(trait_vectors)
  trait_vectors$pval <- fit$vectors$pvals
  
  ## filter for sig traits
  sig_traits <- fit$vectors$pvals < 0.05
  trait_vectors_sig <- trait_vectors[sig_traits, ]
  
  ## shortened labels
  trait_labels <- c(
    "direct.dry_shoot_weight" = "Shoot dry wt",
    "direct.dry_root_weight" = "Root dry wt",
    "direct.dry_total_weight" = "Total dry wt",
    "direct.per_plant_weight" = "Plant wt",
    "direct.shoot_moisture" = "Shoot moist",
    "direct.survival_perc.corr" = "Survival %",
    "direct.germination_perc" = "Germination %",
    "direct.shoot_root_ratio" = "Shoot:root",
    "direct.dry_pod_weight.corr" = "Pod dry wt",
    "direct.flowers" = "Flowers",
    "direct.pods" = "Pods",
    "direct.infected_peas_perc.corr" = "Infected %"
  )
  
  ## rename traits
  trait_vectors_sig$Trait <- trait_labels[trait_vectors_sig$Trait]
  
  ## scale segs
  scale_factor <- 2.5
  
  ## set shape aesthetics for species
  plot_data$species <- factor(plot_data$species, 
                              levels = c("L", 
                                         "R", 
                                         "P"))
  
  ## plot
  p <- ggplot(plot_data, aes(x = CAP1, y = CAP2)) +
    geom_point(aes(colour = treat.rn, 
                   shape = species),
               size = 4, 
               stroke = 2,
               alpha = 0.7) +  # Optional grouping
    geom_segment(data = trait_vectors_sig,
                 aes(x = 0, y = 0, xend = CAP1 * scale_factor, 
                     yend = CAP2 * scale_factor),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text_repel(data = trait_vectors_sig,
                    aes(x = CAP1 * scale_factor, 
                        y = CAP2 * scale_factor, label = Trait),
                    color = "blue", size = 4, max.overlaps = 100) +
    # additional formatting
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    scale_shape_manual(values = c("L" = 16,  # circle
                                  "R" = 15,   # square
                                  "P" = 17),
                       labels = c("Lettuce",
                                  "Radish",
                                  "Pea")) +  # triangle
    labs(color = "Amendment- \n Spiking level",
         shape = "Crop Species",
         x = x_lab,
         y = y_lab) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, 
                                size = 16, 
                                face = "bold"), 
      axis.title.x = element_text(size=14), 
      axis.title.y = element_text(size=14), 
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size= 12), 
      legend.text = element_text(size = 10), 
      legend.title = element_text(size =12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Extract continuous trait results (vectors)
  vector_df <- data.frame(
    Term = rownames(fit$vectors$arrows),
    R_squared = fit$vectors$r,
    P_value = fit$vectors$pvals,
    Time = time,       # assuming time is a single value
    Species = paste(species_in_df, collapse = "_")
    )
  
  ## ANOVA
  # Overall model significance
  aov_all <- anova(dbRDA_model)
  
  # Significance of each predictor
  aov_terms <- anova(dbRDA_model, by = "terms")
  
  # Extract Sum of Squares
  sum_sqs <- aov_terms$"SumOfSqs"
  total_variation <- sum(sum_sqs)
  prop_explained <- sum_sqs / total_variation
  
  # Combine into a data frame
  aov_df <- data.frame(
    Term = rownames(aov_terms),
    SumOfSqs = sum_sqs,
    Proportion = prop_explained,
    F_value = aov_terms$`F`,
    Sig = aov_terms$`Pr(>F)`,
    Df = aov_terms$Df,
    Time = time,
    Species = paste(species_in_df, collapse = "_")
    )
  
  return(list(p, vector_df, aov_df))
}