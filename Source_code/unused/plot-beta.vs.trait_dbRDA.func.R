plot_beta_trait_dbRDA <- function(combs, times, crops, df){
  
  ## combs
  print(combs)
  
  ## set to characters
  time <- as.character(times)
  crop <- as.character(crops)
  
  ## subset to specific timepoints and species
  df.f <- prune_samples(sample_data(df)$timepoint %in% time, df)
  df.ff <- prune_samples(sample_data(df.f)$species %in% crop, df.f)
  
  # Step 1: Calculate Bray-Curtis distance from phyloseq object
  dist_bc <- phyloseq::distance(df.ff, method = "bray")
  
  # Step 2: Extract OTU table and sample data
  otu_table_df <- as.data.frame(phyloseq::otu_table(df.ff))
  if (phyloseq::taxa_are_rows(df.ff)) {
    otu_table_df <- t(otu_table_df)
  }
  
  sample_data_df <- as.data.frame(phyloseq::sample_data(df.ff))
  
  if (crops %in% c("L", "R")){
  
  ## select only traits of interest
  selected_traits <- c(
    "direct.dry_total_weight",
    "direct.shoot_root_ratio",
    "direct.survival_perc.corr")
  }

  else if (crops %in% c("P")){
    
    ## select only traits of interest
    selected_traits <- c(
      "direct.dry_total_weight",
      "direct.shoot_root_ratio",
      "direct.survival_perc.corr",
      ## pea only traits
      "direct.dry_pod_weight.corr",
      "direct.infected_peas_perc.corr",
      "indirect.chlorophyll_content4_mean",
      "indirect.leaf_num2",
      "indirect.dry_total_weight",
      "indirect.shoot_root_ratio",
      "indirect.plant_height4",
      "indirect.shoot_moisture"
    )
  }
  
  ## shortened labels
  trait_labels <- c(
    "direct.dry_shoot_weight" = "Shoot dry wt",
    "direct.dry_total_weight" = "Total dry wt",
    "direct.per_plant_weight" = "Plant wt",
    "direct.shoot_moisture" = "Shoot moist (dir)",
    "direct.survival_perc.corr" = "Survival %",
    "direct.germination_perc" = "Germination %",
    "direct.shoot_root_ratio" = "Shoot:root (dir)",
    "direct.dry_pod_weight.corr" = "Pod dry wt",
    "direct.flowers" = "Flowers",
    "direct.pods" = "Pods",
    "direct.infected_peas_perc.corr" = "Infected %",
    "indirect.plant_height4" = "Height",
    "indirect.shoot_root_ratio" = "Shoot:root (ind)",
    "indirect.shoot_moisture" = "Shoot moist (ind)",
    "indirect.dry_total_weight" = "Total dry wt",
    "indirect.chlorophyll_content4_mean" = "NVDI (ind)",
    "indirect.leaf_num2" = "Leaves (no.)"
  )
  
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
  trait_data_final <- cbind(trait_data_subset, sample_data_df[, c("treat.rn")])
  
  # Step 4: Run dbRDA using capscale with all predictors
  dbRDA_model <- capscale(dist_bc ~ ., data = trait_data_final)
  
  # Extract eigenvalues from capscale model
  eig_vals <- dbRDA_model$CCA$eig  # eigenvalues for constrained axes
  var_explained <- eig_vals / sum(eig_vals) * 100  # percent variation
  print(var_explained)

  # Format axis labels with percent explained (rounded to 1 decimal)
  x_lab <- paste0("CAP1 (", round(var_explained[1], 1), "%)")
  y_lab <- paste0("CAP2 (", round(var_explained[2], 1), "%)")
  
  ## to get consistent results
  set.seed(123)  # Use any fixed number
  
  # Step 5: Fit trait vectors using envfit
  fit <- envfit(dbRDA_model, trait_data_final, permutations = 999)
  
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
  
  ## rename traits
  trait_vectors_sig$Trait <- trait_labels[trait_vectors_sig$Trait]
  
  ## scale segs
  scale_factor <- 1
  
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
    geom_text(data = trait_vectors_sig,
              aes(x = CAP1 * scale_factor, 
                  y = CAP2 * scale_factor, label = Trait),
              color = "blue", vjust = -0.5) +
    # additional formatting
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    labs(color = "Amendment- \n Spiking level",
         x = x_lab,
         y = y_lab) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
      axis.title.x = element_text(size=14), 
      axis.title.y = element_text(size=14), 
      axis.text.x = element_text(size = 12), 
      axis.text.y = element_text(size= 12), 
      legend.text = element_text(size = 10), 
      legend.title = element_text(size =10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Extract continuous trait results (vectors)
  vector_df <- data.frame(
    Term = rownames(fit$vectors$arrows),
    R_squared = fit$vectors$r,
    P_value = fit$vectors$pvals,
    Species = crop,   # assuming crop is a single value
    Time = time       # assuming time is a single value
  )
  
  # Extract categorical trait results (factors)
  factor_df <- data.frame(
    Term = names(fit$factors$pvals),
    R_squared = NA,  # Not applicable for factors
    P_value = fit$factors$pvals,
    Species = crop,
    Time = time
  )
  
  # Combine both
  results_df <- rbind(vector_df, factor_df)
  
  
  ## ANOVA
  # Overall model significance
  aov_all <- anova(dbRDA_model)
  print(aov_all)
  
  # Significance of each predictor
  aov_terms <- anova(dbRDA_model, by = "terms")
  
  # Extract Sum of Squares
  sum_sqs <- aov_terms$"SumOfSqs"
  total_variation <- sum(sum_sqs)
  prop_explained <- sum_sqs / total_variation
  
  # Combine into a data frame
  explained_df <- data.frame(
    Term = rownames(aov_terms),
    SumOfSqs = sum_sqs,
    Proportion = prop_explained,
    Time = time,
    Species = crop
  )
  
  return(list(p, results_df, explained_df))
}