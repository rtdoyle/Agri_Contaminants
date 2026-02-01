mod_beta.vs.trait_dbRDA.pea <- function(times, df){
  
  ## combs
  print(times)
  
  ## set to characters
  time <- as.character(times)
  
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
  
  ## selected traits
  selected_traits <- c(
    "direct.shoot_root_ratio", 
    "direct.per_plant_weight",
    "direct.shoot_moisture",   
    "direct.survival_perc.corr",
    "direct.dry_pod_weight.corr",
    "direct.NDVI_mean",
    "direct.pods",
    "direct.flowers")
  
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
                                               "timepoint",
                                               "species")])
  
  # Create a formula string with numeric predictors
  numeric_predictors <- paste(selected_traits, collapse = " + ")
  
  # Combine with categorical predictors and interactions
  full_formula <- as.formula(
    paste("dist_bc ~ treat.rn + ", numeric_predictors)
  )
  
  # Step 4: Run dbRDA using capscale with all predictors
  dbRDA_model <- capscale(full_formula, data = trait_data_final)
  
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
  
  ## rename traits
  trait_vectors_sig$Trait <- trait_labels[trait_vectors_sig$Trait]
  
  ## scale segs
  scale_factor <- 2.5
  
  ## plot
  p <- ggplot(plot_data, aes(x = CAP1, y = CAP2)) +
    geom_point(aes(color = treat.rn,
                   shape = species), size = 3) +  # Optional grouping
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
    scale_shape_manual(values = c(17)) +
    labs(color = "Amendment- \n Spiking level",
         shape = "Crop species") +
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
      legend.title = element_text(size =10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Extract continuous trait results (vectors)
  vector_df <- data.frame(
    Term = rownames(fit$vectors$arrows),
    R_squared = fit$vectors$r,
    P_value = fit$vectors$pvals,
    Time = time       # assuming time is a single value
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
    Time = time
  )
  
  return(list(p, vector_df, aov_df))
}