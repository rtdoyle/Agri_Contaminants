plot_beta_trait_NMDS <- function(combs, times, crops, df){
  
  ## timepoints
  print(combs)
  
  ## subset to peas for TP1
  df.f <- df %>%
    filter(species == crops &
             timepoint == times) %>%
    droplevels(.)
  
  samples <- unique(df.f$Study_ID)
  
  otu_subset <- otu_table[rownames(otu_table) %in% samples, ]
  
  ## normalize
  rel_abundance <- decostand(otu_subset, method = "total")
  
  # Perform ordination (e.g., NMDS or PCoA)
  nmds <- metaMDS(rel_abundance, k = 2, trymax = 100)
  plot(nmds)
  
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
        "indirect.shoot_moisture")
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
  
  
  ## run through 
  env_data <- df.f[, c(selected_traits, "treat.rn")]
  
  ## to get consistent results
  set.seed(123)  # Use any fixed number
  
  fit <- envfit(nmds, env_data, permutations = 999,
                na.rm = TRUE)
  
  # NMDS sample scores
  nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
  nmds_scores$Study_ID <- rownames(nmds_scores)
  
  # Trait vectors
  trait_vectors <- as.data.frame(scores(fit, display = "vectors"))
  trait_vectors$Trait <- rownames(trait_vectors)
  
  ## filter for sig traits
  sig_traits <- fit$vectors$pvals < 0.05
  trait_vectors_sig <- trait_vectors[sig_traits, ]
  
  ## rename traits
  trait_vectors_sig$Trait <- trait_labels[trait_vectors_sig$Trait]
  
  ## scale segs
  scale_factor <- 1
  
  ## plot
  p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = df.f$treat.rn), size = 3) +  # Optional grouping
    geom_segment(data = trait_vectors_sig,
                 aes(x = 0, y = 0, xend = NMDS1 * scale_factor, 
                     yend = NMDS2 * scale_factor),
                 arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
    geom_text(data = trait_vectors_sig,
              aes(x = NMDS1 * scale_factor, 
                  y = NMDS2 * scale_factor, label = Trait),
              color = "blue", vjust = -0.5) +
    # additional formatting
    scale_color_manual(values = c("gray",
                                  "lightblue",
                                  "cornflowerblue",
                                  "blue",
                                  "#EABD8C",
                                  "#FFAD00",
                                  "#B06500")) +
    labs(color = "Amendment- \n Spiking level") +
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
    Species = paste0(crops),   # assuming crop is a single value
    Time = paste0(times)       # assuming time is a single value
  )
  
  # Extract categorical trait results (factors)
  factor_df <- data.frame(
    Term = names(fit$factors$pvals),
    R_squared = NA,  # Not applicable for factors
    P_value = fit$factors$pvals,
    Species = paste0(crops),
    Time = paste0(times)
  )
  
  # Combine both
  results_df <- rbind(vector_df, factor_df)
  
  return(list(p, results_df))
}