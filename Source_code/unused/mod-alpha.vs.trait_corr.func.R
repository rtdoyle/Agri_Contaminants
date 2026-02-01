run_correlations <- function(data,
                             metric.list, 
                             trait.list, 
                             species.list, 
                             output_dir = "./16S_analyses_files/figs/corr_plots") {
  results <- list()
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (t in trait.list) {
    
    # Determine which species to use based on trait
    current_species <- if (t == "direct.survival_perc.corr") {
      species.list
    } else {
      "P"
    }
    
    for (species in current_species) {
      species_data <- data %>% filter(species == !!species)
      
      for (m in metric.list) {
        plot_df <- species_data %>%
          select(all_of(c(m, t, "treat.rn"))) %>%
          rename(metric = !!m, trait = !!t, treatment = treat.rn)
        
        model <- lm(trait ~ log(metric), data = plot_df)
        tidy_model <- broom::tidy(model)
        r_squared <- summary(model)$r.squared
        corr <- cor.test(plot_df$metric, plot_df$trait)
        
        result_name <- paste(species, m, t, sep = "_vs_")
        results[[result_name]] <- data.frame(
          Species = species,
          Diversity = m,
          Trait = t,
          R_squared = r_squared,
          P_value_regression = tidy_model$p.value[2],
          Pearson_corr = corr$estimate,
          P_value_correlation = corr$p.value
        )
        
        # Plot if correlation is marginally significant
        if (corr$p.value < 0.1) {
          p <- ggplot(plot_df, aes(x = log(metric), y = trait)) +
            geom_hline(aes(yintercept = 0), linetype = 2) +
            geom_vline(aes(xintercept = 0), linetype = 2) +
            geom_point(alpha = 0.6, aes(colour = treatment)) +
            geom_smooth(method = "lm", se = TRUE, color = "blue") +
            labs(title = paste("Correlation:", m, "vs", t, "in", species),
                 subtitle = paste0("R² = ", round(r_squared, 3),
                                   ", r = ", round(corr$estimate, 3),
                                   ", p = ", signif(corr$p.value, 3))) +
            scale_colour_manual(name = "Amendment- \n Spiking Level",
                                values = c("gray", "lightblue", "cornflowerblue",
                                           "blue", "#EABD8C", "#FFAD00", "#B06500"),
                                guide = guide_legend(nrow = 4, ncol = 2)) +
            theme_minimal()
          
          ggsave(filename = file.path(output_dir, paste0(species, "_", m, "_vs_", t, ".png")),
                 plot = p, width = 6, height = 4)
        }
      }
    }
  }
  
  return(results)
}