compute_and_plot_correlations <- function(data, metrics, trait_sets, fig_dir, mod_dir, 
                                          save_csv = TRUE) {
  
  # Trait label lookup
  trait_labels <- c(
    "direct.survival_perc.corr" = "Survival (%)",
    "direct.shoot_root_ratio" = "Shoot:root",
    "direct.dry_total_weight" = "Total weight",
    "direct.per_plant_weight" = "Per plant weight",
    "direct.shoot_moisture" = "Shoot moisture",
    "direct.germination_perc" = "Germination",
    "direct.dry_pod_weight.corr" = "Dry pod weight",
    "direct.per_pod_weight" = "Per pod weight",
    "direct.infected_peas_perc.corr" = "Infection",
    "direct.pods" = "Pods",
    "direct.flowers" = "Flowers",
    "indirect.plant_height1" = "Plant height (wk1)",             
    "indirect.plant_height2" = "Plant height (wk2)",            
    "indirect.plant_height3" = "Plant height (wk3)",             
    "indirect.chlorophyll_content1_mean" = "Leaf Chloro (wk1)",
    "indirect.chlorophyll_content2_mean" = "Leaf Chloro (wk2)",
    "indirect.chlorophyll_content3_mean" = "Leaf Chloro (wk3)",
    "indirect.leaf_num" = "Leaves (wk2)",
    "indirect.plant_height4" = "Plant height (wk4)",            
    "indirect.chlorophyll_content4_mean" = "Leaf Chloro (wk4)",
    "indirect.leaf_num2" = "Leaves (wk4)",               
    "indirect.wet_shoot_weight" = "Wet shoot weight",
    "indirect.dry_shoot_weight" = "Dry shoot weight",         
    "indirect.wet_root_weight" = "Wet root weight",           
    "indirect.dry_root_weight" = "Dry root weight",          
    "indirect.dry_total_weight" = "Dry total weight",
    "indirect.shoot_root_ratio" = "Shoot:root",         
    "indirect.shoot_moisture" = "Shoot moisture"        
  )
  
  metric_labs <- c(
    Shannon = "Shannon",
    InvSimpson = "Inv. Simpson (/100)"
  )
  
  # Use map() to return results
  summary_list <- map(names(trait_sets), function(set_name) {
    traits <- trait_sets[[set_name]]
    
    # Filter data if needed
    data_subset <- if (set_name %in% c("pea", "indirect")) {
      data %>% filter(species == "P")
    } else {
      data
    }
    
    # Compute correlations
    cor_df <- map_df(metrics, function(metric) {
      map_df(unique(data_subset$species), function(sp) {
        map_df(unique(data_subset$timepoint), function(tp) {
          subset_df <- data_subset %>%
            filter(species == sp, timepoint == tp)
          
          if (nrow(subset_df) == 0) return(NULL)
          
          map_df(traits, function(trait) {
            cor_test <- cor.test(subset_df[[metric]], subset_df[[trait]], method = "pearson")
            tibble(
              metric = metric,
              species = sp,
              timepoint = tp,
              trait = trait,
              r = cor_test$estimate,
              p_value = cor_test$p.value
            )
          })
        })
      })
    }) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        p_cat = case_when(
          p_adj < 0.001 ~ "<0.001",
          p_adj < 0.01  ~ "<0.01",
          p_adj < 0.05  ~ "<0.05",
          p_adj < 0.1   ~ "<0.1",
          TRUE ~ "ns"
        ),
        tp_species = paste(timepoint, species, sep = "_"),
        trait_abbrev = dplyr::recode(trait, !!!trait_labels)
      )
    
    # Save summary table as CSV if requested
    if (save_csv) {
      write.csv(cor_df, paste0(mod_dir, "corr-alpha.vs.trait_", set_name, ".csv"), row.names = FALSE)
    }
    
    # Ordering for heatmap
    species_levels <- unique(cor_df$species)
    time_levels <- c("TP1", "TP2")
    ordered_cols <- as.vector(outer(time_levels, species_levels, paste, sep = "_"))
    cor_df$tp_species <- factor(cor_df$tp_species, levels = ordered_cols)
    
    # Heatmap
    htmap <- ggplot(cor_df, aes(x = tp_species, y = trait_abbrev, fill = r)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      geom_text(aes(label = ifelse(p_adj < 0.1, p_cat, "")), size = 3) +
      facet_wrap(~ metric, scales = "free_x",
                 labeller = labeller(metric = metric_labs)) +
      labs(x = "Timepoint-Species", y = "Trait", fill = "Pearson r") +
      theme_bw(base_size = 14) +
      theme(axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
            axis.title.y = element_text(size = 16),
            axis.text.y = element_text(size = 14),
            strip.text = element_text(size = 14, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # Save plot
    ggsave(paste0(fig_dir, "heatmap_", set_name, ".png"),
           width = 8, height = 4, plot = htmap, units = "in")
    
    message("Saved heatmap and summary for trait set: ", set_name)
    
    return(cor_df)  # Return summary table for this set
  })
  
  names(summary_list) <- names(trait_sets)
  return(summary_list)
}