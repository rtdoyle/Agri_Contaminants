plot_cor_heatmaps <- function(data, metrics, trait_sets, fig_dir) {
  
  walk(names(trait_sets), function(set_name) {
    traits <- trait_sets[[set_name]]
    
    # Filter data if trait set is pea or indirect
    data_subset <- if (set_name %in% c("pea", "indirect")) {
      data %>% filter(species == "P")
    } else {
      data
    }
    
    # Create a lookup vector where names match actual trait values
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
      "indirect.plant_height3"= "Plant height (wk3)",             
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
    
    # Compute correlations
    cor_df <- compute_correlations(data_subset, metrics, traits)
    
    # Add labels and ordering
    cor_df <- cor_df %>%
      mutate(tp_species = paste(timepoint, species, sep = "_"),
             trait_abbrev = dplyr::recode(trait, !!!trait_labels))
    
    species_levels <- unique(cor_df$species)
    time_levels <- c("TP1", "TP2")  # adjust if needed
    ordered_cols <- as.vector(outer(time_levels, species_levels, paste, sep = "_"))
    
    cor_df$tp_species <- factor(cor_df$tp_species, levels = ordered_cols)
    
    metric_labs <- c(
      Shannon = "Shannon",
      InvSimpson = "Inv. Simpson (/100)"
    )
    
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
    
    message("Saved heatmap for trait set: ", set_name)
  })
}