run_corr <- function(data, species_col, output_dir, top_n_vars = 5, vars_by_species) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    print(paste("Created output directory:", output_dir))
  } else {
    print(paste("Using existing output directory:", output_dir))
  }
  
  top_vars_by_species <- list()
  skipped_species <- c()  # Track skipped species
  
  for (crop in names(vars_by_species)) {
    print(paste("---- Processing species:", crop, "----"))
    
    subset_data <- if (crop == "ALL") data else data[data[[species_col]] == crop, , drop = FALSE]
    species_vars <- vars_by_species[[crop]]
    
    # Check columns exist
    available_vars <- species_vars[species_vars %in% colnames(subset_data)]
    missing_vars <- setdiff(species_vars, available_vars)
    print(paste("Available vars:", paste(available_vars, collapse = ", ")))
    if (length(missing_vars) > 0) {
      print(paste("Missing vars:", paste(missing_vars, collapse = ", ")))
    }
    
    subset_vars <- subset_data[, available_vars, drop = FALSE]
    species_label <- crop
    
    if (nrow(subset_vars) == 0) {
      print(paste("Skipping", species_label, "- no data available"))
      skipped_species <- c(skipped_species, species_label)
      next
    }
    
    # Replace NAs with 0
    subset_vars[is.na(subset_vars)] <- 0
    print(paste("Rows:", nrow(subset_vars), "Cols:", ncol(subset_vars)))
    
    # Ensure numeric and preserve names
    subset_vars <- as.data.frame(lapply(subset_vars, as.numeric))
    rownames(subset_vars) <- rownames(subset_data)  # keep sample IDs
    colnames(subset_vars) <- available_vars         # keep variable names
    
    # Remove zero-variance columns
    zero_var_cols <- names(which(apply(subset_vars, 2, var) == 0))
    if (length(zero_var_cols) > 0) {
      print(paste("Removing zero-variance columns:", paste(zero_var_cols, collapse = ", ")))
      subset_vars <- subset_vars[, apply(subset_vars, 2, var) > 0, drop = FALSE]
    }
    
    if (ncol(subset_vars) == 0) {
      stop(paste("Error: No valid variables for PCA after removing zero-variance columns in species", species_label))
    }
    
    # Correlation heatmap
    cor_mat <- cor(subset_vars, use = "pairwise.complete.obs", method = "pearson")
    heatmap_file <- file.path(output_dir, paste0("Heatmap_", species_label, ".png"))
    pheatmap(cor_mat,
             color = colorRampPalette(c("blue", "white", "red"))(50),
             display_numbers = TRUE,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             filename = heatmap_file,
             width = 8,
             height = 6)
    
    # PCA
    print("Running PCA...")
    pca_res <- prcomp(subset_vars, scale. = TRUE)
    
    # Scree plot
    scree <- fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 60))
    scree_file <- file.path(output_dir, paste0("scree_", species_label, ".png"))
    ggsave(scree_file, plot = scree, width = 6, height = 5, units = "in", dpi = 400)
    
    # Correlation plot
    corr_plot <- fviz_pca_var(pca_res,
                              col.var = "contrib",
                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                              repel = TRUE,
                              title = paste("Variable Correlation -", species_label)) +
      coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.margin = margin(0, 0, 0, 0),
            aspect.ratio = 1,
            legend.position = "bottom")
    
    corr_file <- file.path(output_dir, paste0("PCA_", species_label, ".png"))
    
    if (crop == "ALL") {
      eig <- (pca_res$sdev)^2 / sum(pca_res$sdev^2) * 100
      pc1_var <- round(eig[1], 1)
      pc2_var <- round(eig[2], 1)
      
      biplot_vars <- fviz_pca_ind(pca_res, axes = c(1, 2),
                                  habillage = subset_data[[species_col]],
                                  addEllipses = TRUE,
                                  repel = TRUE,
                                  label = "none",
                                  palette = c("lettuce" = "darkorange", 
                                              "radish" = "purple", 
                                              "pea" = "green"),
                                  title = "PCA Individuals with Ellipses",
                                  legend.title = "Species") +
        labs(x = paste0("PC1 (", pc1_var, "%)"),
             y = paste0("PC2 (", pc2_var, "%)")) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom")
      
      combined_plot <- cowplot::plot_grid(biplot_vars, corr_plot,
                                          rel_widths = c(1.2, 0.8),
                                          align = "h",
                                          axis = "tb",
                                          ncol = 2)
      
      ggsave(file.path(output_dir, paste0("PCA_", species_label, ".png")),
             plot = combined_plot, width = 12, height = 5, units = "in", dpi = 400)
    } else {
      ggsave(corr_file, plot = corr_plot, width = 6, height = 5, units = "in", dpi = 400)
    }
    
    print(paste("Finished species:", species_label))
  }
  
  print("All species processed.")
  if (length(skipped_species) > 0) {
    print(paste("Skipped species:", paste(skipped_species, collapse = ", ")))
  } else {
    print("No species were skipped.")
  }
  
  return(top_vars_by_species)
}