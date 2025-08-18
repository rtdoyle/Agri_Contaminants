permanova.func <- function(times, df){
  
  ## times
  time <- as.character(times)
  print(time)
  
  ## exclude non-experiment samples
  metadata.f <- df %>%
    filter(timepoint %in% time &
             species != "no_plant") %>%
    droplevels(.)
  
  ## samples to include
  sample_list <- unique(metadata.f$NAME)
  
  ## subset OTU table to these samples
  otu_subset <- otu_table[rownames(otu_table) %in% sample_list, ]
  
  # Bray-Curtis dissimilarity
  bray_dist <- vegdist(otu_subset, method = "bray")
  
  
  # Example with multiple variables
  adonis2_result <- adonis2(bray_dist ~ treat.rn + Species, 
                            data = metadata.f, permutations = 999,
                            by = "terms")
  
  as.data.frame(adonis2_result)
  adonis2_result$timepoint <- paste0(times)
  
  ## check dispersion
  dispersion <- betadisper(bray_dist, metadata.f$treat.rn)
  print(anova(dispersion))
  boxplot(dispersion)
  plot(dispersion)
  permutest(dispersion, permutations = 999)
  ## still significant dispersion between timepoints
  
  return(list(adonis2_result))
    
}