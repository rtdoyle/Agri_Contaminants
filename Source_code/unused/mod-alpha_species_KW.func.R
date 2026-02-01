div_species_KW <- function(combs, times, metrics, df){
  
  ### print traits
  print(combs)
  
  df.f <- df %>%
    filter(timepoint == times) %>%
    droplevels(.)
  
  # dynamically assign the column to a new variable
  df.f$metric_value <- df.f[[metrics]]  
  
  # Kruskal-Wallis test
  kw <- kruskal.test(df.f$metric_value ~ df.f$species)
  
  # Extract Kruskal-Wallis result
  kw_df <- data.frame(
    Test = "Kruskal-Wallis",
    Statistic = kw$statistic,
    P_value = kw$p.value,
    Comparison = "Overall",
    Time = paste(times),
    Metric = paste(metrics)
  )
  
  # Run Dunn's test with Benjamini-Hochberg correction
  dt <- dunnTest(metric_value ~ species, data = df.f, method = "bh")
  
  # Extract Dunn's test results
  dt_df <- as.data.frame(dt$res)
  dt_df$Test <- "Dunn's test"
  dt_df$Statistic <- dt_df$Z
  dt_df$P_value <- dt_df$P.adj
  dt_df$Time = paste(times)
  dt_df$Metric = paste(metrics)
  
  # Combine both
  results_df <- rbind(
    kw_df[, c("Test", "Statistic", "P_value", "Comparison",
              "Time","Metric")],
    dt_df[, c("Test", "Statistic", "P_value", "Comparison",
              "Time","Metric")]
  )
  
  ### return dfs
  return(list(results_df)) 
}