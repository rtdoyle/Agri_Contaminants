div_amend_KW <- function(combs, times, metrics, crops, df){
  
  ### print traits
  print(combs)
  
  df.f <- df %>%
    filter(timepoint == times &
             species == crops) %>%
    droplevels(.)
  
  # dynamically assign the column to a new variable
  df.f$metric_value <- df.f[[metrics]]  
  
  # Kruskal-Wallis test
  kw <- kruskal.test(df.f$metric_value ~ df.f$Amendment)
  
  # Extract Kruskal-Wallis result
  kw_df <- data.frame(
    Test = "Kruskal-Wallis",
    Statistic = kw$statistic,
    P_value = kw$p.value,
    Comparison = "Overall",
    Species = paste(crops),
    Time = paste(times),
    Metric = paste(metrics),
    comb = paste(combs)
  )
  
  # Run Dunn's test with Benjamini-Hochberg correction
  dt <- dunnTest(metric_value ~ Amendment, data = df.f, method = "bh")
  
  # Extract Dunn's test results
  dt_df <- as.data.frame(dt$res)
  dt_df$Test <- "Dunn's test"
  dt_df$Statistic <- dt_df$Z
  dt_df$P_value <- dt_df$P.adj
  dt_df$Species = paste(crops)
  dt_df$Time = paste(times)
  dt_df$Metric = paste(metrics)
  dt_df$comb = paste(combs)
  
  # Combine both
  results_df <- rbind(
    kw_df[, c("Test", "Statistic", "P_value", "Comparison",
              "Species","Time","Metric","comb")],
    dt_df[, c("Test", "Statistic", "P_value", "Comparison",
                "Species","Time","Metric","comb")]
  )
  
  ### return dfs
  return(list(results_df)) 
}