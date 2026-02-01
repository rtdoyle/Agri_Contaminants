div_contam_KW <- function(combs, times, metrics, crops, amends, df){
  
  ### print traits
  print(combs)
  
  ## subset to times, species, amends
  df.f <- df %>%
    filter(timepoint == times &
             species == crops &
             Amendment == amends) %>%
    droplevels(.)
  
  # dynamically assign the column to a new variable
  df.f$metric_value <- df.f[[metrics]]  
  
  # Kruskal-Wallis test
  kw <- kruskal.test(df.f$metric_value ~ df.f$spikeFac)
  
  # Extract Kruskal-Wallis result
  kw_df <- data.frame(
    Test = "Kruskal-Wallis",
    Statistic = kw$statistic,
    P_value = kw$p.value,
    Comparison = "Overall",
    Species = paste(crops),
    Time = paste(times),
    Metric = paste(metrics),
    Amendment = paste(amends)
  )
  
  # Preserve the order for later use
  ordered_levels <- levels(df.f$spikeFac)
  
  # Convert to unordered factor for Dunn's test
  df.f$spikeFac <- factor(df.f$spikeFac, levels = ordered_levels, 
                          ordered = FALSE)
  
  # Run Dunn's test with Benjamini-Hochberg correction
  dt <- dunnTest(metric_value ~ spikeFac, data = df.f, method = "bh")
  
  # Extract Dunn's test results
  dt_df <- as.data.frame(dt$res)
  dt_df$Test <- "Dunn's test"
  dt_df$Statistic <- dt_df$Z
  dt_df$P_value <- dt_df$P.adj
  dt_df$Species = paste(crops)
  dt_df$Time = paste(times)
  dt_df$Amendment = paste(amends)
  dt_df$Metric = paste(metrics)
  
  # Combine both
  results_df <- rbind(
    kw_df[, c("Test", "Statistic", "P_value", "Comparison",
              "Species","Time","Metric","Amendment")],
    dt_df[, c("Test", "Statistic", "P_value", "Comparison",
                "Species","Time","Metric", "Amendment")]
  )
  
  ### return dfs
  return(list(results_df
  )) 
}