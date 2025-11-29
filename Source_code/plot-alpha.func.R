plot_alpha_div <- function(metrics, df){
  
  ### print traits
  print(metrics)
  
  ## change timepoint to an unordered factor
  df$tpFac <- as.factor(df$timepoint)
  
  df_sum <- df %>%
    filter(treat.rn != "Garden Soil") %>%
    group_by(species, treat.rn, tp_treat, species_tp, species_treat, tpFac) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE),
              min = min(get(metrics), na.rm = TRUE),
              max = max(get(metrics), na.rm = TRUE))
  
  df_sum.GS <- df %>%
    filter(treat.rn == "Garden Soil") %>%
    droplevels(.) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE))
  
  df_sum.C <- df %>%
    filter(treat.rn == "control-0") %>%
    droplevels(.) %>%
    group_by(species, treat.rn, tp_treat, species_tp, species_treat, tpFac) %>%
    summarize(mean = mean(get(metrics), na.rm = TRUE))
  
    ### specify species names
  species_names <- c(
    L = 'Lettuce',
    R = 'Radish',
    P = 'Pea'
  )
  
  ## custom labels
  time_labs <- c(
    TP1 = "1 wk PP",
    TP2 = "Harvest"
  )
  
  ### specify species names
  species_tp.names <- c(
    L_TP1 = 'Let: 1wPP',
    L_TP2 = 'Let: Har',
    R_TP1 = 'Rad: 1wPP',
    R_TP2 = 'Rad: Har',
    P_TP1 = 'Pea: 1wPP',
    P_TP2 = 'Pea: Har'
  )
  
  ### plots
  p <- ggplot(df_sum,
              aes(x = treat.rn, 
                  y = mean)) +
    geom_col(
      aes(fill = treat.rn)
    ) +
    geom_beeswarm(data = df %>%
                    filter(treat.rn != "Garden Soil"),
                  aes(x = treat.rn, y = get(metrics))) +
    geom_hline(aes(yintercept = df_sum.GS$mean),
               linetype = 1) +
    geom_hline(data = df_sum.C,
      aes(yintercept = mean),
               linetype = 2) +
    scale_fill_manual(name = "Amendment- \n Spiking Level",
                      values = c("gray",
                                 "lightblue",
                                 "cornflowerblue",
                                 "blue",
                                 "#EABD8C",
                                 "#FFAD00",
                                 "#B06500"),
                      guide = guide_legend(nrow = 4, ncol = 2)) +
    facet_wrap(.~species_tp, scales = "fixed", ncol = 6,
               labeller = 
                 labeller(species_tp = species_tp.names)) +
    labs(y = paste0(metrics), 
         x = "Amendment type + Contaminant level") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),  # removes tick marks
      axis.title = element_text(size = 16),
      strip.text = element_text(size = 14, face = "bold"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ### return
  return(p)
}