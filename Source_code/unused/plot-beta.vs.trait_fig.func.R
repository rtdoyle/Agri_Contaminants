combine_crop_plots <- function(crop_code, crop_name, plots.out) {
  # Access plots
  p1 <- plots.out[[paste0("TP1_", crop_code)]][[1]]
  p2 <- plots.out[[paste0("TP2_", crop_code)]][[1]]
  
  # Build plots to extract axis limits
  p1_built <- ggplot_build(p1)
  p2_built <- ggplot_build(p2)
  
  x1_limits <- range(p1_built$layout$panel_scales_x[[1]]$range$range)
  x2_limits <- range(p2_built$layout$panel_scales_x[[1]]$range$range)
  
  y1_limits <- range(p1_built$layout$panel_scales_y[[1]]$range$range)
  y2_limits <- range(p2_built$layout$panel_scales_y[[1]]$range$range)
  
  # Add a buffer (e.g., 5% of the range)
  x_buffer <- 0.05 * diff(range(c(x1_limits, x2_limits)))
  y_buffer <- 0.05 * diff(range(c(y1_limits, y2_limits)))
  
  x_limits <- c(min(x1_limits, x2_limits) - x_buffer, max(x1_limits, x2_limits) + x_buffer)
  y_limits <- c(min(y1_limits, y2_limits) - y_buffer, max(y1_limits, y2_limits) + y_buffer)
  
  # Combine plots with consistent axis limits
  fig <- plot_grid(
    p1 +
      ggtitle(paste0(crop_name, ": 1 wk PP")) +
      theme(legend.position = "none", axis.title.x = element_blank()) +
      coord_cartesian(xlim = x_limits, ylim = y_limits),
    p2 +
      ggtitle(paste0(crop_name, ": Harvest")) +
      theme(legend.position = "none") +
      coord_cartesian(xlim = x_limits, ylim = y_limits),
    rel_widths = c(1, 1),
    rel_heights = c(1, 1.2),
    ncol = 1,
    nrow = 2,
    align = "v",
    labels = NULL
  )
  
  return(fig)
}