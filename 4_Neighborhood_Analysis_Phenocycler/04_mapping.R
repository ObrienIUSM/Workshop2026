library(dplyr); library(ggplot2)

in_path  <- "C:/Users/ananamat/Desktop/Obrien webinar"
out_path <- "C:/Users/ananamat/Desktop/Obrien webinar"
scaled.subset <- readRDS(file.path(in_path, "scaled_umap_clusters_30um_NEclusters.rds"))

target_samples <- c("24-0165", "24-0168")
highlight_clusters <- c("2")
max_dim <- 3800
point_size <- 0.45
dpi <- 300
width_factor <- 1.5

plot_sample <- function(df, sample_id) {
  df <- df %>%
    filter(Biopsy == sample_id) %>%
    mutate(cluster_plot = ifelse(as.character(NE_cluster) %in% highlight_clusters, as.character(NE_cluster), "other"))
  if (nrow(df) == 0) return(NULL)
  
  xr <- range(df$X, na.rm = TRUE); yr <- range(df$Y, na.rm = TRUE)
  w <- max(diff(xr), 1); h <- max(diff(yr), 1)
  sf <- min(max_dim / w, max_dim / h, 1)
  
  p <- ggplot(df, aes(X, Y, color = cluster_plot)) + geom_point(size = point_size * sqrt(sf), alpha = 0.9) + scale_color_manual(values = c(setNames("red", highlight_clusters), other = "grey80"), breaks = highlight_clusters, name = "cluster") + guides(color = guide_legend(override.aes = list(size = 6))) + coord_fixed(xlim = xr, ylim = yr, expand = FALSE) + scale_y_reverse() + theme_void(base_size = 14) + theme(legend.position = "right", plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA)) + ggtitle(sample_id)
  
  outfile <- file.path(out_path, paste0(sample_id, "_mapping.png"))
  ggsave(outfile, plot = p, width = (w * sf / dpi) * width_factor, height = h * sf / dpi, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)}
lapply(target_samples, \(s) plot_sample(scaled.subset, s))


