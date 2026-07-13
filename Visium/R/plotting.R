cluster_colors <- c(
  `0` = "#1289A7", `1` = "#A3CB38", `2` = "#F79F1F",
  `3` = "#D980FA", `4` = "#B53471", `5` = "#6F1E51",
  `6` = "#ED4C67", `7` = "#12CBC4", `8` = "gold"
)

region_colors <- c(Stroma = "#ffc200", Normal = "#0652DD", `Low-grade` = "#ff5252")
size_colors <- c(Small = "#0652DD", Large = "#ff5252")
moran_colors <- c(
  `Not significant` = "lightgrey",
  `High-High` = "#F5002D",
  `High-Low` = "#FF7A33",
  `Low-High` = "#A30080",
  `Low-Low` = "#000080"
)

save_plot <- function(plot, filename, width, height) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, units = "in")
  invisible(plot)
}

plot_spatial_samples <- function(
    object,
    group_by,
    samples,
    colors,
    ncol = 2,
    point_size = 1.6
) {
  plots <- Seurat::SpatialDimPlot(
    object,
    images = samples,
    group.by = group_by,
    cols = colors,
    combine = FALSE,
    pt.size.factor = point_size
  )
  patchwork::wrap_plots(plots, ncol = ncol, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

sample_level_boxplot <- function(data, x, y, sample_col = "Sample_ID", fill = x) {
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6) +
    ggplot2::geom_line(ggplot2::aes(group = .data[[sample_col]]), alpha = 0.3) +
    ggplot2::geom_point(size = 2) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL)
}

