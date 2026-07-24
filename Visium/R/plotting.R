cluster_colors <- c(
  `0` = "#1289A7", `1` = "#A3CB38", `2` = "#F79F1F",
  `3` = "#D980FA", `4` = "#B53471", `5` = "#6F1E51",
  `6` = "#ED4C67", `7` = "#12CBC4", `8` = "gold"
)

region_colors <- c(Stroma = "#ffc200", Normal = "#0652DD", `Low-grade` = "#ff5252")
size_colors <- c(Small = "#0652DD", Large = "#ff5252")
moran_colors <- c(
  `Not sig` = "lightgrey",
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

# Released implementation of Scripts/Utils/DotPlot_by_class_horizontal.R.
DotPlot_by_class_horizontal <- function(
    object,
    features,
    feature_class,
    cols = c("grey90", "red"),
    col.low = "#441153",
    col.mid = "grey90",
    col.high = "#3CB67B",
    col.min = -2,
    col.max = 2,
    scale = TRUE
) {
  plot_data <- Seurat::DotPlot(
    object,
    features = features,
    cols = cols,
    col.min = col.min,
    col.max = col.max,
    scale = scale
  )$data |>
    dplyr::left_join(feature_class, by = "features.plot") |>
    dplyr::mutate(
      features.plot = factor(.data$features.plot, levels = rev(features)),
      class = factor(.data$class, levels = unique(feature_class$class))
    )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$id, y = .data$features.plot)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$pct.exp, color = .data$avg.exp.scaled)) +
    ggplot2::facet_grid(class ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_color_gradient2(low = col.low, mid = col.mid, high = col.high) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
      strip.text.y = ggplot2::element_text(face = "bold", angle = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Clean function form of the active analysis in Scripts/Utils/3_way_Spatial_plot.R.
.cap_and_normalize <- function(x, quantile_cutoff = 0.99) {
  cap <- stats::quantile(x, quantile_cutoff, na.rm = TRUE, names = FALSE)
  finite <- is.finite(x)
  if (!any(finite) || !is.finite(cap)) return(rep(0, length(x)))
  x[finite & x > cap] <- cap
  value_range <- range(x, na.rm = TRUE)
  if (!all(is.finite(value_range)) || diff(value_range) == 0) {
    return(rep(0, length(x)))
  }
  normalized <- (x - value_range[[1]]) / diff(value_range)
  normalized[!finite] <- 0
  normalized
}

three_way_spatial_plot <- function(
    object,
    features = c("Dysplastic_Epi", "CD8_T_cells", "Normal_Epi"),
    images = NULL,
    colors = c("#E41A1C", "#5ce828", "#377EB8"),
    legend_titles = c("Dysplastic Epi", "T cells", "Normal Epi"),
    quantile_cutoff = 0.99,
    point_size = 1.5
) {
  if (length(features) != 3L || length(colors) != 3L || length(legend_titles) != 3L) {
    stop("features, colors, and legend_titles must each have length three")
  }

  abundance <- Seurat::FetchData(object, vars = features)
  normalized <- lapply(abundance, .cap_and_normalize, quantile_cutoff = quantile_cutoff)
  target_rgb <- grDevices::col2rgb(colors) / 255

  blended_channels <- vapply(seq_len(nrow(target_rgb)), function(channel) {
    contribution <- vapply(seq_along(normalized), function(index) {
      normalized[[index]] * (1 - target_rgb[channel, index])
    }, numeric(nrow(abundance)))
    pmax(0, 1 - rowSums(contribution))
  }, numeric(nrow(abundance)))

  blend_colors <- grDevices::rgb(
    red = blended_channels[, 1],
    green = blended_channels[, 2],
    blue = blended_channels[, 3],
    alpha = 1,
    maxColorValue = 1
  )
  names(blend_colors) <- rownames(abundance)
  object$Coloc_RGB <- blend_colors[colnames(object)]

  spatial_plot <- Seurat::SpatialDimPlot(
    object,
    group.by = "Coloc_RGB",
    pt.size.factor = point_size,
    image.alpha = 1,
    stroke = 0,
    alpha = NA,
    images = images
  ) +
    ggplot2::scale_fill_identity() +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12, color = "black")
    )

  make_gradient_legend <- function(high_color, title) {
    dummy <- data.frame(abundance = c(0, 1))
    legend_plot <- ggplot2::ggplot(
      dummy,
      ggplot2::aes(x = .data$abundance, y = .data$abundance, color = .data$abundance)
    ) +
      ggplot2::geom_point(alpha = 0) +
      ggplot2::scale_color_gradient(
        low = "white",
        high = high_color,
        breaks = c(0, 1),
        labels = c("Min", "Max")
      ) +
      ggplot2::labs(color = title) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold", size = 11, color = "black"),
        legend.text = ggplot2::element_text(size = 9, color = "black"),
        legend.background = ggplot2::element_rect(fill = "transparent", color = NA),
        legend.key.height = grid::unit(0.5, "cm")
      )
    cowplot::get_legend(legend_plot)
  }

  legends <- Map(make_gradient_legend, colors, legend_titles)
  legend_column <- cowplot::plot_grid(plotlist = legends, ncol = 1, align = "v")
  cowplot::plot_grid(spatial_plot, legend_column, rel_widths = c(1, 0.15), align = "h")
}

# Reduction-space counterpart to three_way_spatial_plot(). This reproduces the
# three-channel Cell2location blend used in Figure 3E without requiring SCP.
three_way_reduction_plot <- function(
    object,
    features = c("Dysplastic_Epi", "CD8_T_cells", "Normal_Epi"),
    reduction = "wnn.umap",
    colors = c("#E41A1C", "#5ce828", "#377EB8"),
    legend_titles = c("Dysplastic Epi", "T cells", "Normal Epi"),
    quantile_cutoff = 0.99,
    point_size = 0.5
) {
  if (length(features) != 3L || length(colors) != 3L || length(legend_titles) != 3L) {
    stop("features, colors, and legend_titles must each have length three")
  }
  if (!(reduction %in% names(object@reductions))) {
    stop("Reduction is absent from the object: ", reduction)
  }

  coordinates <- Seurat::Embeddings(object, reduction = reduction)[, 1:2, drop = FALSE]
  abundance <- Seurat::FetchData(object, vars = features)[rownames(coordinates), , drop = FALSE]
  normalized <- lapply(abundance, .cap_and_normalize, quantile_cutoff = quantile_cutoff)
  target_rgb <- grDevices::col2rgb(colors) / 255

  blended_channels <- vapply(seq_len(nrow(target_rgb)), function(channel) {
    contribution <- vapply(seq_along(normalized), function(index) {
      normalized[[index]] * (1 - target_rgb[channel, index])
    }, numeric(nrow(abundance)))
    pmax(0, 1 - rowSums(contribution))
  }, numeric(nrow(abundance)))

  plot_data <- data.frame(
    dim_1 = coordinates[, 1],
    dim_2 = coordinates[, 2],
    blend_color = grDevices::rgb(
      red = blended_channels[, 1],
      green = blended_channels[, 2],
      blue = blended_channels[, 3],
      alpha = 1,
      maxColorValue = 1
    )
  )

  reduction_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$dim_1, y = .data$dim_2, color = .data$blend_color)
  ) +
    ggplot2::geom_point(size = point_size, alpha = 1) +
    ggplot2::scale_color_identity() +
    ggplot2::coord_equal() +
    ggplot2::theme_void()

  make_gradient_legend <- function(high_color, title) {
    dummy <- data.frame(abundance = c(0, 1))
    legend_plot <- ggplot2::ggplot(
      dummy,
      ggplot2::aes(x = .data$abundance, y = .data$abundance, color = .data$abundance)
    ) +
      ggplot2::geom_point(alpha = 0) +
      ggplot2::scale_color_gradient(
        low = "white",
        high = high_color,
        breaks = c(0, 1),
        labels = c("Min", "Max")
      ) +
      ggplot2::labs(color = title) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold", size = 11),
        legend.text = ggplot2::element_text(size = 9),
        legend.key.height = grid::unit(0.5, "cm")
      )
    cowplot::get_legend(legend_plot)
  }

  legends <- Map(make_gradient_legend, colors, legend_titles)
  legend_column <- cowplot::plot_grid(plotlist = legends, ncol = 1, align = "v")
  cowplot::plot_grid(reduction_plot, legend_column, rel_widths = c(1, 0.18), align = "h")
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
