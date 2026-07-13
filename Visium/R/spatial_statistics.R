extract_feature <- function(object, feature, assay = NULL, layer = "data") {
  if (feature %in% colnames(object@meta.data)) {
    out <- object@meta.data[[feature]]
    names(out) <- rownames(object@meta.data)
    return(out)
  }
  if (is.null(assay)) {
    hits <- names(object@assays)[vapply(object@assays, function(x) feature %in% rownames(x), logical(1))]
    if (!length(hits)) stop("Feature not found: ", feature)
    assay <- hits[[1]]
  }
  data <- Seurat::FetchData(object, vars = feature, assay = assay, layer = layer)
  out <- data[[1]]
  names(out) <- rownames(data)
  out
}

spatial_metadata <- function(object, variables = character()) {
  pieces <- lapply(names(object@images), function(image) {
    coordinates <- as.data.frame(Seurat::GetTissueCoordinates(object, image = image)) |>
      tibble::rownames_to_column("barcode")
    metadata <- object@meta.data |>
      tibble::rownames_to_column("barcode") |>
      dplyr::filter(.data$Sample_ID == image)
    dplyr::left_join(metadata, coordinates, by = "barcode")
  })
  out <- dplyr::bind_rows(pieces)
  for (variable in variables) {
    values <- extract_feature(object, variable)
    out[[variable]] <- values[out$barcode]
  }
  out
}

bivariate_local_moran <- function(
    data,
    x,
    y,
    sample_col = "Sample_ID",
    coord_cols = c("x", "y"),
    k = 6,
    permutations = 999,
    seed = 1
) {
  assert_columns(data, c(sample_col, coord_cols, x, y))
  set.seed(seed)
  cluster_labels <- c(
    "Not significant", "High-High", "Low-Low", "Low-High",
    "High-Low", "Undefined", "Isolated"
  )

  data |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::group_modify(function(df, key) {
      complete <- stats::complete.cases(df[, c(coord_cols, x, y)])
      result <- tibble::tibble(
        barcode = df$barcode,
        local_bimoran = NA_real_,
        local_bimoran_p = NA_real_,
        local_bimoran_class = NA_character_
      )
      use <- df[complete, , drop = FALSE]
      if (nrow(use) <= k + 1) return(result)

      geometry <- sf::st_as_sf(use, coords = coord_cols, remove = FALSE)
      weights <- rgeoda::knn_weights(geometry, k = k)
      lisa <- rgeoda::local_bimoran(
        weights,
        use[, c(x, y), drop = FALSE],
        permutations = permutations,
        permutation_method = "complete"
      )
      idx <- match(use$barcode, result$barcode)
      result$local_bimoran[idx] <- rgeoda::lisa_values(lisa)
      result$local_bimoran_p[idx] <- rgeoda::lisa_pvalues(lisa)
      result$local_bimoran_class[idx] <- cluster_labels[rgeoda::lisa_clusters(lisa) + 1]
      result
    }) |>
    dplyr::ungroup()
}

.knn_lag <- function(values, coordinates, k = 6) {
  neighbors <- RANN::nn2(coordinates, coordinates, k = min(k + 1, nrow(coordinates)))$nn.idx[, -1, drop = FALSE]
  apply(neighbors, 1, function(index) mean(values[index], na.rm = TRUE))
}

spatial_lag_correlations <- function(
    data,
    x = "Senepy",
    y = "Stem_Signature",
    adjust_x_for = NULL,
    adjust_y_for = NULL,
    sample_col = "Sample_ID",
    group_col = "LargeSmall",
    coord_cols = c("array_row", "array_col"),
    k = 6,
    min_spots = 50
) {
  needed <- c(sample_col, group_col, coord_cols, x, y, adjust_x_for, adjust_y_for)
  assert_columns(data, needed)

  data |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::group_modify(function(df, key) {
      model_variables <- unique(c(x, y, adjust_x_for, adjust_y_for, coord_cols))
      keep <- stats::complete.cases(df[, model_variables, drop = FALSE])
      df <- df[keep, , drop = FALSE]
      group <- unique(as.character(df[[group_col]]))
      group <- if (length(group) == 1) group else NA_character_
      if (nrow(df) < min_spots) {
        return(tibble::tibble(n = nrow(df), group = group, rho = NA_real_, fisher_z = NA_real_))
      }

      x_values <- df[[x]]
      if (!is.null(adjust_x_for)) {
        model <- stats::lm(stats::reformulate(adjust_x_for, response = x), data = df)
        x_values <- stats::residuals(model)
      }
      y_values <- df[[y]]
      if (!is.null(adjust_y_for)) {
        model <- stats::lm(stats::reformulate(adjust_y_for, response = y), data = df)
        y_values <- stats::residuals(model)
      }
      lag_y <- .knn_lag(y_values, as.matrix(df[, coord_cols]), k = k)
      rho <- suppressWarnings(stats::cor(x_values, lag_y, method = "spearman"))
      rho <- pmin(pmax(rho, -0.999999), 0.999999)
      tibble::tibble(n = nrow(df), group = group, rho = rho, fisher_z = atanh(rho))
    }) |>
    dplyr::ungroup()
}

compare_before_after <- function(
    before,
    after,
    sample_col = "Sample_ID",
    group_col = "group",
    permutations = 999,
    seed = 1
) {
  joined <- before |>
    dplyr::select(
      dplyr::all_of(c(sample_col, group_col)),
      z_before = fisher_z,
      rho_before = rho
    ) |>
    dplyr::inner_join(
      after |>
        dplyr::select(dplyr::all_of(sample_col), z_after = fisher_z, rho_after = rho),
      by = sample_col
    ) |>
    dplyr::filter(is.finite(.data$z_before), is.finite(.data$z_after)) |>
    dplyr::mutate(delta_z = .data$z_after - .data$z_before)

  paired <- stats::wilcox.test(joined$z_after, joined$z_before, paired = TRUE, exact = FALSE)
  set.seed(seed)
  observed <- stats::median(joined$delta_z)
  null <- replicate(permutations, {
    signs <- sample(c(-1, 1), nrow(joined), replace = TRUE)
    stats::median(joined$delta_z * signs)
  })
  permutation_p <- (1 + sum(abs(null) >= abs(observed))) / (permutations + 1)

  stratified <- joined |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::group_modify(function(df, key) {
      if (nrow(df) < 3) {
        return(tibble::tibble(n = nrow(df), median_delta_z = NA_real_, wilcox_p = NA_real_, sign_flip_p = NA_real_))
      }
      group_observed <- stats::median(df$delta_z)
      group_null <- replicate(permutations, {
        signs <- sample(c(-1, 1), nrow(df), replace = TRUE)
        stats::median(df$delta_z * signs)
      })
      tibble::tibble(
        n = nrow(df),
        median_delta_z = group_observed,
        wilcox_p = stats::wilcox.test(df$z_after, df$z_before, paired = TRUE, exact = FALSE)$p.value,
        sign_flip_p = (1 + sum(abs(group_null) >= abs(group_observed))) / (permutations + 1)
      )
    }) |>
    dplyr::ungroup()

  list(table = joined, paired_wilcox = paired, sign_flip_p = permutation_p, stratified = stratified)
}

.zscore <- function(x) {
  value <- as.numeric(scale(x))
  value[!is.finite(value)] <- NA_real_
  value
}

prepare_cd8_neighborhood_data <- function(
    object,
    epithelial_clusters = c("0", "1", "2", "3", "4", "6", "8"),
    rna_assay = "Spatial",
    adt_assay = "ADT_dsb",
    gdf15 = "GDF15",
    cd3 = "ADT-CD3E",
    cd8 = "ADT-CD8A",
    cd45ra = "ADT-PTPRC",
    cd45ro = "ADT-PTPRC.2"
) {
  variables <- c(gdf15, cd3, cd8, cd45ra, cd45ro)
  df <- spatial_metadata(object)
  df[[gdf15]] <- extract_feature(object, gdf15, assay = rna_assay)[df$barcode]
  for (feature in c(cd3, cd8, cd45ra, cd45ro)) {
    df[[feature]] <- extract_feature(object, feature, assay = adt_assay)[df$barcode]
  }
  df <- df |>
    dplyr::filter(
      .data$Adenoma %in% c("Normal", "Low-grade"),
      as.character(.data$combined_cluster) %in% epithelial_clusters
    ) |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::mutate(
      z_gdf15 = .zscore(.data[[gdf15]]),
      z_cd3 = .zscore(.data[[cd3]]),
      z_cd8 = .zscore(.data[[cd8]]),
      z_cd45 = pmax(.zscore(.data[[cd45ra]]), .zscore(.data[[cd45ro]]), na.rm = TRUE),
      CD8_positive = .data$z_cd3 >= 0.5 & .data$z_cd8 >= 0.5 & .data$z_cd45 >= 0
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      is.finite(.data$z_cd3),
      is.finite(.data$z_cd8),
      is.finite(.data$z_cd45)
    ) |>
    dplyr::group_by(.data$Sample_ID, .data$Adenoma) |>
    dplyr::mutate(
      GDF15_group = ifelse(.data[[gdf15]] > stats::median(.data[[gdf15]], na.rm = TRUE), "High", "Low")
    ) |>
    dplyr::ungroup()
  add_tissue_components(df)
}

.spot_spacing <- function(df, coord_cols = c("array_col", "array_row")) {
  coordinates <- as.matrix(df[, coord_cols, drop = FALSE])
  stats::median(RANN::nn2(coordinates, coordinates, k = 2)$nn.dists[, 2], na.rm = TRUE)
}

add_tissue_components <- function(
    data,
    sample_col = "Sample_ID",
    coord_cols = c("array_col", "array_row"),
    radius_multiplier = 1.6,
    k = 6
) {
  data |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::group_modify(function(df, key) {
      coordinates <- as.matrix(df[, coord_cols, drop = FALSE])
      radius <- radius_multiplier * .spot_spacing(df, coord_cols)
      nearest <- RANN::nn2(coordinates, coordinates, k = min(k + 1, nrow(df)))
      edges <- lapply(seq_len(nrow(df)), function(i) {
        index <- nearest$nn.idx[i, -1]
        distance <- nearest$nn.dists[i, -1]
        index <- index[distance <= radius]
        if (!length(index)) return(NULL)
        cbind(i, index)
      })
      edges <- do.call(rbind, edges)
      if (is.null(edges)) {
        df$tissue_island <- seq_len(nrow(df))
      } else {
        graph <- igraph::graph_from_data_frame(
          data.frame(from = edges[, 1], to = edges[, 2]),
          directed = FALSE,
          vertices = data.frame(name = seq_len(nrow(df)))
        )
        membership <- igraph::components(graph)$membership
        df$tissue_island <- membership[as.character(seq_len(nrow(df)))]
      }
      df
    }) |>
    dplyr::ungroup()
}

.source_neighborhood_density <- function(
    df,
    source,
    radius_multiplier = 2,
    coord_cols = c("array_col", "array_row"),
    island_col = "tissue_island",
    target_col = "CD8_positive",
    min_sources = 5
) {
  source_df <- df[source, , drop = FALSE]
  if (nrow(source_df) < min_sources) return(NA_real_)
  radius <- radius_multiplier * .spot_spacing(df, coord_cols)
  density <- rep(NA_real_, nrow(source_df))

  for (i in seq_len(nrow(source_df))) {
    candidates <- df[df[[island_col]] == source_df[[island_col]][i], , drop = FALSE]
    delta <- sweep(as.matrix(candidates[, coord_cols, drop = FALSE]), 2,
                   as.numeric(source_df[i, coord_cols]), "-")
    nearby <- rowSums(delta^2) <= radius^2
    same_spot <- candidates$barcode == source_df$barcode[i]
    nearby <- nearby & !same_spot
    if (sum(nearby)) density[i] <- mean(candidates[[target_col]][nearby], na.rm = TRUE)
  }
  mean(density, na.rm = TRUE)
}

cd8_density_by_source <- function(
    data,
    comparison = c("region", "gdf15"),
    radius_multiplier = 2,
    min_sources = 5
) {
  comparison <- match.arg(comparison)
  if (comparison == "region") {
    return(
      data |>
        dplyr::group_by(.data$Sample_ID, .data$LargeSmall) |>
        dplyr::group_modify(function(df, key) {
          tibble::tibble(
            Normal = .source_neighborhood_density(
              df, df$Adenoma == "Normal" & df$GDF15_group == "High",
              radius_multiplier = radius_multiplier, min_sources = min_sources
            ),
            Low_grade = .source_neighborhood_density(
              df, df$Adenoma == "Low-grade" & df$GDF15_group == "High",
              radius_multiplier = radius_multiplier, min_sources = min_sources
            )
          )
        }) |>
        dplyr::ungroup()
    )
  }

  data |>
    dplyr::group_by(.data$Sample_ID, .data$LargeSmall, .data$Adenoma) |>
    dplyr::group_modify(function(df, key) {
      tibble::tibble(
        GDF15_low = .source_neighborhood_density(
          df, df$GDF15_group == "Low",
          radius_multiplier = radius_multiplier, min_sources = min_sources
        ),
        GDF15_high = .source_neighborhood_density(
          df, df$GDF15_group == "High",
          radius_multiplier = radius_multiplier, min_sources = min_sources
        )
      )
    }) |>
    dplyr::ungroup()
}

distance_to_surface <- function(
    data,
    surface_col = "surface",
    sample_col = "Sample_ID",
    island_col = "tissue_island",
    coord_cols = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    visium_spacing_um = 100
) {
  assert_columns(data, c(surface_col, sample_col, island_col, coord_cols, "barcode"))
  data |>
    dplyr::group_by(.data[[sample_col]], .data[[island_col]]) |>
    dplyr::group_modify(function(df, key) {
      coordinates <- as.matrix(df[, coord_cols, drop = FALSE])
      anchors <- as.logical(df[[surface_col]])
      if (!any(anchors, na.rm = TRUE)) {
        return(tibble::tibble(barcode = df$barcode, surface_distance_um = NA_real_))
      }
      spacing_px <- stats::median(FNN::get.knn(coordinates, k = 1)$nn.dist[, 1], na.rm = TRUE)
      distance_px <- FNN::get.knnx(coordinates[anchors, , drop = FALSE], coordinates, k = 1)$nn.dist[, 1]
      tibble::tibble(
        barcode = df$barcode,
        surface_distance_um = distance_px / spacing_px * visium_spacing_um
      )
    }) |>
    dplyr::ungroup()
}

summarize_distance_bins <- function(
    data,
    score,
    distance_col = "surface_distance_um",
    group_col = "Adenoma",
    surface_col = "surface",
    bin_width = 100,
    min_spots = 5
) {
  data |>
    dplyr::filter(!.data[[surface_col]], is.finite(.data[[distance_col]]), is.finite(.data[[score]])) |>
    dplyr::mutate(distance_bin_um = floor(.data[[distance_col]] / bin_width) * bin_width + bin_width / 2) |>
    dplyr::group_by(.data[[group_col]], .data$distance_bin_um) |>
    dplyr::summarise(
      mean_score = mean(.data[[score]], na.rm = TRUE),
      median_score = stats::median(.data[[score]], na.rm = TRUE),
      n_spots = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_spots >= min_spots)
}

separate_tissue_islands_semla <- function(object, roi_col = "tissue_roi") {
  if (is.null(object@tools$Staffli)) {
    stop("The Seurat object must contain a semla Staffli object in object@tools$Staffli")
  }
  object[[roi_col]] <- "Tissue"
  object <- semla::DisconnectRegions(object, column_name = roi_col, selected_groups = "Tissue")
  object$tissue_island <- object$Tissue_split
  object
}

compare_sample_summaries <- function(
    data,
    score,
    condition,
    sample_col = "Sample_ID",
    levels,
    paired,
    minimum_spots = 30,
    summary = stats::median
) {
  summarized <- data |>
    dplyr::filter(.data[[condition]] %in% levels, is.finite(.data[[score]])) |>
    dplyr::group_by(.data[[sample_col]], .data[[condition]]) |>
    dplyr::summarise(n = dplyr::n(), value = summary(.data[[score]], na.rm = TRUE), .groups = "drop") |>
    dplyr::filter(.data$n >= minimum_spots)

  if (paired) {
    wide <- tidyr::pivot_wider(summarized, names_from = dplyr::all_of(condition), values_from = .data$value) |>
      tidyr::drop_na(dplyr::all_of(levels))
    test <- stats::wilcox.test(wide[[levels[2]]], wide[[levels[1]]], paired = TRUE, exact = FALSE)
    effect <- stats::median(wide[[levels[2]]] - wide[[levels[1]]])
  } else {
    test <- stats::wilcox.test(
      stats::reformulate(condition, response = "value"),
      data = summarized,
      exact = FALSE
    )
    effect <- stats::median(summarized$value[summarized[[condition]] == levels[2]]) -
      stats::median(summarized$value[summarized[[condition]] == levels[1]])
    wide <- NULL
  }
  list(summary = summarized, paired_data = wide, test = test, median_difference = effect)
}

composition_tests <- function(
    data,
    category,
    condition,
    sample_col = "Sample_ID",
    paired = FALSE,
    permutations = 999,
    seed = 1
) {
  set.seed(seed)
  input <- data |>
    dplyr::transmute(
      .sample = as.character(.data[[sample_col]]),
      .condition = as.character(.data[[condition]]),
      .category = as.character(.data[[category]])
    ) |>
    dplyr::filter(!is.na(.data$.sample), !is.na(.data$.condition), !is.na(.data$.category))
  category_levels <- sort(unique(input$.category))
  condition_levels <- sort(unique(input$.condition))

  proportions <- input |>
    dplyr::count(.data$.sample, .data$.condition, .data$.category, name = "count") |>
    dplyr::group_by(.data$.sample, .data$.condition) |>
    tidyr::complete(.category = category_levels, fill = list(count = 0)) |>
    dplyr::mutate(proportion = .data$count / sum(.data$count)) |>
    dplyr::ungroup()
  wide <- proportions |>
    tidyr::pivot_wider(names_from = .data$.category, values_from = .data$proportion, values_fill = 0)
  matrix <- sqrt(as.matrix(wide[, setdiff(colnames(wide), c(".sample", ".condition")), drop = FALSE]))
  condition_vector <- wide$.condition
  permanova <- vegan::adonis2(matrix ~ condition_vector, permutations = permutations)
  dispersion <- vegan::betadisper(stats::dist(matrix), group = condition_vector)
  dispersion_test <- vegan::permutest(dispersion, permutations = permutations)

  posthoc <- proportions |>
    dplyr::group_by(.data$.category) |>
    dplyr::group_modify(function(df, key) {
      if (paired) {
        if (length(condition_levels) != 2) stop("Paired composition tests require two conditions")
        pair <- tidyr::pivot_wider(
          df,
          id_cols = .data$.sample,
          names_from = .data$.condition,
          values_from = .data$proportion
        ) |>
          tidyr::drop_na(dplyr::all_of(condition_levels))
        p <- stats::wilcox.test(
          pair[[condition_levels[2]]], pair[[condition_levels[1]]],
          paired = TRUE, exact = FALSE
        )$p.value
      } else {
        p <- stats::wilcox.test(df$proportion ~ df$.condition, exact = FALSE)$p.value
      }
      tibble::tibble(p_value = p)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      !!category := .data$.category,
      p_adjusted = stats::p.adjust(.data$p_value, method = "BH")
    ) |>
    dplyr::select(-.data$.category)

  list(
    proportions = proportions,
    permanova = permanova,
    dispersion = dispersion_test,
    posthoc = posthoc
  )
}
