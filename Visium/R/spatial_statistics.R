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
  df <- spatial_metadata(object)
  df[[gdf15]] <- extract_feature(object, gdf15, assay = rna_assay)[df$barcode]

  adt_object <- SeuratObject::JoinLayers(object, assay = adt_assay)
  for (feature in c(cd3, cd8, cd45ra, cd45ro)) {
    df[[feature]] <- extract_feature(adt_object, feature, assay = adt_assay)[df$barcode]
  }

  df <- df |>
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(c(gdf15, cd3, cd8, cd45ra, cd45ro)),
        is.finite
      )
    ) |>
    dplyr::group_by(.data$Sample_ID) |>
    dplyr::mutate(
      z_GDF15 = .zscore(.data[[gdf15]]),
      z_CD3E = .zscore(.data[[cd3]]),
      z_CD8A = .zscore(.data[[cd8]]),
      z_PTPRC = pmax(.zscore(.data[[cd45ra]]), .zscore(.data[[cd45ro]]), na.rm = TRUE),
      CD8T = .data$z_CD3E >= 0.5 & .data$z_CD8A >= 0.5 & .data$z_PTPRC >= 0
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$Sample_ID, .data$Adenoma) |>
    dplyr::mutate(
      GDF15_hi = .data$z_GDF15 > stats::median(.data$z_GDF15, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$Adenoma %in% c("Normal", "Low-grade"),
      as.character(.data$combined_cluster) %in% epithelial_clusters,
      is.finite(.data$z_CD3E),
      is.finite(.data$z_CD8A),
      is.finite(.data$z_PTPRC)
    )

  df <- add_tissue_components(df, radius_multiplier = 1.6, k = 6)
  df |>
    dplyr::group_by(.data$Sample_ID, .data$Adenoma) |>
    dplyr::mutate(
      median_z_GDF15 = stats::median(.data$z_GDF15, na.rm = TRUE),
      GDF15_group = dplyr::case_when(
        .data$z_GDF15 > .data$median_z_GDF15 ~ "High",
        TRUE ~ "Low"
      )
    ) |>
    dplyr::ungroup()
}

.spot_spacing <- function(df, coord_cols = c("array_col", "array_row")) {
  if (nrow(df) < 2) return(NA_real_)
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
      if (nrow(df) < 2) {
        df$tissue_island <- seq_len(nrow(df))
        return(df)
      }
      coordinates <- as.matrix(df[, coord_cols, drop = FALSE])
      radius <- radius_multiplier * .spot_spacing(df, coord_cols)
      nearest <- RANN::nn2(coordinates, coordinates, k = min(k, nrow(df)))
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
    target_col = "CD8T",
    min_sources = 5,
    min_component_size = 5
) {
  source <- !is.na(source) & as.logical(source)
  source_df <- df[source, , drop = FALSE]
  empty_result <- list(
    value = NA_real_,
    n_src_total = nrow(source_df),
    n_src_used = 0L,
    n_tgt = sum(!is.na(df[[target_col]]) & as.logical(df[[target_col]]))
  )
  if (nrow(source_df) < min_sources) return(empty_result)

  radius <- radius_multiplier * .spot_spacing(df, coord_cols)
  if (!is.finite(radius)) return(empty_result)
  density <- rep(NA_real_, nrow(source_df))

  for (component in unique(source_df[[island_col]])) {
    source_component <- source_df[source_df[[island_col]] == component, , drop = FALSE]
    candidates <- df[df[[island_col]] == component, , drop = FALSE]
    if (nrow(candidates) < min_component_size) next

    source_coordinates <- as.matrix(source_component[, coord_cols, drop = FALSE])
    candidate_coordinates <- as.matrix(candidates[, coord_cols, drop = FALSE])
    delta_x <- outer(source_coordinates[, 1], candidate_coordinates[, 1], "-")
    delta_y <- outer(source_coordinates[, 2], candidate_coordinates[, 2], "-")
    distances <- sqrt(delta_x^2 + delta_y^2)

    nearby_count <- rowSums(distances <= radius) - 1L
    cd8_candidates <- !is.na(candidates[[target_col]]) & as.logical(candidates[[target_col]])
    if (any(cd8_candidates)) {
      cd8_coordinates <- candidate_coordinates[cd8_candidates, , drop = FALSE]
      cd8_delta_x <- outer(source_coordinates[, 1], cd8_coordinates[, 1], "-")
      cd8_delta_y <- outer(source_coordinates[, 2], cd8_coordinates[, 2], "-")
      cd8_distances <- sqrt(cd8_delta_x^2 + cd8_delta_y^2)
      cd8_count <- rowSums(cd8_distances <= radius)
      cd8_count <- cd8_count - as.integer(source_component[[target_col]] %in% TRUE)
    } else {
      cd8_count <- rep(0L, nrow(source_component))
    }

    component_density <- ifelse(nearby_count > 0, cd8_count / nearby_count, NA_real_)
    density[source_df[[island_col]] == component] <- component_density
  }

  if (!any(is.finite(density))) return(empty_result)
  list(
    value = mean(density, na.rm = TRUE),
    n_src_total = nrow(source_df),
    n_src_used = sum(is.finite(density)),
    n_tgt = empty_result$n_tgt
  )
}

cd8_density_by_source <- function(
    data,
    comparison = c("region", "gdf15"),
    radius_multiplier = 2,
    min_sources = 5,
    exclude_sample_pattern = "Pilot"
) {
  comparison <- match.arg(comparison)
  if (comparison == "region") {
    return(
      data |>
        dplyr::group_by(.data$Sample_ID, .data$LargeSmall) |>
        dplyr::group_modify(function(df, key) {
          low_grade <- .source_neighborhood_density(
            df,
            df$Adenoma == "Low-grade" & df$GDF15_hi,
            radius_multiplier = radius_multiplier,
            min_sources = min_sources,
            min_component_size = 5
          )
          normal <- .source_neighborhood_density(
            df,
            df$Adenoma == "Normal" & df$GDF15_hi,
            radius_multiplier = radius_multiplier,
            min_sources = min_sources,
            min_component_size = 5
          )
          tibble::tibble(
            dens_LGD = low_grade$value,
            dens_Normal = normal$value,
            delta_dens = low_grade$value - normal$value,
            n_src_LGD = low_grade$n_src_total,
            n_src_LGD_used = low_grade$n_src_used,
            n_src_Normal = normal$n_src_total,
            n_src_Normal_used = normal$n_src_used,
            n_tgt_epi = max(low_grade$n_tgt, normal$n_tgt)
          )
        }) |>
        dplyr::ungroup()
    )
  }

  data |>
    dplyr::filter(!grepl(exclude_sample_pattern, .data$Sample_ID)) |>
    dplyr::group_by(.data$Sample_ID, .data$LargeSmall, .data$Adenoma) |>
    dplyr::group_modify(function(df, key) {
      high <- .source_neighborhood_density(
        df,
        df$GDF15_group == "High",
        radius_multiplier = radius_multiplier,
        min_sources = min_sources,
        min_component_size = 3
      )
      low <- .source_neighborhood_density(
        df,
        df$GDF15_group == "Low",
        radius_multiplier = radius_multiplier,
        min_sources = min_sources,
        min_component_size = 3
      )
      tibble::tibble(
        dens_High = high$value,
        dens_Low = low$value,
        delta_HL = high$value - low$value
      )
    }) |>
    dplyr::ungroup()
}

.first_existing_file <- function(paths) {
  matches <- paths[file.exists(paths)]
  if (length(matches)) matches[[1]] else NA_character_
}

.read_spaceranger_coordinates <- function(path, sample_id) {
  has_header <- grepl("^barcode", readLines(path, n = 1))
  if (has_header) {
    coordinates <- readr::read_csv(path, show_col_types = FALSE)
  } else {
    coordinates <- readr::read_csv(
      path,
      col_names = c(
        "barcode", "in_tissue", "array_row", "array_col",
        "pxl_row_in_fullres", "pxl_col_in_fullres"
      ),
      show_col_types = FALSE
    )
  }
  dplyr::mutate(coordinates, Sample_ID = sample_id)
}

build_staffli_from_spaceranger <- function(
    object,
    manifest,
    image_height = 300,
    load_images = TRUE,
    verbose = TRUE
) {
  assert_columns(object@meta.data, "Sample_ID")
  assert_columns(manifest, c("sample_id", "space_ranger_outs"))

  sample_order <- unique(as.character(object$Sample_ID))
  inputs <- manifest |>
    dplyr::filter(.data$sample_id %in% sample_order) |>
    dplyr::distinct(.data$sample_id, .keep_all = TRUE) |>
    dplyr::arrange(match(.data$sample_id, sample_order)) |>
    dplyr::mutate(
      spatial_dir = file.path(.data$space_ranger_outs, "spatial"),
      coordinate_file = purrr::map_chr(
        .data$spatial_dir,
        ~ .first_existing_file(file.path(.x, c("tissue_positions.csv", "tissue_positions_list.csv")))
      ),
      image_file = purrr::map_chr(
        .data$spatial_dir,
        ~ .first_existing_file(file.path(
          .x,
          c(
            "tissue_lowres_image.png", "tissue_hires_image.png",
            "tissue_lowres_image.jpg", "tissue_hires_image.jpg"
          )
        ))
      ),
      scalefactor_file = file.path(.data$spatial_dir, "scalefactors_json.json")
    )

  missing_samples <- setdiff(sample_order, inputs$sample_id)
  if (length(missing_samples)) {
    stop("The manifest is missing samples in the Seurat object: ", paste(missing_samples, collapse = ", "))
  }
  required_files <- c(inputs$coordinate_file, inputs$image_file, inputs$scalefactor_file)
  if (anyNA(required_files) || any(!file.exists(required_files))) {
    stop("Staffli construction requires Space Ranger coordinates, an H&E image, and scalefactors for every sample")
  }

  raw_coordinates <- purrr::map2_dfr(
    inputs$coordinate_file,
    inputs$sample_id,
    .read_spaceranger_coordinates
  )
  object_key <- object@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::transmute(
      barcode = .data$barcode,
      Sample_ID = as.character(.data$Sample_ID),
      raw_barcode = stringr::str_extract(.data$barcode, "[ACGT]{16}-[0-9]+")
    )
  if (anyNA(object_key$raw_barcode)) {
    stop("Could not extract a raw 10x barcode from every Seurat spot name")
  }

  staffli_coordinates <- object_key |>
    dplyr::left_join(
      dplyr::rename(raw_coordinates, raw_barcode = barcode),
      by = c("Sample_ID", "raw_barcode")
    )
  if (anyNA(staffli_coordinates$pxl_row_in_fullres) ||
      anyNA(staffli_coordinates$pxl_col_in_fullres)) {
    unmatched <- staffli_coordinates |>
      dplyr::filter(is.na(.data$pxl_row_in_fullres) | is.na(.data$pxl_col_in_fullres)) |>
      dplyr::count(.data$Sample_ID)
    stop(
      "Some Seurat spots could not be matched to Space Ranger coordinates:\n",
      paste0(unmatched$Sample_ID, ": ", unmatched$n, collapse = "\n")
    )
  }

  staffli_coordinates <- staffli_coordinates |>
    dplyr::transmute(
      barcode = .data$barcode,
      x = .data$array_col,
      y = .data$array_row,
      pxl_col_in_fullres = .data$pxl_col_in_fullres,
      pxl_row_in_fullres = .data$pxl_row_in_fullres,
      sampleID = as.integer(match(.data$Sample_ID, inputs$sample_id)),
      Sample_ID = .data$Sample_ID,
      raw_barcode = .data$raw_barcode,
      in_tissue = .data$in_tissue,
      array_row = .data$array_row,
      array_col = .data$array_col
    )
  staffli_coordinates <- staffli_coordinates[
    match(colnames(object), staffli_coordinates$barcode),
    ,
    drop = FALSE
  ]
  stopifnot(identical(staffli_coordinates$barcode, colnames(object)))

  image_info <- semla::LoadImageInfo(inputs$image_file)
  scalefactors <- semla::LoadScaleFactors(inputs$scalefactor_file)
  image_info <- semla::UpdateImageInfo(image_info, scalefactors)
  if (nrow(image_info) != nrow(inputs) || nrow(scalefactors) != nrow(inputs)) {
    stop("semla returned an unexpected number of image-information or scalefactor rows")
  }
  image_info$sampleID <- as.character(seq_len(nrow(inputs)))
  scalefactors$sampleID <- as.character(seq_len(nrow(inputs)))

  object@tools$Staffli <- semla::CreateStaffliObject(
    imgs = inputs$image_file,
    meta_data = tibble::as_tibble(staffli_coordinates),
    image_info = image_info,
    scalefactors = scalefactors
  )
  stopifnot(identical(object@tools$Staffli@meta_data$barcode, colnames(object)))

  if (isTRUE(load_images)) {
    object <- semla::LoadImages(
      object,
      image_height = image_height,
      verbose = verbose
    )
  }
  object
}

distance_to_surface <- function(
    data,
    surface_col = "surface_label",
    surface_value = "All",
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
      anchors <- !is.na(df[[surface_col]]) & as.character(df[[surface_col]]) == surface_value
      if (nrow(df) < 2 || !any(anchors)) {
        return(tibble::tibble(
          barcode = df$barcode,
          tract_surface_depth_px = NA_real_,
          tract_surface_depth_um = NA_real_,
          spot_spacing_px = NA_real_,
          n_surface_spots = sum(anchors)
        ))
      }
      spacing_px <- stats::median(FNN::get.knn(coordinates, k = 1)$nn.dist[, 1], na.rm = TRUE)
      distance_px <- FNN::get.knnx(coordinates[anchors, , drop = FALSE], coordinates, k = 1)$nn.dist[, 1]
      tibble::tibble(
        barcode = df$barcode,
        tract_surface_depth_px = distance_px,
        tract_surface_depth_um = distance_px / spacing_px * visium_spacing_um,
        spot_spacing_px = spacing_px,
        n_surface_spots = sum(anchors)
      )
    }) |>
    dplyr::ungroup()
}

summarize_distance_bins <- function(
    data,
    score,
    distance_col = "tract_surface_depth_um",
    group_col = "Adenoma",
    surface_col = "surface_label",
    surface_value = "All",
    bin_width = 100,
    minimum_distance = 100,
    min_spots = 5
) {
  assert_columns(data, c(score, distance_col, group_col, surface_col))
  data |>
    dplyr::filter(
      as.character(.data[[surface_col]]) != surface_value,
      is.finite(.data[[distance_col]]),
      .data[[distance_col]] > minimum_distance
    ) |>
    dplyr::mutate(
      distance_bin_um = ceiling(.data[[distance_col]] / bin_width) * bin_width - bin_width / 2
    ) |>
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
  if (!"Tissue_split" %in% colnames(object@meta.data)) {
    stop("semla::DisconnectRegions() did not create the expected Tissue_split metadata")
  }
  object$tissue_island <- object$Tissue_split
  object
}

compare_score_clusters_paired <- function(
    md,
    score_col = "Senepy",
    sample_col = "Sample_ID",
    condition_col,
    level_a,
    level_b,
    condition_level = c("spot", "sample"),
    agg = c("median", "mean"),
    min_spots = 30,
    drop_na_score = TRUE
) {
  condition_level <- match.arg(condition_level)
  agg <- match.arg(agg)
  df <- md
  if (drop_na_score) {
    df <- dplyr::filter(df, !is.na(.data[[score_col]]))
  }

  if (condition_level == "spot") {
    df_sum <- df |>
      dplyr::filter(.data[[condition_col]] %in% c(level_a, level_b)) |>
      dplyr::group_by(.data[[sample_col]], .data[[condition_col]]) |>
      dplyr::summarise(
        n = dplyr::n(),
        score = if (agg == "median") {
          stats::median(.data[[score_col]], na.rm = TRUE)
        } else {
          mean(.data[[score_col]], na.rm = TRUE)
        },
        .groups = "drop"
      ) |>
      dplyr::filter(.data$n >= min_spots)
    wide <- df_sum |>
      dplyr::select(
        dplyr::all_of(sample_col),
        dplyr::all_of(condition_col),
        dplyr::all_of("score")
      ) |>
      tidyr::pivot_wider(
        names_from = dplyr::all_of(condition_col),
        values_from = "score"
      ) |>
      dplyr::filter(!is.na(.data[[level_a]]), !is.na(.data[[level_b]]))
    test <- stats::wilcox.test(
      wide[[level_a]],
      wide[[level_b]],
      paired = TRUE,
      exact = FALSE
    )
    return(list(
      mode = "paired_within_sample",
      wide = wide,
      n_samples = nrow(wide),
      test = test,
      df_sum = df_sum,
      ave_median = c(
        stats::median(wide[[level_a]]),
        stats::median(wide[[level_b]])
      ),
      diff_median = stats::median(wide[[level_a]] - wide[[level_b]], na.rm = TRUE)
    ))
  }

  sample_labels <- df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      condition = unique(stats::na.omit(as.character(.data[[condition_col]]))),
      .groups = "drop"
    )
  invalid_samples <- sample_labels |>
    dplyr::filter(lengths(.data$condition) != 1)
  if (nrow(invalid_samples) > 0) {
    stop(
      "Some samples have multiple values in ", condition_col,
      ". Example sample(s): ",
      paste(utils::head(invalid_samples[[sample_col]], 5), collapse = ", ")
    )
  }
  sample_labels$condition <- vapply(sample_labels$condition, `[`, character(1), 1)

  df_sum <- df |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      n = dplyr::n(),
      score = if (agg == "median") {
        stats::median(.data[[score_col]], na.rm = TRUE)
      } else {
        mean(.data[[score_col]], na.rm = TRUE)
      },
      .groups = "drop"
    ) |>
    dplyr::left_join(sample_labels, by = stats::setNames(sample_col, sample_col)) |>
    dplyr::filter(.data$condition %in% c(level_a, level_b), .data$n >= min_spots)
  group_a <- dplyr::pull(dplyr::filter(df_sum, .data$condition == level_a), "score")
  group_b <- dplyr::pull(dplyr::filter(df_sum, .data$condition == level_b), "score")
  test <- stats::wilcox.test(group_a, group_b, paired = FALSE, exact = FALSE)
  list(
    mode = "unpaired_between_samples",
    per_sample = df_sum,
    n_a = sum(df_sum$condition == level_a),
    n_b = sum(df_sum$condition == level_b),
    test = test,
    ave_median = c(
      stats::median(group_a, na.rm = TRUE),
      stats::median(group_b, na.rm = TRUE)
    ),
    diff_median = stats::median(group_a, na.rm = TRUE) -
      stats::median(group_b, na.rm = TRUE)
  )
}

run_comp_tests <- function(
    meta,
    sample_col,
    condition_col,
    category_col,
    condition_level = c("sample", "spot"),
    condition_keep = NULL,
    category_levels = NULL,
    n_perm = 999,
    seed = 1,
    control_level = NULL,
    case_level = NULL,
    paired_only = TRUE,
    drop_empty_strata = TRUE,
    na_category_label = NULL
) {
  condition_level <- match.arg(condition_level)
  set.seed(seed)

  df <- meta |>
    dplyr::mutate(
      .sample = as.character(.data[[sample_col]]),
      .cond = as.character(.data[[condition_col]]),
      .cat = as.character(.data[[category_col]])
    ) |>
    dplyr::filter(!is.na(.data$.sample), !is.na(.data$.cond))

  if (!is.null(condition_keep)) {
    df <- dplyr::filter(df, .data$.cond %in% condition_keep)
  }
  if (!is.null(na_category_label)) {
    missing_category <- is.na(df$.cat) | df$.cat == "NA" | df$.cat == ""
    df$.cat[missing_category] <- na_category_label
  }
  df <- dplyr::filter(df, !is.na(.data$.cat), .data$.cat != "NA", .data$.cat != "")
  if (is.null(category_levels)) category_levels <- sort(unique(df$.cat))

  tab <- dplyr::count(df, .data$.sample, .data$.cond, .data$.cat, name = "k")
  if (condition_level == "sample") {
    tab <- tab |>
      dplyr::group_by(.data$.sample, .data$.cond) |>
      tidyr::complete(.cat = category_levels, fill = list(k = 0)) |>
      dplyr::mutate(
        n = sum(.data$k),
        prop = ifelse(.data$n > 0, .data$k / .data$n, NA_real_)
      ) |>
      dplyr::ungroup()
  } else {
    condition_levels <- sort(unique(df$.cond))
    tab <- tab |>
      dplyr::group_by(.data$.sample) |>
      tidyr::complete(
        .cond = condition_levels,
        .cat = category_levels,
        fill = list(k = 0)
      ) |>
      dplyr::group_by(.data$.sample, .data$.cond) |>
      dplyr::mutate(
        n = sum(.data$k),
        prop = ifelse(.data$n > 0, .data$k / .data$n, NA_real_)
      ) |>
      dplyr::ungroup()
  }

  tab <- tab |>
    dplyr::rename(
      !!sample_col := .data$.sample,
      !!condition_col := .data$.cond,
      !!category_col := .data$.cat
    ) |>
    dplyr::mutate(
      !!category_col := factor(.data[[category_col]], levels = category_levels)
    )

  strata_n <- tab |>
    dplyr::group_by(.data[[sample_col]], .data[[condition_col]]) |>
    dplyr::summarise(n = max(.data$n, na.rm = TRUE), .groups = "drop")
  if (drop_empty_strata) {
    keep_strata <- dplyr::filter(strata_n, !is.na(.data$n), .data$n > 0)
    tab <- dplyr::inner_join(
      tab,
      dplyr::select(keep_strata, dplyr::all_of(c(sample_col, condition_col))),
      by = c(sample_col, condition_col)
    )
  }

  if (condition_level == "spot" && paired_only) {
    condition_levels <- sort(unique(as.character(tab[[condition_col]])))
    if (is.null(control_level) || is.null(case_level)) {
      if (all(c("Normal", "Low-grade") %in% condition_levels)) {
        control_level <- "Normal"
        case_level <- "Low-grade"
      } else {
        control_level <- condition_levels[[1]]
        case_level <- condition_levels[[2]]
      }
    }
    have_both <- strata_n |>
      dplyr::filter(
        .data[[condition_col]] %in% c(control_level, case_level),
        !is.na(.data$n),
        .data$n > 0
      ) |>
      dplyr::count(.data[[sample_col]], name = "n_conditions") |>
      dplyr::filter(.data$n_conditions == 2) |>
      dplyr::select(dplyr::all_of(sample_col))
    tab <- dplyr::semi_join(tab, have_both, by = sample_col)
  }

  wide <- tab |>
    dplyr::select(
      dplyr::all_of(sample_col),
      dplyr::all_of(condition_col),
      dplyr::all_of(category_col),
      dplyr::all_of("prop")
    ) |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(category_col),
      values_from = "prop",
      values_fill = 0
    )
  condition_levels_present <- sort(unique(as.character(wide[[condition_col]])))
  if (length(condition_levels_present) != 2) {
    stop(
      "After filtering, condition levels != 2. Current levels: ",
      paste(condition_levels_present, collapse = ", "),
      ". Consider paired_only=FALSE or na_category_label='Unclassified'."
    )
  }

  composition_matrix <- as.matrix(
    dplyr::select(wide, -dplyr::all_of(c(sample_col, condition_col)))
  )
  X_hell <- sqrt(composition_matrix)
  permanova <- vegan::adonis2(
    X_hell ~ wide[[condition_col]],
    permutations = n_perm
  )
  distance <- vegan::vegdist(X_hell, method = "euclidean")
  betadisper <- vegan::betadisper(distance, group = wide[[condition_col]])
  betadisper_permutest <- vegan::permutest(betadisper, permutations = n_perm)

  if (condition_level == "sample") {
    group_1 <- condition_levels_present[[1]]
    group_2 <- condition_levels_present[[2]]
    posthoc <- tab |>
      dplyr::group_by(.data[[category_col]]) |>
      dplyr::summarise(
        p = stats::wilcox.test(
          .data$prop[.data[[condition_col]] == group_1],
          .data$prop[.data[[condition_col]] == group_2]
        )$p.value,
        diff_med = stats::median(
          .data$prop[.data[[condition_col]] == group_2],
          na.rm = TRUE
        ) - stats::median(
          .data$prop[.data[[condition_col]] == group_1],
          na.rm = TRUE
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(p_adj = stats::p.adjust(.data$p, method = "BH")) |>
      dplyr::arrange(.data$p_adj)
  } else {
    posthoc <- lapply(levels(tab[[category_col]]), function(category) {
      paired <- tab |>
        dplyr::filter(.data[[category_col]] == category) |>
        dplyr::select(
          dplyr::all_of(sample_col),
          dplyr::all_of(condition_col),
          dplyr::all_of("prop")
        ) |>
        tidyr::pivot_wider(
          names_from = dplyr::all_of(condition_col),
          values_from = "prop",
          values_fill = 0
        )
      case <- paired[[case_level]]
      control <- paired[[control_level]]
      test <- stats::wilcox.test(case, control, paired = TRUE)
      tibble::tibble(
        !!category_col := category,
        p = test$p.value,
        diff_med = stats::median(case - control, na.rm = TRUE)
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(p_adj = stats::p.adjust(.data$p, method = "BH")) |>
      dplyr::arrange(.data$p_adj)
  }

  list(
    tab = tab,
    wide = wide,
    X_hell = X_hell,
    permanova = permanova,
    betadisper = betadisper,
    betadisper_permutest = betadisper_permutest,
    posthoc = posthoc,
    condition_levels = condition_levels_present
  )
}
