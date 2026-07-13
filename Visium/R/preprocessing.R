assert_columns <- function(x, columns, object_name = deparse(substitute(x))) {
  missing <- setdiff(columns, colnames(x))
  if (length(missing)) {
    stop(object_name, " is missing: ", paste(missing, collapse = ", "))
  }
  invisible(x)
}

read_sample_manifest <- function(path) {
  manifest <- readr::read_csv(path, show_col_types = FALSE)
  assert_columns(manifest, c("sample_id", "space_ranger_outs", "has_adt", "include"))
  manifest |>
    dplyr::filter(.data$include) |>
    dplyr::mutate(
      sample_id = as.character(.data$sample_id),
      space_ranger_outs = as.character(.data$space_ranger_outs),
      has_adt = as.logical(.data$has_adt)
    )
}

.layer_or_null <- function(object, assay, layer) {
  if (!layer %in% SeuratObject::Layers(object[[assay]])) return(NULL)
  SeuratObject::LayerData(object, assay = assay, layer = layer)
}

split_rna_adt <- function(object, has_adt) {
  spatial_layers <- SeuratObject::Layers(object[["Spatial"]])
  rna_layer <- if ("counts.Gene Expression" %in% spatial_layers) {
    "counts.Gene Expression"
  } else {
    grep("^counts", spatial_layers, value = TRUE)[1]
  }
  rna <- SeuratObject::LayerData(object, assay = "Spatial", layer = rna_layer)
  adt <- NULL
  if (isTRUE(has_adt)) {
    adt_layer <- "counts.Antibody Capture"
    if (!adt_layer %in% spatial_layers) {
      stop("Sample ", unique(object$Sample_ID), " is marked has_adt but has no antibody layer")
    }
    adt <- .layer_or_null(object, "Spatial", adt_layer)
  }

  object[["Spatial"]] <- SeuratObject::CreateAssay5Object(counts = rna)
  if (isTRUE(has_adt)) {
    adt_names <- sub("\\.1$", "", rownames(adt))
    adt_names <- make.unique(paste0("ADT-", sub("^ADT[-_]", "", adt_names)))
    rownames(adt) <- adt_names
    object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt)
  }
  object
}

load_visium_samples <- function(manifest) {
  assert_columns(manifest, c("sample_id", "space_ranger_outs", "has_adt"))
  missing <- manifest$space_ranger_outs[!dir.exists(manifest$space_ranger_outs)]
  if (length(missing)) {
    stop("Missing Space Ranger directories:\n", paste(missing, collapse = "\n"))
  }

  objects <- purrr::pmap(
    manifest[c("sample_id", "space_ranger_outs", "has_adt")],
    function(sample_id, space_ranger_outs, has_adt) {
      object <- Seurat::Load10X_Spatial(
        data.dir = space_ranger_outs,
        filename = "filtered_feature_bc_matrix.h5",
        assay = "Spatial",
        slice = sample_id,
        filter.matrix = TRUE
      )
      object$Sample_ID <- sample_id
      object$LargeSmall <- factor(
        ifelse(grepl("large", sample_id, ignore.case = TRUE), "Large", "Small"),
        levels = c("Small", "Large")
      )
      object$has_ADT <- has_adt
      object <- split_rna_adt(object, has_adt)
      object$percent.mt <- Seurat::PercentageFeatureSet(object, assay = "Spatial", pattern = "^MT-")
      object
    }
  )
  names(objects) <- manifest$sample_id
  objects
}

qc_rna <- function(object, min_umi = 300, max_mito = 40) {
  Seurat::DefaultAssay(object) <- "Spatial"
  counts <- SeuratObject::LayerData(object, assay = "Spatial", layer = "counts")
  object$nCount_Spatial <- Matrix::colSums(counts)
  object$nFeature_Spatial <- Matrix::colSums(counts > 0)
  keep <- object$nCount_Spatial > min_umi & object$percent.mt < max_mito
  subset(object, cells = colnames(object)[keep])
}

add_pathology_annotations <- function(objects, annotation_file) {
  annotations <- readr::read_csv(annotation_file, show_col_types = FALSE)
  assert_columns(annotations, c("sample_id", "barcode", "Adenoma"))

  purrr::imap(objects, function(object, sample_id) {
    ann <- annotations |>
      dplyr::filter(.data$sample_id == sample_id) |>
      dplyr::select(.data$barcode, .data$Adenoma)
    raw_barcode <- sub("^.*_([ACGT]+-[0-9]+)$", "\\1", colnames(object))
    idx <- match(raw_barcode, ann$barcode)
    object$Adenoma <- ann$Adenoma[idx]
    object
  })
}

merge_visium_samples <- function(objects) {
  if (!length(objects)) stop("No Visium objects supplied")
  Seurat::merge(
    x = objects[[1]],
    y = objects[-1],
    add.cell.ids = names(objects),
    merge.data = FALSE
  )
}

prepare_adt <- function(object, min_antibodies = 20, max_isotype_fraction = 0.20) {
  adt_object <- SeuratObject::JoinLayers(object, assay = "ADT")
  counts <- SeuratObject::LayerData(adt_object, assay = "ADT", layer = "counts")
  isotypes <- grep("IgG", rownames(counts), value = TRUE, ignore.case = TRUE)
  targets <- setdiff(rownames(counts), isotypes)
  if (!length(isotypes)) stop("No IgG isotype controls were found in the ADT assay")

  detected <- Matrix::colSums(counts[targets, , drop = FALSE] > 0)
  total <- Matrix::colSums(counts)
  isotype_fraction <- Matrix::colSums(counts[isotypes, , drop = FALSE]) / pmax(total, 1)
  keep <- detected >= min_antibodies & isotype_fraction < max_isotype_fraction

  object$ADT_n_detected <- detected[colnames(object)]
  object$ADT_isotype_fraction <- isotype_fraction[colnames(object)]
  object$ADT_pass <- keep[colnames(object)]

  list(
    object = object,
    keep_cells = names(keep)[keep],
    targets = targets,
    isotypes = isotypes
  )
}

build_multimodal_reference <- function(
    merged,
    rna_dims = 30,
    adt_dims = 15,
    resolution = 0.4,
    seed = 1
) {
  set.seed(seed)
  adt_cells <- colnames(merged)[merged$has_ADT %in% TRUE]
  adt_all <- subset(merged, cells = adt_cells)
  adt_qc <- prepare_adt(adt_all)
  reference <- subset(adt_qc$object, cells = adt_qc$keep_cells)

  Seurat::DefaultAssay(reference) <- "Spatial"
  reference <- reference |>
    Seurat::NormalizeData(verbose = FALSE) |>
    Seurat::FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
    Seurat::ScaleData(verbose = FALSE) |>
    Seurat::RunPCA(npcs = rna_dims, reduction.name = "pca", verbose = FALSE)
  reference <- Seurat::IntegrateLayers(
    reference,
    method = SeuratWrappers::FastMNNIntegration,
    orig.reduction = "pca",
    new.reduction = "mnn",
    verbose = FALSE
  )

  adt_joined <- SeuratObject::JoinLayers(reference, assay = "ADT")
  adt_counts <- SeuratObject::LayerData(adt_joined, assay = "ADT", layer = "counts")
  cells_by_slide <- split(colnames(adt_counts), adt_joined$Sample_ID)
  adt_dsb <- lapply(cells_by_slide, function(cells) {
    dsb::ModelNegativeADTnorm(
      cell_protein_matrix = as.matrix(adt_counts[, cells, drop = FALSE]),
      denoise.counts = TRUE,
      use.isotype.control = TRUE,
      isotype.control.name.vec = adt_qc$isotypes,
      fast.km = TRUE
    )
  })
  adt_dsb <- do.call(cbind, adt_dsb)
  adt_dsb <- adt_dsb[, colnames(reference), drop = FALSE]
  reference[["ADT_dsb"]] <- SeuratObject::CreateAssay5Object(
    counts = adt_counts[, colnames(reference), drop = FALSE],
    data = adt_dsb
  )

  Seurat::DefaultAssay(reference) <- "ADT_dsb"
  SeuratObject::VariableFeatures(reference) <- adt_qc$targets
  reference <- Seurat::ScaleData(reference, features = adt_qc$targets, verbose = FALSE)
  reference <- Seurat::RunPCA(
    reference,
    features = adt_qc$targets,
    npcs = adt_dims,
    reduction.name = "apca",
    verbose = FALSE
  )
  reference[["ADT_dsb"]] <- split(reference[["ADT_dsb"]], f = reference$Sample_ID)
  reference <- Seurat::IntegrateLayers(
    reference,
    method = SeuratWrappers::HarmonyIntegration,
    orig.reduction = "apca",
    new.reduction = "adt_harmony",
    verbose = FALSE
  )

  reference <- Seurat::FindMultiModalNeighbors(
    reference,
    reduction.list = list("mnn", "adt_harmony"),
    dims.list = list(seq_len(rna_dims), seq_len(adt_dims)),
    modality.weight.name = c("RNA.weight", "ADT.weight")
  )
  reference <- Seurat::RunUMAP(
    reference,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    return.model = TRUE,
    seed.use = seed
  )
  reference <- Seurat::FindClusters(
    reference,
    graph.name = "wsnn",
    algorithm = 3,
    resolution = resolution,
    random.seed = seed,
    verbose = FALSE
  )
  reference$combined_cluster <- as.character(Seurat::Idents(reference))

  query_cells <- setdiff(colnames(merged), colnames(reference))
  list(reference = reference, query = subset(merged, cells = query_cells))
}

map_rna_only_spots <- function(reference, query, dims = 30, seed = 1) {
  set.seed(seed)
  Seurat::DefaultAssay(reference) <- "Spatial"
  Seurat::DefaultAssay(query) <- "Spatial"
  features <- SeuratObject::VariableFeatures(reference)

  query <- query |>
    Seurat::NormalizeData(verbose = FALSE) |>
    Seurat::ScaleData(features = features, verbose = FALSE) |>
    Seurat::RunPCA(features = features, npcs = dims, verbose = FALSE)

  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = query,
    features = features,
    reference.reduction = "pca",
    reduction = "rpca",
    dims = seq_len(dims),
    normalization.method = "LogNormalize"
  )
  mapped <- Seurat::MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = list(combined_cluster = "combined_cluster"),
    reference.reduction = "pca",
    reference.dims = seq_len(dims),
    reduction.model = "wnn.umap",
    new.reduction.name = "ref.pca"
  )
  mapped$combined_cluster <- as.character(mapped$predicted.combined_cluster)
  mapped[["wnn.umap"]] <- mapped[["ref.umap"]]
  mapped
}

merge_reference_and_query <- function(reference, mapped_query) {
  final <- Seurat::merge(reference, mapped_query, merge.data = FALSE)
  embeddings <- rbind(
    Seurat::Embeddings(reference, "wnn.umap"),
    Seurat::Embeddings(mapped_query, "wnn.umap")
  )
  embeddings <- embeddings[colnames(final), , drop = FALSE]
  final[["wnn.umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "wnnUMAP_",
    assay = "Spatial"
  )
  final@images <- c(reference@images, mapped_query@images[setdiff(names(mapped_query@images), names(reference@images))])
  Seurat::Idents(final) <- factor(final$combined_cluster)
  SeuratObject::UpdateSeuratObject(final)
}
