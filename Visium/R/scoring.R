run_irgsea <- function(object, gene_sets, group_by = "combined_cluster", assay = "Spatial") {
  methods <- c("AUCell", "UCell", "singscore", "JASMINE")
  object <- irGSEA::irGSEA.score(
    object = object,
    assay = assay,
    slot = "data",
    seeds = 1,
    ncores = 1,
    min.cells = 3,
    custom = TRUE,
    geneset = gene_sets,
    msigdb = FALSE,
    maxGSSize = 2000,
    method = methods
  )
  integrated <- irGSEA::irGSEA.integrate(
    object = object,
    group.by = group_by,
    method = methods
  )
  object@misc$irGSEA_RRA <- integrated
  object
}

copy_ucell_scores <- function(object, score_map) {
  ucell <- SeuratObject::LayerData(object, assay = "UCell", layer = "scale.data")
  normalize_name <- function(x) toupper(gsub("_", "-", x, fixed = TRUE))
  available <- stats::setNames(rownames(ucell), normalize_name(rownames(ucell)))
  resolved <- unname(available[normalize_name(unname(score_map))])
  missing <- unname(score_map)[is.na(resolved)]
  if (length(missing)) stop("Missing UCell signatures: ", paste(missing, collapse = ", "))
  for (new_name in names(score_map)) {
    source_name <- resolved[match(new_name, names(score_map))]
    object[[new_name]] <- as.numeric(ucell[source_name, colnames(object)])
  }
  object
}

classify_senepy_outliers <- function(
    object,
    score_col = "Senepy",
    output_col = "SenePyhigh"
) {
  score <- object@meta.data[[score_col]]
  cutoff <- mean(score, na.rm = TRUE) + 2 * stats::sd(score, na.rm = TRUE)
  object[[output_col]] <- score > cutoff
  object@misc$SenePy_outlier_cutoff <- cutoff
  object
}

add_cell_cycle_scores <- function(object, assay = "Spatial") {
  Seurat::DefaultAssay(object) <- assay
  object <- Seurat::CellCycleScoring(
    object,
    s.features = Seurat::cc.genes.updated.2019$s.genes,
    g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )
  object$cell_cycle_seurat <- factor(
    object$Phase,
    levels = c("G1", "S", "G2M")
  )
  object$Phase <- NULL
  object <- ccAFv2::PredictCellCycle(
    object,
    species = "human",
    gene_id = "symbol",
    spatial = TRUE,
    threshold = 0.5,
    do_sctransform = FALSE,
    assay = assay
  )
  object
}

prepare_cytotrace2_input <- function(
    object,
    clusters = c("0", "1", "2", "3", "4", "6"),
    cluster_col = "combined_cluster",
    source_assay = "Spatial",
    output_assay = "RNA",
    reduction = "wnn.umap"
) {
  if (!(cluster_col %in% colnames(object@meta.data))) {
    stop("CytoTRACE2 cluster column is absent from metadata: ", cluster_col)
  }
  if (!(source_assay %in% names(object@assays))) {
    stop("CytoTRACE2 source assay is absent: ", source_assay)
  }
  cells <- rownames(object@meta.data)[
    as.character(object@meta.data[[cluster_col]]) %in% as.character(clusters)
  ]
  if (!length(cells)) stop("No cells remain after applying the CytoTRACE2 cluster subset")

  input <- subset(object, cells = cells)
  input[[output_assay]] <- methods::as(input[[source_assay]], "Assay")
  Seurat::DefaultAssay(input) <- output_assay
  if (!(reduction %in% names(input@reductions))) {
    stop("CytoTRACE2 reduction is absent: ", reduction)
  }
  Seurat::DietSeurat(
    input,
    assays = output_assay,
    dimreducs = reduction
  )
}

run_cytotrace2 <- function(object, assay = "RNA", ncores = 1) {
  Seurat::DefaultAssay(object) <- assay
  CytoTRACE2::cytotrace2(
    object,
    species = "human",
    is_seurat = TRUE,
    slot_type = "counts",
    ncores = ncores
  )
}

add_diffusion_map <- function(
    object,
    reduction = "mnn",
    dims = 1:20,
    sigma = 1,
    components = c(1, 3),
    orient_by = "dpt_pseudotime",
    orient_component = 1,
    name = "diffmap"
) {
  embedding <- Seurat::Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  old_matrix_warning <- getOption("Matrix.warnDeprecatedCoerce")
  options(Matrix.warnDeprecatedCoerce = 0)
  on.exit(options(Matrix.warnDeprecatedCoerce = old_matrix_warning), add = TRUE)
  dm <- destiny::DiffusionMap(embedding, sigma = sigma)
  dc <- destiny::eigenvectors(dm)[, components, drop = FALSE]
  colnames(dc) <- paste0("DC", seq_len(ncol(dc)))
  rownames(dc) <- rownames(embedding)

  if (!is.null(orient_by)) {
    if (!orient_by %in% colnames(object@meta.data)) {
      stop("Diffusion-map orientation variable is absent from metadata: ", orient_by)
    }
    if (orient_component < 1 || orient_component > ncol(dc)) {
      stop("orient_component must identify a retained diffusion component")
    }
    orientation <- stats::cor(
      dc[, orient_component],
      object@meta.data[rownames(dc), orient_by],
      use = "complete.obs"
    )
    if (is.finite(orientation) && orientation < 0) {
      dc[, orient_component] <- -dc[, orient_component]
    }
  }

  object[[name]] <- SeuratObject::CreateDimReducObject(
    embeddings = dc,
    key = "DC_",
    assay = Seurat::DefaultAssay(object)
  )
  object
}

add_metadata_by_barcode <- function(object, metadata, barcode_col = "barcode") {
  stopifnot(barcode_col %in% colnames(metadata))
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata[[barcode_col]]
  metadata[[barcode_col]] <- NULL
  common <- intersect(colnames(object), rownames(metadata))
  object <- SeuratObject::AddMetaData(object, metadata[common, , drop = FALSE])
  object
}
