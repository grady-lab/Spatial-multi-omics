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

run_cytotrace2 <- function(object, assay = "Spatial", ncores = 1) {
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
    components = c(1, 2),
    name = "diffmap"
) {
  embedding <- Seurat::Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  dm <- destiny::DiffusionMap(embedding, sigma = sigma)
  dc <- destiny::eigenvectors(dm)[, components, drop = FALSE]
  colnames(dc) <- paste0("DC", seq_len(ncol(dc)))
  rownames(dc) <- rownames(embedding)
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
