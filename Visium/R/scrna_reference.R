marteau_reference_excluded_annotations <- function() {
  c(
    "Dysplastic_Senescent_metaplastic",
    "Dysplastic_Enteroendocrine_like",
    "B cell naive",
    "Plasma IgM",
    "Plasmablast",
    "Schwann cell",
    "DC mature",
    "DC3",
    "Pericyte",
    "pDC",
    "Monocyte non-classical",
    "ILC"
  )
}

curate_marteau_reference <- function(object, study_pmid = "34910928") {
  metadata <- object@meta.data
  assert_columns(
    metadata,
    c(
      "study_pmid",
      "study_id",
      "sample_id",
      "cell_type_fine",
      "my_annotation",
      "my_annotationV2",
      "_index",
      "nCount_RNA",
      "nFeature_RNA"
    ),
    "Marteau metadata"
  )

  cells <- rownames(metadata)[
    !is.na(metadata$study_pmid) & as.character(metadata$study_pmid) == study_pmid
  ]
  if (!length(cells)) stop("No cells found for study_pmid ", study_pmid)
  reference <- subset(object, cells = cells)
  metadata <- reference@meta.data

  annotation <- as.character(metadata$my_annotation)
  missing_annotation <- is.na(annotation) | annotation == ""
  annotation[missing_annotation] <- as.character(
    metadata$cell_type_fine[missing_annotation]
  )

  replace_matching <- function(pattern, replacement) {
    matched <- !is.na(annotation) & grepl(pattern, annotation)
    annotation[matched] <<- replacement
  }
  replace_matching("Fibro", "Fibroblast")
  replace_matching("Endothelial", "Endothelial")
  replace_matching("CD4 cycling", "CD4.T")
  replace_matching("CD8 cycling", "CD8.T")
  replace_matching("Macrophage cycling", "Macrophage_Monocyte")
  replace_matching("cDC", "cDC")
  replace_matching("Monocyte", "Macrophage_Monocyte")

  keep <- !(annotation %in% marteau_reference_excluded_annotations())
  reference <- subset(reference, cells = rownames(metadata)[keep])
  metadata <- reference@meta.data
  annotation <- annotation[keep]

  dysplastic <- !is.na(annotation) & grepl("Dysplastic", annotation)
  refined_annotation <- as.character(metadata$my_annotationV2)
  if (any(is.na(refined_annotation[dysplastic]) | refined_annotation[dysplastic] == "")) {
    stop("Dysplastic cells are missing the curated my_annotationV2 label")
  }
  annotation[dysplastic] <- refined_annotation[dysplastic]

  keep <- !is.na(annotation) & annotation != "Dysplastic_Normal"
  reference <- subset(reference, cells = rownames(metadata)[keep])
  reference$my_annotation <- make.names(annotation[keep])
  reference
}

expected_marteau_reference_counts <- function() {
  c(
    B.cell.activated = 66L,
    CD4.T = 168L,
    CD8.T = 923L,
    Dysplastic_Absorptive_like = 2642L,
    Dysplastic_Paneth_like = 286L,
    Dysplastic_Proliferating = 1702L,
    Dysplastic_Senescent_like = 2105L,
    Dysplastic_Stem_like = 5493L,
    Dysplastic_Tuft_like = 1367L,
    Endothelial = 65L,
    Enteroendocrine = 119L,
    Fibroblast = 74L,
    GC.B.cell = 21L,
    Macrophage = 76L,
    Macrophage_Monocyte = 21L,
    Mast.cell = 160L,
    NK = 28L,
    Normal.Colonocyte = 18529L,
    Normal.Colonocyte.BEST4 = 667L,
    Normal.Crypt.cell = 2171L,
    Normal.Goblet = 10759L,
    Normal.TA = 4007L,
    Plasma.IgA = 820L,
    Plasma.IgG = 41L,
    Treg = 39L,
    Tuft = 428L,
    cDC = 68L,
    gamma.delta = 26L
  )
}

validate_marteau_reference <- function(reference) {
  metadata <- reference@meta.data
  assert_columns(metadata, c("my_annotation", "sample_id", "study_id"))
  expected <- expected_marteau_reference_counts()
  observed <- table(as.character(metadata$my_annotation))
  observed_aligned <- stats::setNames(
    as.integer(observed[names(expected)]),
    names(expected)
  )

  unexpected <- setdiff(names(observed), names(expected))
  if (
    ncol(reference) != 52871L ||
      nrow(reference) != 24593L ||
      length(unique(metadata$sample_id)) != 66L ||
      !identical(unique(as.character(metadata$study_id)), "Chen_2021_Cell") ||
      length(unexpected) ||
      anyNA(observed_aligned) ||
      !identical(unname(observed_aligned), unname(expected))
  ) {
    stop(
      "The curated Marteau reference does not match the realized Cell2location ",
      "reference (52,871 cells, 24,593 genes, 66 samples, and 28 annotations)"
    )
  }

  tibble::tibble(
    my_annotation = names(expected),
    n_cells = unname(expected)
  )
}

export_marteau_reference_h5ad <- function(
    reference,
    path,
    assay = "RNA"
) {
  validate_marteau_reference(reference)
  counts <- Seurat::GetAssayData(reference, assay = assay, layer = "counts")
  metadata <- reference@meta.data |>
    dplyr::transmute(
      nCount_RNA = .data$nCount_RNA,
      nFeature_RNA = .data$nFeature_RNA,
      original_barcode_index = .data[["_index"]],
      sample_id = .data$sample_id,
      my_annotation = .data$my_annotation,
      study_id = .data$study_id
    )

  anndata <- reticulate::import("anndata", convert = FALSE)
  pandas <- reticulate::import("pandas", convert = FALSE)
  counts_dense <- as.matrix(t(counts))
  adata <- anndata$AnnData(
    X = reticulate::r_to_py(counts_dense),
    obs = reticulate::r_to_py(metadata),
    var = pandas$DataFrame(index = reticulate::r_to_py(rownames(counts)))
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  adata$write_h5ad(normalizePath(path, mustWork = FALSE))
  invisible(path)
}
