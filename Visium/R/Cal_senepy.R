# Clean, project-portable version of Scripts/Utils/Cal_senepy.R.
Calculate_senepy <- function(
    object,
    assay = "Spatial",
    layer = "data",
    species = "Human",
    tissue = "intestine",
    cell_type = "epithelial cell",
    epithelial_modules = 0:1,
    conda_env = "adenoma-senepy",
    binarize = FALSE,
    importance = TRUE
) {
  reticulate::use_condaenv(conda_env, required = TRUE)

  joined <- SeuratObject::JoinLayers(object, assay = assay)
  expression <- SeuratObject::LayerData(joined, assay = assay, layer = layer)
  expression_py <- reticulate::r_to_py(Matrix::t(expression), convert = FALSE)

  anndata <- reticulate::import("anndata", convert = FALSE)
  pandas <- reticulate::import("pandas", convert = FALSE)
  senepy <- reticulate::import("senepy", convert = FALSE)

  adata <- anndata$AnnData(
    X = expression_py,
    obs = reticulate::r_to_py(object@meta.data),
    var = pandas$DataFrame(index = reticulate::r_to_py(rownames(expression)))
  )
  adata$obs_names <- reticulate::r_to_py(colnames(expression))

  hubs <- senepy$load_hubs(species = species)
  hubs$merge_hubs(
    hubs$metadata,
    new_name = "Universal",
    calculate_thresh = TRUE,
    p_thres = 0.01
  )
  translator <- senepy$translator(hub = hubs$hubs, data = adata)

  score_names <- paste0("epithelial_cell_score_", epithelial_modules)
  for (i in seq_along(epithelial_modules)) {
    module_key <- reticulate::tuple(
      tissue,
      cell_type,
      as.integer(epithelial_modules[[i]])
    )
    adata$obs[score_names[[i]]] <- senepy$score_hub(
      adata,
      hubs$hubs[[module_key]],
      translator = translator,
      binarize = binarize,
      importance = importance
    )
  }
  adata$obs["Universal_score"] <- senepy$score_hub(
    adata,
    hubs$hubs[["Universal"]],
    translator = translator,
    binarize = binarize,
    importance = importance
  )

  scores <- reticulate::py_to_r(adata$obs)[, c(score_names, "Universal_score"), drop = FALSE]
  SeuratObject::AddMetaData(object, scores[colnames(object), , drop = FALSE])
}

