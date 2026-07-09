# Internal Seurat v4/v5 compatibility helpers

.sc_seurat_v5 <- function() {
  utils::packageVersion("SeuratObject") >= "5.0.0" || utils::packageVersion("Seurat") >= "5.0.0"
}

.sc_get_assay_data <- function(object, assay = NULL, slot = "data") {
  args <- list(object = object)
  if (!is.null(assay)) args$assay <- assay
  if (.sc_seurat_v5()) {
    args$layer <- slot
  } else {
    args$slot <- slot
  }
  do.call(Seurat::GetAssayData, args)
}

.sc_set_assay_data <- function(object, new.data, assay = NULL, slot = "data") {
  args <- list(object = object, new.data = new.data)
  if (!is.null(assay)) args$assay <- assay
  if (.sc_seurat_v5()) {
    args$layer <- slot
    return(do.call(Seurat::SetAssayData, args))
  }
  if (!(slot %in% c("counts", "data", "scale.data"))) {
    assay <- assay %||% Seurat::DefaultAssay(object)
    object[[assay]]@misc[[slot]] <- new.data
    return(object)
  }
  args$slot <- slot
  do.call(Seurat::SetAssayData, args)
}

.sc_fetch_data <- function(object, vars, cells = NULL, slot = NULL, ...) {
  args <- list(object = object, vars = vars, ...)
  if (!is.null(cells)) args$cells <- cells
  if (!is.null(slot)) {
    if (.sc_seurat_v5()) args$layer <- slot else args$slot <- slot
  }
  do.call(Seurat::FetchData, args)
}

.sc_create_assay_object <- function(counts = NULL, data = NULL) {
  args <- list()
  if (!is.null(counts)) args$counts <- counts
  if (!is.null(data)) args$data <- data
  if (.sc_seurat_v5() && exists("CreateAssay5Object", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    do.call(SeuratObject::CreateAssay5Object, args)
  } else {
    do.call(Seurat::CreateAssayObject, args)
  }
}

.sc_assay_has_slot <- function(object, assay, slot = "counts") {
  data <- .sc_get_assay_data(object = object, assay = assay, slot = slot)
  !SeuratObject::IsMatrixEmpty(data)
}

.sc_layers <- function(object, search = NULL) {
  if (.sc_seurat_v5() && exists("Layers", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    SeuratObject::Layers(object = object, search = search)
  } else {
    search[search %in% c("counts", "data", "scale.data")]
  }
}
