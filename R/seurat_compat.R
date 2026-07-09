# Internal Seurat v4/v5 compatibility helpers

.sc_seurat_v5 <- function() {
  utils::packageVersion("SeuratObject") >= "5.0.0" || utils::packageVersion("Seurat") >= "5.0.0"
}

.sc_assay_object <- function(object, assay = NULL) {
  if (inherits(object, what = c("Assay", "Assay5", "StdAssay", "ChromatinAssay"))) {
    return(object)
  }
  assay <- assay %||% Seurat::DefaultAssay(object)
  object[[assay]]
}

.sc_uses_layers <- function(object, assay = NULL) {
  .sc_seurat_v5() && inherits(.sc_assay_object(object, assay), what = c("Assay5", "StdAssay"))
}

.sc_get_assay_data <- function(object, assay = NULL, slot = "data") {
  args <- list(object = object)
  if (!is.null(assay)) args$assay <- assay
  if (.sc_uses_layers(object, assay)) {
    args$layer <- slot
  } else {
    args$slot <- slot
  }
  do.call(Seurat::GetAssayData, args)
}

.sc_set_assay_data <- function(object, new.data, assay = NULL, slot = "data") {
  args <- list(object = object, new.data = new.data)
  if (!is.null(assay)) args$assay <- assay
  if (.sc_uses_layers(object, assay)) {
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
    if (.sc_uses_layers(object)) args$layer <- slot else args$slot <- slot
  }
  do.call(Seurat::FetchData, args)
}

.sc_create_assay_object <- function(counts = NULL, data = NULL, assay.version = "v4") {
  args <- list()
  if (!is.null(counts)) args$counts <- counts
  if (!is.null(data)) args$data <- data
  if (identical(assay.version, "v5") && .sc_seurat_v5() &&
      exists("CreateAssay5Object", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    do.call(SeuratObject::CreateAssay5Object, args)
  } else {
    do.call(Seurat::CreateAssayObject, args)
  }
}

.sc_create_seurat_object <- function(..., assay.version = "v4") {
  if (!.sc_seurat_v5()) {
    return(Seurat::CreateSeuratObject(...))
  }
  old.option <- getOption("Seurat.object.assay.version")
  on.exit(options(Seurat.object.assay.version = old.option), add = TRUE)
  options(Seurat.object.assay.version = assay.version)
  Seurat::CreateSeuratObject(...)
}

.sc_assay_has_slot <- function(object, assay, slot = "counts") {
  data <- .sc_get_assay_data(object = object, assay = assay, slot = slot)
  !SeuratObject::IsMatrixEmpty(data)
}

.sc_layers <- function(object, search = NULL) {
  if (.sc_uses_layers(object) && exists("Layers", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    SeuratObject::Layers(object = object, search = search)
  } else {
    search[search %in% c("counts", "data", "scale.data")]
  }
}

.SCFilterAssays <- function(object, classes.keep = c("Assay", "Assay5", "StdAssay", "ChromatinAssay")) {
  ns_so <- asNamespace("SeuratObject")
  ns_seurat <- asNamespace("Seurat")
  filter_objects <- NULL

  if (exists(".FilterObjects", envir = ns_so, mode = "function", inherits = FALSE)) {
    filter_objects <- get(".FilterObjects", envir = ns_so, inherits = FALSE)
  } else if (exists(".FilterObjects", envir = ns_seurat, mode = "function", inherits = FALSE)) {
    filter_objects <- get(".FilterObjects", envir = ns_seurat, inherits = FALSE)
  } else if (exists("FilterObjects", envir = ns_so, mode = "function", inherits = FALSE)) {
    filter_objects <- get("FilterObjects", envir = ns_so, inherits = FALSE)
  }

  if (!is.null(filter_objects)) {
    filtered <- tryCatch(
      filter_objects(object = object, classes.keep = classes.keep),
      error = function(e) NULL
    )
    if (!is.null(filtered)) return(filtered)
  }

  assays <- Seurat::Assays(object = object)
  assays[vapply(assays, function(assay) {
    inherits(object[[assay]], classes.keep)
  }, logical(1))]
}

.SCResolveNbCl <- function(nb_cl = NULL) {
  if (!is.null(nb_cl)) {
    nb <- suppressWarnings(as.integer(nb_cl[1]))
    if (is.na(nb) || nb < 1L) nb <- 1L
    return(nb)
  }

  env <- Sys.getenv(c("SUPER_CELL_WORKERS", "REGCOMPASS_WORKERS", "SLURM_CPUS_PER_TASK", "NSLOTS"), unset = NA)
  env <- suppressWarnings(as.integer(env))
  env <- env[is.finite(env) & env >= 1L]
  if (length(env)) return(env[[1]])

  if (requireNamespace("future", quietly = TRUE)) {
    fc <- tryCatch(future::availableCores(), error = function(e) NA_integer_)
    if (is.finite(fc) && fc >= 1L) return(as.integer(fc))
  }

  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || cores < 1L) cores <- 1L
  max(1L, as.integer(cores) - 2L)
}

.SCFormatMetacellNames <- function(membership, prefixMC = "") {
  out <- paste0(prefixMC %||% "", membership)
  names(out) <- names(membership)
  as.character(out)
}

.sc_filter_assays <- .SCFilterAssays
