#' Validate a SuperCell2 metacell Seurat output
#'
#' @param object A Seurat metacell object returned by SuperCell2.
#' @param rna_assay Optional RNA assay name retained for downstream callers.
#' @param atac_assay Optional ATAC assay name retained for downstream callers.
#' @param require_fragments Whether to require fragment files and indexes from the manifest.
#' @param require_complete_fragment_coverage Whether every object cell must be represented in the fragment cell map when fragments are required.
#' @return Invisibly returns TRUE when the output contract is valid.
#' @export
validate_metacell_output <- function(object,
                                     rna_assay = NULL,
                                     atac_assay = NULL,
                                     require_fragments = FALSE,
                                     require_complete_fragment_coverage = FALSE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }

  membership_table <- tryCatch(object@misc$membership_table, error = function(e) NULL)
  if (!is.data.frame(membership_table) || !all(c("cell_id", "metacell_id") %in% colnames(membership_table))) {
    stop("Missing valid `misc$membership_table`.", call. = FALSE)
  }

  mc_ids <- as.character(colnames(object))
  if (!setequal(unique(membership_table$metacell_id), mc_ids)) {
    stop("Membership metacell IDs do not match object colnames.", call. = FALSE)
  }

  if (anyDuplicated(membership_table$cell_id)) {
    stop("Duplicated single-cell IDs in membership table.", call. = FALSE)
  }

  assays <- names(object@assays)
  for (assay in assays) {
    assay_ids <- as.character(colnames(.sc_get_assay_data(object, assay = assay, slot = "counts")))
    if (!identical(assay_ids, mc_ids)) {
      stop("Assay column order differs from object metacell IDs: ", assay, call. = FALSE)
    }
  }

  if (isTRUE(require_fragments)) {
    manifest <- object@misc$fragment_manifest
    if (!is.data.frame(manifest) || nrow(manifest) == 0L) {
      stop("Fragment manifest is missing.", call. = FALSE)
    }
    required_manifest_cols <- c("fragment_file", "index_file", "n_metacells")
    if (!all(required_manifest_cols %in% colnames(manifest))) {
      stop("Fragment manifest is missing required columns: ",
           paste(setdiff(required_manifest_cols, colnames(manifest)), collapse = ", "),
           call. = FALSE)
    }
    if (any(!nzchar(manifest$fragment_file)) || any(!nzchar(manifest$index_file)) ||
        any(!file.exists(manifest$fragment_file)) || any(!file.exists(manifest$index_file)) ||
        any(file.info(manifest$fragment_file)$size == 0) || any(file.info(manifest$index_file)$size == 0)) {
      stop("One or more fragment files or indexes are missing or empty.", call. = FALSE)
    }

    fragment_cell_map <- object@misc$fragment_cell_map
    required_map_cols <- c("fragment_file", "object_cell", "fragment_barcode")
    if (!is.data.frame(fragment_cell_map) || nrow(fragment_cell_map) == 0L) {
      stop("Fragment cell map is missing.", call. = FALSE)
    }
    if (!all(required_map_cols %in% colnames(fragment_cell_map))) {
      stop("Fragment cell map is missing required columns: ",
           paste(setdiff(required_map_cols, colnames(fragment_cell_map)), collapse = ", "),
           call. = FALSE)
    }
    fragment_cell_map$fragment_file <- as.character(fragment_cell_map$fragment_file)
    fragment_cell_map$object_cell <- as.character(fragment_cell_map$object_cell)
    fragment_cell_map$fragment_barcode <- as.character(fragment_cell_map$fragment_barcode)

    unknown_cells <- setdiff(fragment_cell_map$object_cell, mc_ids)
    if (length(unknown_cells)) {
      stop("Fragment cell map contains object cells absent from the Seurat object: ",
           paste(utils::head(unknown_cells, 10L), collapse = ", "), call. = FALSE)
    }
    if (anyDuplicated(fragment_cell_map$object_cell)) {
      stop("One or more object cells are assigned to multiple fragment files.", call. = FALSE)
    }
    if (isTRUE(require_complete_fragment_coverage) && !setequal(fragment_cell_map$object_cell, mc_ids)) {
      stop("Fragment cell map does not cover all object cells.", call. = FALSE)
    }

    manifest_files <- as.character(manifest$fragment_file)
    mapped_files <- unique(fragment_cell_map$fragment_file)
    missing_from_map <- setdiff(manifest_files, mapped_files)
    if (length(missing_from_map)) {
      stop("One or more manifest fragment files have no cell map rows: ",
           paste(utils::head(missing_from_map, 10L), collapse = ", "), call. = FALSE)
    }
    missing_from_manifest <- setdiff(mapped_files, manifest_files)
    if (length(missing_from_manifest)) {
      stop("Fragment cell map contains files absent from the manifest: ",
           paste(utils::head(missing_from_manifest, 10L), collapse = ", "), call. = FALSE)
    }
    mapped_counts <- tapply(fragment_cell_map$object_cell, fragment_cell_map$fragment_file,
                            function(x) length(unique(x)))
    expected_counts <- stats::setNames(as.integer(manifest$n_metacells), manifest_files)
    for (fragment_file in manifest_files) {
      if (!identical(unname(mapped_counts[[fragment_file]]), unname(expected_counts[[fragment_file]]))) {
        stop("Fragment manifest n_metacells does not match fragment cell map for: ", fragment_file, call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}
