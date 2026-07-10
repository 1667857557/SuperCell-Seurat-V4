#' Validate a SuperCell2 metacell Seurat output
#'
#' @param object A Seurat metacell object returned by SuperCell2.
#' @param rna_assay Optional RNA assay name retained for downstream callers.
#' @param atac_assay Optional ATAC assay name retained for downstream callers.
#' @param require_fragments Whether to require fragment files and indexes from the manifest.
#' @return Invisibly returns TRUE when the output contract is valid.
#' @export
validate_metacell_output <- function(object,
                                     rna_assay = NULL,
                                     atac_assay = NULL,
                                     require_fragments = FALSE) {
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
    if (any(!file.exists(manifest$fragment_file)) || any(!file.exists(manifest$index_file))) {
      stop("One or more fragment files or indexes are missing.", call. = FALSE)
    }
  }

  invisible(TRUE)
}
