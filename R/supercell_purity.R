#' Compute purity of super-cells
#'
#'
#' @param clusters vector of clustering assignment (reference assignment)
#' @param supercell_membership vector of assignment of single-cell data to super-cells (membership field of \link{SCimplify} function output)
#' @param method method to compute super-cell purity.
#' \code{"max_proportion"} if the purity is defined as a proportion of the most abundant cluster (cell type) within super-cell or
#' \code{"entropy"} if the purity is defined as the Shanon entropy of the cell types super-cell consists of.
#'
#' @return a named vector of super-cell purity. Super-cells whose cells all have
#'   missing cluster labels are retained and assigned `NA`.
#'
#' @export
#'

supercell_purity <- function(
  clusters,
  supercell_membership,
  method = c("max_proportion", "entropy")[1]
){
  if(!(method %in% c("max_proportion", "entropy"))){
    stop(paste("Method", method, "is not known. The available methods are: max_proportion, entropy"))
  }
  if (length(clusters) != length(supercell_membership)) {
    stop("`clusters` and `supercell_membership` must have the same length.", call. = FALSE)
  }

  membership_values <- as.character(supercell_membership)
  membership_levels <- sort(unique(membership_values[!is.na(membership_values) & nzchar(membership_values)]))
  if (!length(membership_levels)) {
    return(stats::setNames(numeric(), character()))
  }
  membership_factor <- factor(membership_values, levels = membership_levels)

  cluster_values <- as.character(clusters)
  cluster_levels <- sort(unique(cluster_values[!is.na(cluster_values)]))
  result <- stats::setNames(rep(NA_real_, length(membership_levels)), membership_levels)
  if (!length(cluster_levels)) {
    return(result)
  }
  cluster_factor <- factor(cluster_values, levels = cluster_levels)
  cl.gr <- table(cluster_factor, membership_factor)
  observed <- colSums(cl.gr) > 0
  if (!any(observed)) {
    return(result)
  }
  cl.use <- cl.gr[, observed, drop = FALSE]

  values <- switch(
    method,
    entropy = apply(cl.use, 2, entropy::entropy),
    max_proportion = {
      group.size <- colSums(cl.use)
      apply(sweep(cl.use, 2, group.size, "/"), 2, max)
    }
  )
  result[names(values)] <- as.numeric(values)
  result
}
