#' Assign super-cells to the most aboundant cluster
#'
#'
#' @param clusters a vector of clustering assignment
#' @param supercell_membership a vector of assignment of single-cell data to super-cells (membership field of \link{SCimplify} function output)
#' @param method method to define the most abuldant cell cluster within super-cells. Available: "jaccard" (default), "relative", "absolute".
#' \itemize{
#'   \item jaccard - assignes super-cell to cluster with the maximum jaccard coefficient (recommended)
#'   \item relative - assignes super-cell to cluster with the maximum relative abundance (normalized by cluster size), may result in assignment of super-cells to poorly represented (small) cluser due to normalizetaion
#'   \item absolute - assignes super-cell to cluster with the maximum absolute abundance within super-cell, may result in disappearence of poorly represented (small) clusters
#' }
#'
#' @return a named vector of super-cell assignments. Super-cells whose cells all
#'   have missing cluster labels are retained and assigned `NA`.
#'
#' @export
#'


supercell_assign <- function(clusters, supercell_membership, method = c("jaccard", "relative", "absolute")){
  method <- method[[1L]]
  if (is.null(method) || is.na(method) || !(method %in% c("jaccard", "relative", "absolute"))) {
    stop("Unknown assignment method; use jaccard, relative or absolute.", call. = FALSE)
  }
  if (length(clusters) != length(supercell_membership)) {
    stop("`clusters` and `supercell_membership` must have the same length.", call. = FALSE)
  }

  membership_values <- as.character(supercell_membership)
  membership_levels <- sort(unique(membership_values[!is.na(membership_values) & nzchar(membership_values)]))
  if (!length(membership_levels)) {
    return(stats::setNames(character(), character()))
  }
  membership_factor <- factor(membership_values, levels = membership_levels)

  cluster_values <- as.character(clusters)
  cluster_levels <- sort(unique(cluster_values[!is.na(cluster_values)]))
  result <- stats::setNames(rep(NA_character_, length(membership_levels)), membership_levels)
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
  if(method == "jaccard"){
    cl.use <- as.matrix(cl.use)
    jaccard.mtrx <- cl.use
    for(i in rownames(cl.use)){
      for(j in colnames(cl.use)){
        denominator <- sum(cl.use[i,]) + sum(cl.use[,j]) - cl.use[i,j]
        jaccard.mtrx[i,j] <- if (denominator > 0) cl.use[i,j] / denominator else 0
      }
    }
    assigned <- apply(jaccard.mtrx, 2, function(x) names(x)[which.max(x)])
  } else if(method == "relative"){
    cluster.size <- rowSums(cl.gr)
    cl.use <- sweep(cl.use, 1, cluster.size, "/")
    assigned <- apply(cl.use, 2, function(x) names(x)[which.max(x)])
  } else {
    assigned <- apply(cl.use, 2, function(x) names(x)[which.max(x)])
  }

  result[names(assigned)] <- as.character(assigned)
  result
}
