.SCAlignGraphGrouping <- function(x, cell_ids, argument, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop("`", argument, "` must be supplied.", call. = FALSE)
  }
  if (length(x) != length(cell_ids)) {
    stop("`", argument, "` must contain one value per row of `X`.", call. = FALSE)
  }
  if (!is.null(names(x))) {
    if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
      stop("Named `", argument, "` values require unique, non-empty cell IDs.", call. = FALSE)
    }
    index <- match(cell_ids, names(x))
    if (anyNA(index)) {
      stop("Named `", argument, "` does not cover every row name of `X`.", call. = FALSE)
    }
    x <- x[index]
  }
  x <- trimws(as.character(x))
  if (anyNA(x) || any(!nzchar(x))) {
    stop("`", argument, "` cannot contain missing or empty values.", call. = FALSE)
  }
  stats::setNames(x, cell_ids)
}

.SCDisjointGraphUnion <- function(graphs) {
  graphs <- Filter(function(x) inherits(x, "igraph"), graphs)
  if (!length(graphs)) return(NULL)
  if (length(graphs) == 1L) return(graphs[[1L]])
  do.call(igraph::disjoint_union, graphs)
}

#' Build independent cell graphs by annotation while pooling conditions
#'
#' Builds one k-nearest-neighbour graph for each value of `cell.graph.group`.
#' All conditions within a graph group are included in the same graph. The
#' `cell.split.condition` vector is applied only after graph construction and
#' graph clustering, so returned metacells remain condition-pure without fitting
#' separate condition-specific graphs.
#'
#' @param X Low-dimensional embedding with cells in rows and components in
#'   columns. Row names must be unique cell IDs.
#' @param cell.graph.group Cell-level graph grouping vector, typically a broad
#'   cell-type annotation. A separate graph is built for each value.
#' @param cell.split.condition Optional condition vector. Conditions are pooled
#'   during graph construction and used only to split mixed memberships after
#'   clustering.
#' @param gamma Graining level passed to [SCimplify_from_embedding()] within each
#'   graph group.
#' @param k.knn Number of neighbours used within each graph group.
#' @param n.pc Components of `X` used for graph construction.
#' @param do.approx Whether to use approximate coarse-graining within each graph
#'   group.
#' @param approx.N Approximate presample size within each graph group.
#' @param block.size Mapping block size for approximate coarse-graining.
#' @param seed Base random seed. Group `i` uses `seed + i - 1`.
#' @param igraph.clustering Clustering method passed to
#'   [SCimplify_from_embedding()].
#' @param return.singlecell.NW Return the disjoint union of the independently
#'   built single-cell graphs.
#' @param return.hierarchical.structure Return a named list of group-specific
#'   hierarchical structures.
#' @param return.group.results Retain complete group-specific SuperCell results.
#' @param ... Additional arguments passed to [SCimplify_from_embedding()].
#'
#' @return A SuperCell-style list with globally unique membership IDs, combined
#'   graph objects, graph-group provenance, and optional group-specific results.
#' @export
SCimplify_by_graph_group_from_embedding <- function(
    X,
    cell.graph.group,
    cell.split.condition = NULL,
    gamma = 10,
    k.knn = 5,
    n.pc = 10,
    do.approx = FALSE,
    approx.N = 20000,
    block.size = 10000,
    seed = 12345,
    igraph.clustering = c("walktrap", "louvain"),
    return.singlecell.NW = TRUE,
    return.hierarchical.structure = TRUE,
    return.group.results = FALSE,
    ...) {
  if (is.null(dim(X)) || length(dim(X)) != 2L) {
    stop("`X` must be a two-dimensional embedding matrix.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (is.null(rownames(X)) || anyNA(rownames(X)) ||
      any(!nzchar(rownames(X))) || anyDuplicated(rownames(X))) {
    stop("`X` must have unique, non-empty cell IDs in row names.", call. = FALSE)
  }
  cell_ids <- rownames(X)
  graph_group <- .SCAlignGraphGrouping(
    cell.graph.group, cell_ids, "cell.graph.group"
  )
  split_condition <- .SCAlignGraphGrouping(
    cell.split.condition, cell_ids, "cell.split.condition", allow_null = TRUE
  )
  group_levels <- unique(unname(graph_group))
  group_sizes <- table(factor(graph_group, levels = group_levels))
  if (any(group_sizes < 2L)) {
    stop(
      "Every `cell.graph.group` level must contain at least two cells: ",
      paste(names(group_sizes)[group_sizes < 2L], collapse = ", "),
      call. = FALSE
    )
  }
  seed <- as.integer(seed)[1L]
  if (!is.finite(seed)) stop("`seed` must be a finite integer.", call. = FALSE)

  group_results <- vector("list", length(group_levels))
  names(group_results) <- group_levels
  membership <- stats::setNames(integer(length(cell_ids)), cell_ids)
  membership_offset <- 0L

  for (i in seq_along(group_levels)) {
    graph_level <- group_levels[[i]]
    cells_i <- cell_ids[graph_group == graph_level]
    k_i <- min(as.integer(k.knn), length(cells_i) - 1L)
    if (k_i < 1L) {
      stop("Graph group `", graph_level, "` cannot support kNN construction.",
           call. = FALSE)
    }
    gamma_i <- min(as.numeric(gamma), length(cells_i))
    condition_i <- if (is.null(split_condition)) NULL else {
      unname(split_condition[cells_i])
    }
    result_i <- SCimplify_from_embedding(
      X = X[cells_i, , drop = FALSE],
      cell.annotation = NULL,
      cell.split.condition = condition_i,
      gamma = gamma_i,
      k.knn = k_i,
      n.pc = n.pc,
      do.approx = do.approx,
      approx.N = min(as.integer(approx.N), length(cells_i)),
      block.size = block.size,
      seed = seed + i - 1L,
      igraph.clustering = igraph.clustering,
      return.singlecell.NW = return.singlecell.NW,
      return.hierarchical.structure = return.hierarchical.structure,
      ...
    )
    local_membership <- as.character(result_i$membership)
    if (!is.null(names(result_i$membership))) {
      index <- match(cells_i, names(result_i$membership))
      if (anyNA(index)) {
        stop("SuperCell membership cannot be aligned in graph group `",
             graph_level, "`.", call. = FALSE)
      }
      local_membership <- local_membership[index]
    }
    local_levels <- unique(local_membership)
    local_map <- stats::setNames(
      seq.int(membership_offset + 1L,
              membership_offset + length(local_levels)),
      local_levels
    )
    membership[cells_i] <- unname(local_map[local_membership])
    membership_offset <- membership_offset + length(local_levels)
    result_i$graph.group <- graph_level
    result_i$input.cells <- cells_i
    group_results[[i]] <- result_i
  }

  membership <- membership[cell_ids]
  metacell_rows <- split(cell_ids, membership)
  mc_graph_group <- vapply(metacell_rows, function(ids) {
    values <- unique(unname(graph_group[ids]))
    if (length(values) != 1L) {
      stop("A returned metacell spans multiple graph groups.", call. = FALSE)
    }
    values[[1L]]
  }, character(1))
  mc_condition <- NULL
  if (!is.null(split_condition)) {
    mc_condition <- vapply(metacell_rows, function(ids) {
      values <- unique(unname(split_condition[ids]))
      if (length(values) != 1L) {
        stop("A returned metacell spans multiple conditions.", call. = FALSE)
      }
      values[[1L]]
    }, character(1))
  }

  result <- list(
    graph.supercells = .SCDisjointGraphUnion(
      lapply(group_results, `[[`, "graph.supercells")
    ),
    gamma = gamma,
    N.SC = length(metacell_rows),
    membership = membership,
    supercell_size = as.integer(vapply(metacell_rows, length, integer(1))),
    genes.use = NA,
    simplification.algo = igraph.clustering[[1L]],
    do.approx = do.approx,
    n.pc = n.pc,
    k.knn = k.knn,
    sc.cell.graph.group. = graph_group,
    sc.cell.split.condition. = split_condition,
    SC.cell.graph.group. = mc_graph_group,
    SC.cell.split.condition. = mc_condition,
    graph_group_levels = group_levels,
    graph_scope = "independent_by_cell.graph.group",
    condition_scope = "joint_within_graph_group_then_membership_split"
  )
  if (isTRUE(return.singlecell.NW)) {
    result$graph.singlecell <- .SCDisjointGraphUnion(
      lapply(group_results, `[[`, "graph.singlecell")
    )
  }
  if (isTRUE(return.hierarchical.structure)) {
    result$h_membership <- lapply(group_results, `[[`, "h_membership")
  }
  if (isTRUE(return.group.results)) result$group.results <- group_results
  result
}
