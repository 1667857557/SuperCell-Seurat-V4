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

.SCOrderedMembershipLevels <- function(membership) {
  membership <- as.character(membership)
  numeric_membership <- suppressWarnings(as.integer(membership))
  if (!anyNA(numeric_membership) &&
      identical(as.character(numeric_membership), membership)) {
    return(as.character(sort(unique(numeric_membership))))
  }
  unique(membership)
}

.SCRemapSupercellGraph <- function(graph, local_map, graph_level) {
  if (!inherits(graph, "igraph")) {
    return(list(
      graph = NULL,
      aligned = FALSE,
      reason = "missing_supercell_graph"
    ))
  }
  n_vertices <- igraph::vcount(graph)
  if (n_vertices != length(local_map)) {
    return(list(
      graph = NULL,
      aligned = FALSE,
      reason = paste0(
        "graph_vertex_count_", n_vertices,
        "_differs_from_final_membership_count_", length(local_map)
      )
    ))
  }
  vertex_names <- igraph::V(graph)$name
  if (is.null(vertex_names)) {
    vertex_names <- names(local_map)
  } else {
    vertex_names <- as.character(vertex_names)
    if (!all(vertex_names %in% names(local_map))) {
      vertex_names <- names(local_map)
    }
  }
  mapped <- unname(local_map[vertex_names])
  if (anyNA(mapped) || anyDuplicated(mapped)) {
    stop(
      "Supercell graph has ambiguous membership IDs in graph group `",
      graph_level, "`.", call. = FALSE
    )
  }
  igraph::V(graph)$name <- as.character(mapped)
  list(graph = graph, aligned = TRUE, reason = NA_character_)
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
#' @return A SuperCell-style list with globally unique membership IDs,
#'   graph-group provenance, and combined graph objects when their vertices can
#'   be aligned exactly to final membership IDs. Approximate construction may
#'   return `graph.supercells = NULL` when post-hoc condition splitting creates
#'   final memberships that are absent from the contracted presample graph.
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
  if (!ncol(X)) stop("`X` must contain at least one component.", call. = FALSE)
  if (any(!is.finite(X))) {
    stop("`X` cannot contain non-finite values.", call. = FALSE)
  }
  gamma <- as.numeric(gamma)[1L]
  k.knn <- as.integer(k.knn)[1L]
  approx.N <- as.integer(approx.N)[1L]
  block.size <- as.integer(block.size)[1L]
  seed <- as.integer(seed)[1L]
  n.pc <- as.integer(n.pc)
  if (!is.finite(gamma) || gamma <= 0) {
    stop("`gamma` must be a positive finite number.", call. = FALSE)
  }
  if (is.na(k.knn) || k.knn < 1L) {
    stop("`k.knn` must be a positive integer.", call. = FALSE)
  }
  if (is.na(approx.N) || approx.N < 1L ||
      is.na(block.size) || block.size < 1L) {
    stop("`approx.N` and `block.size` must be positive integers.", call. = FALSE)
  }
  if (!length(n.pc) || anyNA(n.pc) || any(n.pc < 1L) ||
      any(n.pc > ncol(X)) || anyDuplicated(n.pc)) {
    stop("`n.pc` must contain unique valid component indices of `X`.", call. = FALSE)
  }
  if (!is.finite(seed)) stop("`seed` must be a finite integer.", call. = FALSE)
  selected_components <- n.pc
  X <- X[, selected_components, drop = FALSE]
  internal_components <- seq_len(ncol(X))

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

  group_results <- vector("list", length(group_levels))
  names(group_results) <- group_levels
  membership <- stats::setNames(integer(length(cell_ids)), cell_ids)
  membership_offset <- 0L

  for (i in seq_along(group_levels)) {
    graph_level <- group_levels[[i]]
    cells_i <- cell_ids[graph_group == graph_level]
    k_i <- min(k.knn, length(cells_i) - 1L)
    gamma_i <- min(gamma, length(cells_i))
    condition_i <- if (is.null(split_condition)) NULL else {
      unname(split_condition[cells_i])
    }
    result_i <- SCimplify_from_embedding(
      X = X[cells_i, , drop = FALSE],
      cell.annotation = NULL,
      cell.split.condition = condition_i,
      gamma = gamma_i,
      k.knn = k_i,
      n.pc = internal_components,
      do.approx = do.approx,
      approx.N = min(approx.N, length(cells_i)),
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
    local_levels <- .SCOrderedMembershipLevels(local_membership)
    local_map <- stats::setNames(
      seq.int(membership_offset + 1L,
              membership_offset + length(local_levels)),
      local_levels
    )
    membership[cells_i] <- unname(local_map[local_membership])
    membership_offset <- membership_offset + length(local_levels)
    graph_alignment <- .SCRemapSupercellGraph(
      result_i$graph.supercells, local_map, graph_level
    )
    if (!isTRUE(graph_alignment$aligned)) {
      result_i$graph.supercells.unaligned <- result_i$graph.supercells
    }
    result_i$graph.supercells <- graph_alignment$graph
    result_i$graph.supercells.aligned <- graph_alignment$aligned
    result_i$graph.supercells.unavailable.reason <- graph_alignment$reason
    result_i$graph.group <- graph_level
    result_i$input.cells <- cells_i
    result_i$global.membership.map <- local_map
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
  supercell_size <- vapply(metacell_rows, length, integer(1))
  graph_alignment <- vapply(
    group_results, `[[`, logical(1), "graph.supercells.aligned"
  )
  graph_supercells_available <- all(graph_alignment)

  result <- list(
    graph.supercells = if (graph_supercells_available) {
      .SCDisjointGraphUnion(lapply(group_results, `[[`, "graph.supercells"))
    } else {
      NULL
    },
    gamma = gamma,
    N.SC = length(metacell_rows),
    membership = membership,
    supercell_size = supercell_size,
    genes.use = NA,
    simplification.algo = igraph.clustering[[1L]],
    do.approx = do.approx,
    n.pc = selected_components,
    k.knn = k.knn,
    sc.cell.graph.group. = graph_group,
    sc.cell.split.condition. = split_condition,
    SC.cell.graph.group. = mc_graph_group,
    SC.cell.split.condition. = mc_condition,
    graph_group_levels = group_levels,
    graph_scope = "independent_by_cell.graph.group",
    condition_scope = "joint_within_graph_group_then_membership_split",
    graph_supercells_available = graph_supercells_available,
    graph_supercells_unavailable_groups = names(graph_alignment)[!graph_alignment]
  )
  if (graph_supercells_available) {
    graph_ids <- as.character(igraph::V(result$graph.supercells)$name)
    if (!setequal(graph_ids, names(metacell_rows)) || anyDuplicated(graph_ids)) {
      stop("Combined supercell graph is not aligned to global membership IDs.",
           call. = FALSE)
    }
  }
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
