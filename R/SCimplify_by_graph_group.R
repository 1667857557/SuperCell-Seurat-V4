.SCAlignGraphGrouping <- function(x, cell_ids, argument, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop("`", argument, "` must be supplied.", call. = FALSE)
  }
  if (length(x) != length(cell_ids)) {
    stop("`", argument, "` must contain one value per input cell.",
         call. = FALSE)
  }
  if (!is.null(names(x))) {
    if (anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
      stop("Named `", argument,
           "` values require unique, non-empty cell IDs.", call. = FALSE)
    }
    index <- match(cell_ids, names(x))
    if (anyNA(index)) {
      stop("Named `", argument,
           "` does not cover every input cell.", call. = FALSE)
    }
    x <- x[index]
  }
  x <- trimws(as.character(x))
  if (anyNA(x) || any(!nzchar(x))) {
    stop("`", argument, "` cannot contain missing or empty values.",
         call. = FALSE)
  }
  stats::setNames(x, cell_ids)
}

.SCStableMetacellIds <- function(keys) {
  levels <- unique(as.character(keys))
  width <- max(3L, nchar(length(levels)))
  ids <- paste0("MC", sprintf(paste0("%0", width, "d"), seq_along(levels)))
  stats::setNames(ids, levels)
}

.SCBuildGraphMembership <- function(
    seurat, k.knn, kith, kernel, gamma, graph.name,
    assay, reduction, dims, seed, verbose) {
  set.seed(as.integer(seed))
  if (length(assay) == 1L) {
    if (is.null(graph.name)) graph.name <- "nn"
    graph <- ComputeUnimodalKnn(
      seurat = seurat,
      k.knn = k.knn,
      kith = kith,
      kernel = kernel,
      graph.name = graph.name,
      assay = assay,
      reduction = reduction,
      dims = dims,
      verbose = verbose
    )
  } else {
    if (is.null(graph.name)) graph.name <- "knn"
    graph <- ComputeMultimodalKnn(
      seurat = seurat,
      k.knn = k.knn,
      kith = kith,
      kernel = kernel,
      graph.name = graph.name,
      assay = assay,
      reduction = reduction,
      dims = dims,
      verbose = verbose
    )
  }
  walktrap <- igraph::cluster_walktrap(graph)
  n_target <- .SCResolveTargetMetacells(n_cells = ncol(seurat), gamma = gamma)
  membership <- igraph::cut_at(walktrap, no = n_target)
  names(membership) <- colnames(seurat)
  list(graph = graph, walktrap = walktrap, membership = membership)
}

.SCSharedGraphMatrix <- function(graph, cell_ids) {
  if (!inherits(graph, "igraph")) {
    stop("Small-metacell repair requires the original igraph WNN.",
         call. = FALSE)
  }
  vertex_ids <- as.character(igraph::V(graph)$name)
  if (!length(vertex_ids) || anyNA(vertex_ids) || any(!nzchar(vertex_ids)) ||
      anyDuplicated(vertex_ids) || !setequal(vertex_ids, cell_ids)) {
    stop("Original WNN vertex IDs do not match the graph-group cells.",
         call. = FALSE)
  }
  weight <- igraph::edge_attr(graph, "weight")
  if (!is.null(weight) && (any(!is.finite(weight)) || any(weight < 0))) {
    stop("Original WNN weights must be finite and non-negative.",
         call. = FALSE)
  }
  W <- igraph::as_adjacency_matrix(
    graph,
    attr = if (is.null(weight)) NULL else "weight",
    sparse = TRUE
  )
  W <- W[cell_ids, cell_ids, drop = FALSE]
  W <- (W + Matrix::t(W)) / 2
  if (length(W@x) && (any(!is.finite(W@x)) || any(W@x < 0))) {
    stop("Symmetrized WNN weights must be finite and non-negative.",
         call. = FALSE)
  }
  diagonal <- Matrix::diag(W)
  if (any(diagonal != 0)) {
    W <- W - Matrix::Diagonal(x = diagonal)
  }
  Matrix::drop0(W)
}

.SCMetacellAffinity <- function(W, degree, membership, source, target) {
  source_cells <- names(membership)[membership == source]
  target_cells <- names(membership)[membership == target]
  if (!length(source_cells) || !length(target_cells)) return(NA_real_)
  cross_weight <- sum(W[source_cells, target_cells, drop = FALSE])
  if (!is.finite(cross_weight) || cross_weight <= 0) return(NA_real_)
  source_volume <- sum(degree[source_cells])
  target_volume <- sum(degree[target_cells])
  if (!is.finite(source_volume) || !is.finite(target_volume) ||
      source_volume <= 0 || target_volume <= 0) {
    return(NA_real_)
  }
  cross_weight / sqrt(source_volume * target_volume)
}

.SCRepairSplitMembership <- function(
    graph, provisional_membership, stratum,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L,
    min_merge_affinity = NULL,
    unresolved_small_policy = c("error", "keep")) {
  unresolved_small_policy <- match.arg(unresolved_small_policy)
  cell_ids <- names(provisional_membership)
  if (!length(cell_ids) || anyNA(cell_ids) || any(!nzchar(cell_ids)) ||
      anyDuplicated(cell_ids)) {
    stop("Provisional membership requires unique named cell IDs.",
         call. = FALSE)
  }
  membership <- trimws(as.character(provisional_membership))
  names(membership) <- cell_ids
  if (anyNA(membership) || any(!nzchar(membership))) {
    stop("Provisional metacell IDs cannot be missing or empty.",
         call. = FALSE)
  }
  stratum <- .SCAlignGraphGrouping(stratum, cell_ids, "repair stratum")
  min_metacell_size <- as.integer(min_metacell_size)[1L]
  min_metacells_per_stratum <- as.integer(min_metacells_per_stratum)[1L]
  if (is.na(min_metacell_size) || min_metacell_size < 1L ||
      is.na(min_metacells_per_stratum) || min_metacells_per_stratum < 1L) {
    stop("Metacell repair size/count controls must be positive integers.",
         call. = FALSE)
  }
  if (min_metacell_size == 1L) {
    return(list(
      membership = membership,
      merge_diagnostics = data.frame(),
      unresolved = data.frame(),
      symmetrization = "(W+t(W))/2",
      algorithm = "shared_wnn_condition_split_affinity_v1"
    ))
  }
  if (is.null(min_merge_affinity)) {
    stop(
      "`min_merge_affinity` must be supplied explicitly when ",
      "`min_metacell_size > 1`.", call. = FALSE
    )
  }
  min_merge_affinity <- suppressWarnings(as.numeric(min_merge_affinity)[1L])
  if (!is.finite(min_merge_affinity) || min_merge_affinity < 0 ||
      min_merge_affinity > 1) {
    stop("`min_merge_affinity` must be one finite value in [0, 1].",
         call. = FALSE)
  }

  stratum_sizes <- table(stratum)
  required_cells <- min_metacell_size * min_metacells_per_stratum
  impossible <- names(stratum_sizes)[stratum_sizes < required_cells]
  if (length(impossible)) {
    stop(
      "Metacell repair is infeasible before graph constraints for strata: ",
      paste(impossible, collapse = ", "),
      "; each needs at least ", required_cells, " cells.",
      call. = FALSE
    )
  }

  W <- .SCSharedGraphMatrix(graph, cell_ids)
  degree <- Matrix::rowSums(W)
  names(degree) <- cell_ids
  merge_rows <- list()

  repeat {
    group_cells <- split(cell_ids, membership)
    group_size <- vapply(group_cells, length, integer(1))
    small <- names(group_size)[group_size < min_metacell_size]
    if (!length(small)) break
    small <- small[order(group_size[small], small)]
    merged <- FALSE

    for (source in small) {
      source_cells <- cell_ids[membership == source]
      if (!length(source_cells) || length(source_cells) >= min_metacell_size) next
      source_stratum <- unique(unname(stratum[source_cells]))
      if (length(source_stratum) != 1L) {
        stop("A provisional metacell spans multiple repair strata.",
             call. = FALSE)
      }
      stratum_cells <- cell_ids[unname(stratum) == source_stratum]
      current_ids <- sort(unique(membership[stratum_cells]))
      if (length(current_ids) <= min_metacells_per_stratum) next
      candidates <- setdiff(current_ids, source)
      if (!length(candidates)) next
      affinity <- vapply(
        candidates,
        function(target) .SCMetacellAffinity(
          W, degree, membership, source, target
        ),
        numeric(1)
      )
      valid <- is.finite(affinity) & affinity >= min_merge_affinity
      if (!any(valid)) next
      candidates <- candidates[valid]
      affinity <- affinity[valid]
      best <- max(affinity)
      tied <- sort(candidates[abs(affinity - best) <= 1e-12])
      target <- tied[[1L]]
      source_size <- sum(membership == source)
      target_size <- sum(membership == target)
      membership[membership == source] <- target
      merge_rows[[length(merge_rows) + 1L]] <- data.frame(
        source_metacell_id = source,
        target_metacell_id = target,
        stratum = source_stratum,
        affinity = best,
        source_size = source_size,
        target_size_before = target_size,
        target_size_after = source_size + target_size,
        stringsAsFactors = FALSE
      )
      merged <- TRUE
      break
    }
    if (!merged) break
  }

  group_cells <- split(cell_ids, membership)
  group_size <- vapply(group_cells, length, integer(1))
  small <- names(group_size)[group_size < min_metacell_size]
  unresolved <- if (length(small)) {
    do.call(rbind, lapply(sort(small), function(id) {
      cells <- group_cells[[id]]
      data.frame(
        metacell_id = id,
        stratum = unique(unname(stratum[cells]))[[1L]],
        n_cells = length(cells),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      metacell_id = character(), stratum = character(), n_cells = integer(),
      stringsAsFactors = FALSE
    )
  }
  if (nrow(unresolved) && identical(unresolved_small_policy, "error")) {
    stop(
      "Small metacell repair left unresolved groups below `min_metacell_size`: ",
      paste(
        paste0(unresolved$metacell_id, "[", unresolved$n_cells, "]"),
        collapse = ", "
      ),
      ". No same-stratum original-WNN neighbor met the affinity/count constraints.",
      call. = FALSE
    )
  }

  final_counts <- vapply(
    split(cell_ids, stratum),
    function(cells) length(unique(membership[cells])),
    integer(1)
  )
  if (any(final_counts < min_metacells_per_stratum)) {
    stop("Small-metacell repair violated `min_metacells_per_stratum`.",
         call. = FALSE)
  }

  list(
    membership = membership,
    merge_diagnostics = if (length(merge_rows)) {
      do.call(rbind, merge_rows)
    } else {
      data.frame(
        source_metacell_id = character(), target_metacell_id = character(),
        stratum = character(), affinity = numeric(), source_size = integer(),
        target_size_before = integer(), target_size_after = integer(),
        stringsAsFactors = FALSE
      )
    },
    unresolved = unresolved,
    symmetrization = "(W+t(W))/2",
    algorithm = "shared_wnn_condition_split_affinity_v1"
  )
}

#' Build cell-type-scoped multimodal metacells
#'
#' Builds one independent multimodal WNN graph for each value of
#' `cell.graph.group`. All conditions inside a graph group are included in the
#' same WNN construction and Walktrap clustering. `cell.split.condition` is
#' applied only after clustering. Optionally, condition-split metacells smaller
#' than `min_metacell_size` are repaired using affinity from that exact original
#' shared WNN, never by rebuilding a condition-specific graph.
#'
#' @param seurat A preprocessed Seurat object containing the requested assays
#'   and reductions.
#' @param cell.graph.group Cell-level graph grouping vector, typically a broad
#'   cell-type annotation. A separate WNN graph is built for each value.
#' @param cell.split.condition Optional condition vector. Conditions are pooled
#'   during WNN construction and used only to split memberships after graph
#'   clustering.
#' @param k.knn Number of neighbours in the WNN graph.
#' @param kith Optional neighbourhood rank passed to graph construction.
#' @param kernel Whether to use the SuperCell kernel-weighted graph.
#' @param gamma Target average number of cells per metacell.
#' @param graph.name Optional graph name used during graph construction.
#' @param assay One or two assays used for graph construction.
#' @param reduction Corresponding dimensional reductions.
#' @param dims Corresponding dimension vectors.
#' @param seed Base random seed. Graph group `i` uses `seed + i - 1`.
#' @param min_metacell_size Minimum final metacell size. Values above one enable
#'   post-condition-split repair.
#' @param min_metacells_per_stratum Minimum number of final metacells retained
#'   in each condition/graph-group stratum.
#' @param min_merge_affinity Explicit normalized original-WNN affinity threshold
#'   for repair. Required when `min_metacell_size > 1`.
#' @param unresolved_small_policy Either `error` or `keep` for small metacells
#'   with no legal sufficiently affine merge candidate.
#' @param return.group.results Retain compact group-specific clustering results.
#' @param verbose Forward progress messages to SuperCell graph construction.
#'
#' @return A list containing globally unique condition-pure membership IDs,
#'   original parent Walktrap-cluster memberships, repaired membership, sizes,
#'   hierarchies and repair diagnostics. Original cell-level parent IDs are
#'   retained even if one final repaired metacell contains multiple parents.
#' @export
SCimplify_by_graph_group <- function(
    seurat,
    cell.graph.group,
    cell.split.condition = NULL,
    k.knn = 30,
    kith = NULL,
    kernel = TRUE,
    gamma = 30,
    graph.name = NULL,
    assay = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims = list(1:30, 2:30),
    seed = 12345L,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L,
    min_merge_affinity = NULL,
    unresolved_small_policy = c("error", "keep"),
    return.group.results = FALSE,
    verbose = FALSE) {
  if (!inherits(seurat, "Seurat")) {
    stop("`seurat` must inherit from Seurat.", call. = FALSE)
  }
  cell_ids <- as.character(colnames(seurat))
  if (!length(cell_ids) || anyNA(cell_ids) || any(!nzchar(cell_ids)) ||
      anyDuplicated(cell_ids)) {
    stop("The Seurat object must contain unique, non-empty cell IDs.",
         call. = FALSE)
  }
  if (length(assay) < 1L || length(assay) > 2L ||
      length(reduction) != length(assay) || length(dims) != length(assay)) {
    stop("Graph-group construction requires one or two aligned assays, reductions, and dimension vectors.",
         call. = FALSE)
  }
  assay <- as.character(assay)
  reduction <- as.list(reduction)
  dims <- lapply(dims, as.integer)
  if (anyNA(assay) || any(!nzchar(assay)) || anyDuplicated(assay) ||
      any(!vapply(dims, length, integer(1))) ||
      any(vapply(dims, function(x) anyNA(x) || any(x < 1L), logical(1)))) {
    stop("Assay, reduction, and dimension specifications are invalid.",
         call. = FALSE)
  }
  missing_assays <- setdiff(assay, names(seurat@assays))
  missing_reductions <- setdiff(unlist(reduction, use.names = FALSE),
                                names(seurat@reductions))
  if (length(missing_assays)) {
    stop("Missing assay(s): ", paste(missing_assays, collapse = ", "), ".",
         call. = FALSE)
  }
  if (length(missing_reductions)) {
    stop("Missing reduction(s): ",
         paste(missing_reductions, collapse = ", "), ".", call. = FALSE)
  }
  k.knn <- as.integer(k.knn)[1L]
  gamma <- suppressWarnings(as.numeric(gamma)[1L])
  seed <- as.integer(seed)[1L]
  min_metacell_size <- as.integer(min_metacell_size)[1L]
  min_metacells_per_stratum <- as.integer(min_metacells_per_stratum)[1L]
  unresolved_small_policy <- match.arg(unresolved_small_policy)
  if (is.na(k.knn) || k.knn < 1L) {
    stop("`k.knn` must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(gamma) || gamma <= 0) {
    stop("`gamma` must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(seed)) {
    stop("`seed` must be a finite integer.", call. = FALSE)
  }
  if (is.na(min_metacell_size) || min_metacell_size < 1L ||
      is.na(min_metacells_per_stratum) || min_metacells_per_stratum < 1L) {
    stop("Metacell size/count controls must be positive integers.",
         call. = FALSE)
  }
  if (min_metacell_size > 1L && is.null(min_merge_affinity)) {
    stop("`min_merge_affinity` is required when `min_metacell_size > 1`.",
         call. = FALSE)
  }
  if (!is.null(min_merge_affinity)) {
    min_merge_affinity <- suppressWarnings(as.numeric(min_merge_affinity)[1L])
    if (!is.finite(min_merge_affinity) || min_merge_affinity < 0 ||
        min_merge_affinity > 1) {
      stop("`min_merge_affinity` must be one finite value in [0, 1].",
           call. = FALSE)
    }
  }
  if (!is.logical(kernel) || length(kernel) != 1L || is.na(kernel) ||
      !is.logical(return.group.results) || length(return.group.results) != 1L ||
      is.na(return.group.results) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`kernel`, `return.group.results`, and `verbose` must be TRUE or FALSE.",
         call. = FALSE)
  }

  graph_group <- .SCAlignGraphGrouping(
    cell.graph.group, cell_ids, "cell.graph.group"
  )
  split_condition <- .SCAlignGraphGrouping(
    cell.split.condition, cell_ids, "cell.split.condition", allow_null = TRUE
  )
  repair_stratum <- if (is.null(split_condition)) {
    graph_group
  } else {
    stats::setNames(
      paste(unname(graph_group), unname(split_condition), sep = "\001"),
      cell_ids
    )
  }
  required_cells <- min_metacell_size * min_metacells_per_stratum
  stratum_sizes <- table(repair_stratum)
  impossible <- names(stratum_sizes)[stratum_sizes < required_cells]
  if (length(impossible)) {
    stop(
      "Condition/graph-group strata cannot satisfy requested metacell constraints: ",
      paste(impossible, collapse = ", "),
      "; each needs at least ", required_cells, " cells.", call. = FALSE
    )
  }

  group_levels <- unique(unname(graph_group))
  group_sizes <- table(factor(graph_group, levels = group_levels))
  if (any(group_sizes < 2L)) {
    stop(
      "Every `cell.graph.group` level must contain at least two cells: ",
      paste(names(group_sizes)[group_sizes < 2L], collapse = ", "),
      call. = FALSE
    )
  }

  old_random_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                                   inherits = FALSE)
  old_random_seed <- if (old_random_seed_exists) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (old_random_seed_exists) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  parent_membership <- stats::setNames(character(length(cell_ids)), cell_ids)
  final_key <- stats::setNames(character(length(cell_ids)), cell_ids)
  group_results <- vector("list", length(group_levels))
  names(group_results) <- group_levels
  hierarchies <- vector("list", length(group_levels))
  names(hierarchies) <- group_levels
  repair_rows <- list()
  unresolved_rows <- list()

  for (i in seq_along(group_levels)) {
    graph_level <- group_levels[[i]]
    cells_i <- cell_ids[graph_group == graph_level]
    object_i <- subset(seurat, cells = cells_i)
    graph_name_i <- if (is.null(graph.name)) NULL else {
      paste0(as.character(graph.name)[1L], "_", i)
    }
    built <- .SCBuildGraphMembership(
      seurat = object_i,
      k.knn = min(k.knn, length(cells_i) - 1L),
      kith = kith,
      kernel = kernel,
      gamma = gamma,
      graph.name = graph_name_i,
      assay = assay,
      reduction = reduction,
      dims = dims,
      seed = seed + i - 1L,
      verbose = verbose
    )
    local <- as.character(built$membership[cells_i])
    if (length(local) != length(cells_i) || anyNA(local) ||
        any(!nzchar(local))) {
      stop("SuperCell returned invalid memberships in graph group `",
           graph_level, "`.", call. = FALSE)
    }
    parent <- paste0("G", i, "::", local)
    names(parent) <- cells_i
    parent_membership[cells_i] <- parent
    provisional <- if (is.null(split_condition)) {
      parent
    } else {
      stats::setNames(
        paste(parent, split_condition[cells_i], sep = "\001"),
        cells_i
      )
    }
    repaired <- .SCRepairSplitMembership(
      graph = built$graph,
      provisional_membership = provisional,
      stratum = repair_stratum[cells_i],
      min_metacell_size = min_metacell_size,
      min_metacells_per_stratum = min_metacells_per_stratum,
      min_merge_affinity = min_merge_affinity,
      unresolved_small_policy = unresolved_small_policy
    )
    final_key[cells_i] <- repaired$membership[cells_i]
    hierarchies[[i]] <- built$walktrap
    if (nrow(repaired$merge_diagnostics)) {
      tab <- repaired$merge_diagnostics
      tab$graph_group <- graph_level
      repair_rows[[length(repair_rows) + 1L]] <- tab
    }
    if (nrow(repaired$unresolved)) {
      tab <- repaired$unresolved
      tab$graph_group <- graph_level
      unresolved_rows[[length(unresolved_rows) + 1L]] <- tab
    }
    if (isTRUE(return.group.results)) {
      group_results[[i]] <- list(
        membership = stats::setNames(local, cells_i),
        repaired_membership = repaired$membership,
        supercell_size = as.integer(table(local)),
        h_membership = built$walktrap,
        graph.group = graph_level,
        input.cells = cells_i,
        repair = repaired[c("merge_diagnostics", "unresolved",
                            "symmetrization", "algorithm")]
      )
    }
    built$graph <- NULL
    object_i <- NULL
    invisible(gc(verbose = FALSE, full = TRUE))
  }

  if (any(!nzchar(parent_membership)) || any(!nzchar(final_key))) {
    stop("SuperCell did not assign every input cell.", call. = FALSE)
  }
  id_map <- .SCStableMetacellIds(final_key[cell_ids])
  membership <- stats::setNames(unname(id_map[final_key[cell_ids]]), cell_ids)
  membership_groups <- split(cell_ids, membership)

  metacell_graph_group <- vapply(membership_groups, function(ids) {
    values <- unique(unname(graph_group[ids]))
    if (length(values) != 1L) {
      stop("A returned metacell spans multiple graph groups.", call. = FALSE)
    }
    values[[1L]]
  }, character(1))
  metacell_condition <- NULL
  if (!is.null(split_condition)) {
    metacell_condition <- vapply(membership_groups, function(ids) {
      values <- unique(unname(split_condition[ids]))
      if (length(values) != 1L) {
        stop("A returned metacell spans multiple conditions.", call. = FALSE)
      }
      values[[1L]]
    }, character(1))
  }
  final_stratum_count <- table(vapply(membership_groups, function(ids) {
    unique(unname(repair_stratum[ids]))[[1L]]
  }, character(1)))
  if (any(final_stratum_count < min_metacells_per_stratum)) {
    stop("Final membership violates `min_metacells_per_stratum`.",
         call. = FALSE)
  }

  membership_table <- data.frame(
    cell_id = cell_ids,
    metacell_id = unname(membership[cell_ids]),
    parent_metacell_id = unname(parent_membership[cell_ids]),
    graph_group = unname(graph_group[cell_ids]),
    stringsAsFactors = FALSE
  )
  if (!is.null(split_condition)) {
    membership_table$condition <- unname(split_condition[cell_ids])
  }
  supercell_size <- vapply(membership_groups, length, integer(1))
  repair_diagnostics <- if (length(repair_rows)) {
    do.call(rbind, repair_rows)
  } else data.frame()
  unresolved_small <- if (length(unresolved_rows)) {
    do.call(rbind, unresolved_rows)
  } else data.frame()

  result <- list(
    membership = membership,
    membership_table = membership_table,
    parent_membership = parent_membership,
    supercell_size = supercell_size,
    N.SC = length(membership_groups),
    gamma = gamma,
    k.knn = k.knn,
    graph_group_levels = group_levels,
    sc.cell.graph.group. = graph_group,
    sc.cell.split.condition. = split_condition,
    SC.cell.graph.group. = metacell_graph_group,
    SC.cell.split.condition. = metacell_condition,
    h_membership = hierarchies,
    repair_diagnostics = repair_diagnostics,
    unresolved_small_metacells = unresolved_small,
    repair_contract = list(
      algorithm = "shared_wnn_condition_split_affinity_v1",
      min_metacell_size = min_metacell_size,
      min_metacells_per_stratum = min_metacells_per_stratum,
      min_merge_affinity = min_merge_affinity,
      unresolved_small_policy = unresolved_small_policy,
      affinity = "sum(W_MN)/sqrt(vol(M)*vol(N))",
      symmetrization = "(W+t(W))/2",
      candidate_scope = "same_condition_same_graph_group_original_WNN",
      parent_provenance = "cell_level_original_walktrap_parent_preserved"
    ),
    graph_scope = "independent_WNN_by_cell.graph.group",
    condition_scope = "joint_within_graph_group_then_membership_split",
    graph_method = "SCimplify_for_Seurat_equivalent_multimodal_WNN_walktrap",
    modality_weighting = "adaptive_WNN_within_graph_group"
  )
  if (isTRUE(return.group.results)) result$group.results <- group_results
  result
}
