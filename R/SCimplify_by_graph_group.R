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

#' Build cell-type-scoped multimodal metacells
#'
#' Builds one independent multimodal WNN graph and one Walktrap hierarchy for
#' each value of `cell.graph.group`. All conditions inside a graph group are
#' included in the same graph and hierarchy. The legacy policy cuts the shared
#' hierarchy at the global `gamma` target and then splits memberships by
#' condition. The hierarchy-constrained policy instead selects a
#' condition-specific gamma-resolution cut of that same hierarchy, then repairs
#' only condition-specific metacells below `min.metacell.size` by following the
#' existing Walktrap merge hierarchy upward to the nearest same-condition
#' sibling subtree. Already-valid metacells never merge with one another. The
#' policy never rebuilds a condition-specific graph or Walktrap hierarchy.
#'
#' @param seurat A preprocessed Seurat object containing the requested assays
#'   and reductions.
#' @param cell.graph.group Cell-level graph grouping vector, typically a broad
#'   cell-type annotation. A separate WNN graph is built for each value.
#' @param cell.split.condition Optional condition vector. Conditions are pooled
#'   during WNN and Walktrap construction.
#' @param k.knn Number of neighbours passed to [SCimplify_for_Seurat()].
#' @param kith Optional neighbourhood rank passed to
#'   [SCimplify_for_Seurat()].
#' @param kernel Whether to use the SuperCell kernel-weighted graph.
#' @param gamma Target average number of cells per initial condition-specific
#'   metacell before local minimum-size repair in hierarchy-constrained mode.
#' @param condition.partition Final condition partition policy. The default
#'   `"legacy_post_split"` preserves historical behavior. Use
#'   `"hierarchy_constrained"` for a condition-specific gamma cut followed by
#'   local minimum-size repair on the one shared Walktrap hierarchy.
#' @param min.metacell.size Minimum final metacell size in hierarchy-constrained
#'   mode. This is a local repair threshold and does not coarsen the initial
#'   gamma-resolution cut globally.
#' @param min.metacells.per.condition Hard minimum final metacells per condition
#'   in hierarchy-constrained mode.
#' @param graph.name Optional graph name passed to
#'   [SCimplify_for_Seurat()].
#' @param assay Two assays used for multimodal WNN construction, typically RNA
#'   and ATAC.
#' @param reduction Two corresponding dimensional reductions.
#' @param dims Two corresponding dimension vectors.
#' @param seed Base random seed. Graph group `i` uses `seed + i - 1`.
#' @param return.group.results Retain the group-specific
#'   [SCimplify_for_Seurat()] results.
#' @param verbose Forward progress messages to SuperCell graph construction.
#'
#' @return A list containing globally unique condition-pure membership IDs,
#'   initial shared-hierarchy provenance, local repair provenance, a canonical
#'   cell-level membership table, metacell sizes, partition diagnostics,
#'   group-specific hierarchies, and graph-construction provenance.
#' @export
SCimplify_by_graph_group <- function(
    seurat,
    cell.graph.group,
    cell.split.condition = NULL,
    k.knn = 30,
    kith = NULL,
    kernel = TRUE,
    gamma = 30,
    condition.partition = c("legacy_post_split", "hierarchy_constrained"),
    min.metacell.size = 1L,
    min.metacells.per.condition = 1L,
    graph.name = NULL,
    assay = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims = list(1:30, 2:30),
    seed = 12345L,
    return.group.results = FALSE,
    verbose = FALSE) {
  if (!inherits(seurat, "Seurat")) {
    stop("`seurat` must inherit from Seurat.", call. = FALSE)
  }
  condition.partition <- match.arg(condition.partition)
  cell_ids <- as.character(colnames(seurat))
  if (!length(cell_ids) || anyNA(cell_ids) || any(!nzchar(cell_ids)) ||
      anyDuplicated(cell_ids)) {
    stop("The Seurat object must contain unique, non-empty cell IDs.",
         call. = FALSE)
  }
  if (length(assay) != 2L || length(reduction) != 2L || length(dims) != 2L) {
    stop("Multimodal graph-group construction requires exactly two assays, reductions, and dimension vectors.",
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
  if (!is.numeric(min.metacell.size) || length(min.metacell.size) != 1L ||
      is.na(min.metacell.size) || !is.finite(min.metacell.size) ||
      min.metacell.size < 1 || min.metacell.size != floor(min.metacell.size) ||
      min.metacell.size > .Machine$integer.max) {
    stop("`min.metacell.size` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(min.metacells.per.condition) ||
      length(min.metacells.per.condition) != 1L ||
      is.na(min.metacells.per.condition) ||
      !is.finite(min.metacells.per.condition) ||
      min.metacells.per.condition < 1 ||
      min.metacells.per.condition != floor(min.metacells.per.condition) ||
      min.metacells.per.condition > .Machine$integer.max) {
    stop("`min.metacells.per.condition` must be a positive integer.",
         call. = FALSE)
  }
  k.knn <- as.integer(k.knn)[1L]
  gamma <- suppressWarnings(as.numeric(gamma)[1L])
  seed <- as.integer(seed)[1L]
  min.metacell.size <- as.integer(min.metacell.size)
  min.metacells.per.condition <- as.integer(min.metacells.per.condition)
  if (is.na(k.knn) || k.knn < 1L) {
    stop("`k.knn` must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(gamma) || gamma <= 0) {
    stop("`gamma` must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(seed)) {
    stop("`seed` must be a finite integer.", call. = FALSE)
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
  group_levels <- unique(unname(graph_group))
  group_sizes <- table(factor(graph_group, levels = group_levels))
  if (any(group_sizes < 2L)) {
    stop(
      "Every `cell.graph.group` level must contain at least two cells: ",
      paste(names(group_sizes)[group_sizes < 2L], collapse = ", "),
      call. = FALSE
    )
  }

  parent_membership <- stats::setNames(character(length(cell_ids)), cell_ids)
  final_key <- stats::setNames(character(length(cell_ids)), cell_ids)
  shared_community_id <- stats::setNames(character(length(cell_ids)), cell_ids)
  shared_cut_k <- stats::setNames(rep(NA_integer_, length(cell_ids)), cell_ids)
  partition_policy_cell <- stats::setNames(character(length(cell_ids)), cell_ids)
  group_results <- vector("list", length(group_levels))
  names(group_results) <- group_levels
  hierarchies <- vector("list", length(group_levels))
  names(hierarchies) <- group_levels
  partition_diagnostics <- list()
  partition_repairs <- list()

  for (i in seq_along(group_levels)) {
    graph_level <- group_levels[[i]]
    cells_i <- cell_ids[graph_group == graph_level]
    object_i <- subset(seurat, cells = cells_i)
    graph_name_i <- if (is.null(graph.name)) NULL else {
      paste0(as.character(graph.name)[1L], "_", i)
    }
    result_i <- SCimplify_for_Seurat(
      seurat = object_i,
      seurat.mc = NULL,
      k.knn = min(k.knn, length(cells_i) - 1L),
      kith = kith,
      kernel = kernel,
      gamma = gamma,
      graph.name = graph_name_i,
      assay = assay,
      reduction = reduction,
      dims = dims,
      membership = NULL,
      metacellNormalization = FALSE,
      avg.in.data = FALSE,
      fragmentFiles = NULL,
      seed = seed + i - 1L,
      prefixMC = "",
      label = NULL,
      return.seurat = FALSE,
      verbose = verbose
    )
    local <- as.character(result_i$membership)
    if (length(local) != length(cells_i)) {
      stop("SuperCell membership length differs from graph group `",
           graph_level, "`.", call. = FALSE)
    }
    if (!is.null(names(result_i$membership))) {
      index <- match(cells_i, names(result_i$membership))
      if (anyNA(index)) {
        stop("SuperCell membership does not cover graph group `",
             graph_level, "`.", call. = FALSE)
      }
      local <- local[index]
    }
    if (anyNA(local) || any(!nzchar(local))) {
      stop("SuperCell returned missing memberships in graph group `",
           graph_level, "`.", call. = FALSE)
    }
    names(local) <- cells_i
    hierarchies[[i]] <- result_i$h_membership

    if (identical(condition.partition, "hierarchy_constrained") &&
        !is.null(split_condition)) {
      if (is.null(result_i$h_membership)) {
        stop("Hierarchy-constrained partitioning requires the shared Walktrap hierarchy.",
             call. = FALSE)
      }
      partition_i <- .SCPartitionSharedHierarchyByCondition(
        hierarchy = result_i$h_membership,
        cell_ids = cells_i,
        condition = split_condition[cells_i],
        gamma = gamma,
        min.metacell.size = min.metacell.size,
        min.metacells.per.condition = min.metacells.per.condition
      )
      parent_membership[cells_i] <- paste0(
        "G", i, "::", unname(partition_i$parent_key[cells_i])
      )
      final_key[cells_i] <- paste0(
        "G", i, "::", unname(partition_i$final_key[cells_i])
      )
      shared_community_id[cells_i] <-
        unname(partition_i$shared_community_id[cells_i])
      shared_cut_k[cells_i] <- unname(partition_i$shared_cut_k[cells_i])
      partition_policy_cell[cells_i] <- "hierarchy_constrained"
      diagnostic <- partition_i$diagnostics
      diagnostic$graph_group <- graph_level
      partition_diagnostics[[length(partition_diagnostics) + 1L]] <- diagnostic
      if (is.data.frame(partition_i$repairs) && nrow(partition_i$repairs)) {
        repairs_i <- partition_i$repairs
        repairs_i$graph_group <- graph_level
        partition_repairs[[length(partition_repairs) + 1L]] <- repairs_i
      }
    } else {
      parent_membership[cells_i] <- paste0("G", i, "::", local[cells_i])
      final_key[cells_i] <- if (is.null(split_condition)) {
        parent_membership[cells_i]
      } else {
        paste(parent_membership[cells_i], split_condition[cells_i], sep = "\001")
      }
      shared_community_id[cells_i] <- unname(local[cells_i])
      if (is.null(split_condition)) {
        shared_cut_k[cells_i] <- length(unique(local))
        partition_policy_cell[cells_i] <- "global_hierarchy_cut_no_split"
      } else {
        partition_policy_cell[cells_i] <- "legacy_post_split"
      }
    }

    if (isTRUE(return.group.results)) {
      result_i$graph.group <- graph_level
      result_i$input.cells <- cells_i
      group_results[[i]] <- result_i
    }
  }

  if (any(!nzchar(parent_membership)) || any(!nzchar(final_key))) {
    stop("SuperCell did not assign every input cell.", call. = FALSE)
  }
  id_map <- if (identical(condition.partition, "hierarchy_constrained") &&
                !is.null(split_condition)) {
    .SCStableMetacellIds(sort(unique(unname(final_key))))
  } else {
    .SCStableMetacellIds(final_key)
  }
  membership <- stats::setNames(unname(id_map[final_key]), cell_ids)
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

  membership_table <- data.frame(
    cell_id = cell_ids,
    metacell_id = unname(membership[cell_ids]),
    parent_metacell_id = unname(parent_membership[cell_ids]),
    graph_group = unname(graph_group[cell_ids]),
    partition_policy = unname(partition_policy_cell[cell_ids]),
    shared_cut_k = as.integer(shared_cut_k[cell_ids]),
    shared_community_id = unname(shared_community_id[cell_ids]),
    stringsAsFactors = FALSE
  )
  if (!is.null(split_condition)) {
    membership_table$condition <- unname(split_condition[cell_ids])
  }
  supercell_size <- vapply(membership_groups, length, integer(1))
  if (identical(condition.partition, "hierarchy_constrained") &&
      !is.null(split_condition) && any(supercell_size < min.metacell.size)) {
    stop("Hierarchy-constrained partition returned a metacell below the hard minimum size.",
         call. = FALSE)
  }

  partition_policy <- if (is.null(split_condition)) {
    "global_hierarchy_cut_no_split"
  } else {
    condition.partition
  }
  partition_schema_version <- switch(
    partition_policy,
    hierarchy_constrained = "shared_walktrap_condition_local_repair_v2",
    legacy_post_split = "legacy_post_split_v1",
    global_hierarchy_cut_no_split = "global_native_gamma_v1"
  )
  condition_scope <- switch(
    partition_policy,
    hierarchy_constrained =
      "joint_within_graph_group_shared_hierarchy_gamma_cut_local_tree_repair",
    legacy_post_split = "joint_within_graph_group_then_membership_split",
    global_hierarchy_cut_no_split = "not_applicable_no_condition_split"
  )

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
    partition_policy = partition_policy,
    partition_schema_version = partition_schema_version,
    partition_diagnostics = if (length(partition_diagnostics)) {
      do.call(rbind, partition_diagnostics)
    } else {
      data.frame()
    },
    partition_repairs = if (length(partition_repairs)) {
      do.call(rbind, partition_repairs)
    } else {
      data.frame()
    },
    min_metacell_size = min.metacell.size,
    min_metacells_per_condition = min.metacells.per.condition,
    graph_scope = "independent_WNN_by_cell.graph.group",
    condition_scope = condition_scope,
    graph_method = "SCimplify_for_Seurat_multimodal_WNN_walktrap",
    modality_weighting = "adaptive_WNN_within_graph_group"
  )
  if (isTRUE(return.group.results)) result$group.results <- group_results
  result
}
