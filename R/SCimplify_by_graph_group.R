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
#' Builds one independent multimodal WNN graph for each value of
#' `cell.graph.group`. All conditions inside a graph group are included in the
#' same WNN construction and Walktrap clustering. `cell.split.condition` is
#' applied only after clustering, so final memberships are condition-pure while
#' retaining a shared cross-condition graph within each broad cell type.
#'
#' @param seurat A preprocessed Seurat object containing the requested assays
#'   and reductions.
#' @param cell.graph.group Cell-level graph grouping vector, typically a broad
#'   cell-type annotation. A separate WNN graph is built for each value.
#' @param cell.split.condition Optional condition vector. Conditions are pooled
#'   during WNN construction and used only to split memberships after graph
#'   clustering.
#' @param k.knn Number of neighbours passed to [SCimplify_for_Seurat()].
#' @param kith Optional neighbourhood rank passed to
#'   [SCimplify_for_Seurat()].
#' @param kernel Whether to use the SuperCell kernel-weighted graph.
#' @param gamma Target average number of cells per metacell. The canonical
#'   default is 30.
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
#'   parent WNN-cluster memberships, a cell-level membership table, metacell
#'   sizes, group-specific hierarchies, and graph-construction provenance.
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
  k.knn <- as.integer(k.knn)[1L]
  gamma <- suppressWarnings(as.numeric(gamma)[1L])
  seed <- as.integer(seed)[1L]
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
  group_results <- vector("list", length(group_levels))
  names(group_results) <- group_levels
  hierarchies <- vector("list", length(group_levels))
  names(hierarchies) <- group_levels

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
    parent_membership[cells_i] <- paste0("G", i, "::", local)
    hierarchies[[i]] <- result_i$h_membership
    if (isTRUE(return.group.results)) {
      result_i$graph.group <- graph_level
      result_i$input.cells <- cells_i
      group_results[[i]] <- result_i
    }
  }

  if (any(!nzchar(parent_membership))) {
    stop("SuperCell did not assign every input cell.", call. = FALSE)
  }
  final_key <- if (is.null(split_condition)) {
    parent_membership[cell_ids]
  } else {
    paste(parent_membership[cell_ids], split_condition[cell_ids], sep = "\001")
  }
  id_map <- .SCStableMetacellIds(final_key)
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
    stringsAsFactors = FALSE
  )
  if (!is.null(split_condition)) {
    membership_table$condition <- unname(split_condition[cell_ids])
  }
  supercell_size <- vapply(membership_groups, length, integer(1))

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
    graph_scope = "independent_WNN_by_cell.graph.group",
    condition_scope = "joint_within_graph_group_then_membership_split",
    graph_method = "SCimplify_for_Seurat_multimodal_WNN_walktrap",
    modality_weighting = "adaptive_WNN_within_graph_group"
  )
  if (isTRUE(return.group.results)) result$group.results <- group_results
  result
}
