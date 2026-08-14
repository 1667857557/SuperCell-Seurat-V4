.SCConditionHierarchyCut <- function(hierarchy, cell_ids, k) {
  k <- as.integer(k)[1L]
  if (is.na(k) || k < 1L || k > length(cell_ids)) {
    stop("Shared hierarchy cut size is outside the valid cell range.",
         call. = FALSE)
  }
  cut <- tryCatch(
    igraph::cut_at(hierarchy, no = k),
    error = function(error) {
      stop("Walktrap hierarchy cannot be cut at k=", k, ": ",
           conditionMessage(error), call. = FALSE)
    }
  )
  if (length(cut) != length(cell_ids)) {
    stop("Walktrap hierarchy membership length does not match graph-group cells.",
         call. = FALSE)
  }
  if (!is.null(names(cut))) {
    index <- match(cell_ids, names(cut))
    if (anyNA(index)) {
      stop("Walktrap hierarchy membership cannot be aligned to graph-group cells.",
           call. = FALSE)
    }
    cut <- cut[index]
  }
  cut <- as.character(cut)
  if (anyNA(cut) || any(!nzchar(cut))) {
    stop("Walktrap hierarchy returned missing community memberships.",
         call. = FALSE)
  }
  stats::setNames(cut, cell_ids)
}

.SCConditionCutSummary <- function(
    hierarchy, cell_ids, condition, condition_value, k) {
  cut <- .SCConditionHierarchyCut(hierarchy, cell_ids, k)
  selected <- cell_ids[condition[cell_ids] == condition_value]
  sizes <- table(cut[selected])
  list(
    k = as.integer(k),
    membership = cut,
    n_metacells = as.integer(length(sizes)),
    min_size = if (length(sizes)) as.integer(min(sizes)) else 0L,
    max_size = if (length(sizes)) as.integer(max(sizes)) else 0L,
    mean_size = if (length(sizes)) mean(as.numeric(sizes)) else NA_real_
  )
}

.SCFindConditionHierarchyCut <- function(
    hierarchy, cell_ids, condition, condition_value, gamma,
    min.metacell.size = 1L, min.metacells.per.condition = 1L) {
  selected <- cell_ids[condition[cell_ids] == condition_value]
  n_cells <- length(selected)
  if (!n_cells) {
    stop("Condition `", condition_value,
         "` has no cells in the graph group.", call. = FALSE)
  }
  min.metacell.size <- as.integer(min.metacell.size)[1L]
  min.metacells.per.condition <- as.integer(min.metacells.per.condition)[1L]
  if (n_cells < min.metacell.size * min.metacells.per.condition) {
    stop(
      "Condition `", condition_value, "` is infeasible: ", n_cells,
      " cells cannot form ", min.metacells.per.condition,
      " metacells with at least ", min.metacell.size, " cells each.",
      call. = FALSE
    )
  }

  nominal_target <- max(1L, as.integer(floor(n_cells / gamma)))
  size_limited_target <- as.integer(floor(n_cells / min.metacell.size))
  target_metacells <- min(
    size_limited_target,
    max(min.metacells.per.condition, nominal_target)
  )

  evaluate <- function(k) {
    tryCatch(
      .SCConditionCutSummary(
        hierarchy = hierarchy,
        cell_ids = cell_ids,
        condition = condition,
        condition_value = condition_value,
        k = k
      ),
      error = function(error) NULL
    )
  }
  feasible <- function(summary) {
    !is.null(summary) &&
      summary$n_metacells <= target_metacells &&
      summary$min_size >= min.metacell.size
  }

  low <- 1L
  high <- length(cell_ids)
  base <- evaluate(low)
  if (!feasible(base)) {
    stop(
      "Shared Walktrap hierarchy cannot satisfy the minimum metacell size for ",
      "condition `", condition_value, "` even at its coarsest cut.",
      call. = FALSE
    )
  }
  while (low < high) {
    mid <- as.integer(floor((low + high + 1L) / 2L))
    candidate <- evaluate(mid)
    if (feasible(candidate)) {
      low <- mid
    } else {
      high <- mid - 1L
    }
  }
  chosen <- .SCConditionCutSummary(
    hierarchy = hierarchy,
    cell_ids = cell_ids,
    condition = condition,
    condition_value = condition_value,
    k = low
  )
  if (chosen$n_metacells < min.metacells.per.condition) {
    stop(
      "Shared Walktrap hierarchy is hierarchy-infeasible for condition `",
      condition_value, "`: the finest legal cut yields ",
      chosen$n_metacells, " metacells, below min.metacells.per.condition=",
      min.metacells.per.condition, ".", call. = FALSE
    )
  }

  sizes <- table(chosen$membership[selected])
  diagnostic <- data.frame(
    condition = condition_value,
    n_cells = as.integer(n_cells),
    requested_gamma = as.numeric(gamma),
    nominal_target_metacells = nominal_target,
    size_limited_target_metacells = size_limited_target,
    target_metacells = as.integer(target_metacells),
    selected_shared_cut_k = chosen$k,
    realized_metacells = chosen$n_metacells,
    min_realized_metacell_size = as.integer(min(sizes)),
    max_realized_metacell_size = as.integer(max(sizes)),
    mean_realized_metacell_size = mean(as.numeric(sizes)),
    min_metacell_size = min.metacell.size,
    min_metacells_per_condition = min.metacells.per.condition,
    feasibility_status = "ok",
    stringsAsFactors = FALSE
  )
  list(
    k = chosen$k,
    membership = chosen$membership,
    diagnostic = diagnostic
  )
}

.SCPartitionSharedHierarchyByCondition <- function(
    hierarchy, cell_ids, condition, gamma,
    min.metacell.size = 1L, min.metacells.per.condition = 1L) {
  condition <- .SCAlignGraphGrouping(
    condition, cell_ids, "cell.split.condition"
  )
  condition_levels <- sort(unique(unname(condition)))
  final_key <- parent_key <- shared_community <- stats::setNames(
    character(length(cell_ids)), cell_ids
  )
  shared_cut_k <- stats::setNames(integer(length(cell_ids)), cell_ids)
  diagnostics <- vector("list", length(condition_levels))
  names(diagnostics) <- condition_levels

  for (i in seq_along(condition_levels)) {
    condition_value <- condition_levels[[i]]
    selected <- .SCFindConditionHierarchyCut(
      hierarchy = hierarchy,
      cell_ids = cell_ids,
      condition = condition,
      condition_value = condition_value,
      gamma = gamma,
      min.metacell.size = min.metacell.size,
      min.metacells.per.condition = min.metacells.per.condition
    )
    cells_here <- cell_ids[condition[cell_ids] == condition_value]
    community <- unname(selected$membership[cells_here])
    final_key[cells_here] <- paste(condition_value, community, sep = "\001")
    parent_key[cells_here] <- paste0(
      "K", selected$k, "::C", community
    )
    shared_community[cells_here] <- community
    shared_cut_k[cells_here] <- selected$k
    diagnostics[[i]] <- selected$diagnostic
  }
  if (any(!nzchar(final_key)) || any(!nzchar(parent_key)) ||
      any(!nzchar(shared_community)) || any(shared_cut_k < 1L)) {
    stop("Condition-constrained hierarchy partition did not assign every cell.",
         call. = FALSE)
  }
  list(
    final_key = final_key,
    parent_key = parent_key,
    shared_cut_k = shared_cut_k,
    shared_community_id = shared_community,
    diagnostics = do.call(rbind, diagnostics)
  )
}
