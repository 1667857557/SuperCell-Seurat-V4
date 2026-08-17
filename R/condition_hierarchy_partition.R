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

.SCConditionHierarchyRange <- function(hierarchy, n_vertices) {
  if (!isTRUE(igraph::is_hierarchical(hierarchy))) {
    stop("Condition-constrained partitioning requires a hierarchical community result.",
         call. = FALSE)
  }
  merge_matrix <- igraph::merges(hierarchy)
  n_merges <- if (is.null(merge_matrix)) 0L else nrow(merge_matrix)
  min_k <- as.integer(n_vertices - n_merges)
  if (is.na(min_k) || min_k < 1L || min_k > n_vertices) {
    stop("Walktrap hierarchy exposes an invalid merge range.", call. = FALSE)
  }
  c(min = min_k, max = as.integer(n_vertices))
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
    n_cells,
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
    !is.null(summary) && summary$n_metacells <= target_metacells
  }

  hierarchy_range <- .SCConditionHierarchyRange(hierarchy, length(cell_ids))
  low <- unname(hierarchy_range[["min"]])
  high <- unname(hierarchy_range[["max"]])
  base <- evaluate(low)
  if (!feasible(base)) {
    stop(
      "Shared Walktrap hierarchy cannot satisfy the requested condition-specific ",
      "metacell target for condition `", condition_value, "` at its coarsest ",
      "valid cut (k=", low, ").", call. = FALSE
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
      condition_value, "`: the gamma-resolution cut yields ",
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
    hierarchy_min_cut_k = as.integer(hierarchy_range[["min"]]),
    hierarchy_max_cut_k = as.integer(hierarchy_range[["max"]]),
    selected_shared_cut_k = chosen$k,
    initial_metacells = chosen$n_metacells,
    initial_min_metacell_size = as.integer(min(sizes)),
    initial_max_metacell_size = as.integer(max(sizes)),
    initial_mean_metacell_size = mean(as.numeric(sizes)),
    min_metacell_size = min.metacell.size,
    min_metacells_per_condition = min.metacells.per.condition,
    stringsAsFactors = FALSE
  )
  list(
    k = chosen$k,
    membership = chosen$membership,
    diagnostic = diagnostic
  )
}

.SCConditionHierarchyInitialNodes <- function(
    hierarchy, cell_ids, initial_k) {
  n_vertices <- length(cell_ids)
  merge_matrix <- igraph::merges(hierarchy)
  if (is.null(merge_matrix)) {
    merge_matrix <- matrix(integer(), nrow = 0L, ncol = 2L)
  }
  merge_matrix <- matrix(
    as.integer(merge_matrix), nrow = nrow(merge_matrix), ncol = 2L
  )
  n_merges <- nrow(merge_matrix)
  steps <- as.integer(n_vertices - initial_k)
  if (is.na(steps) || steps < 0L || steps > n_merges) {
    stop("Initial shared hierarchy cut is outside the available merge range.",
         call. = FALSE)
  }

  parent <- integer(n_vertices + steps)
  if (steps > 0L) {
    for (r in seq_len(steps)) {
      children <- merge_matrix[r, ]
      parent_node <- n_vertices + r
      if (anyNA(children) || any(children < 1L) ||
          any(children >= parent_node)) {
        stop("Walktrap merge matrix uses an unsupported community encoding.",
             call. = FALSE)
      }
      parent[children] <- parent_node
    }
  }
  active_node <- vapply(seq_len(n_vertices), function(vertex) {
    node <- vertex
    while (node <= length(parent) && parent[[node]] > 0L) {
      node <- parent[[node]]
    }
    as.integer(node)
  }, integer(1))
  stats::setNames(active_node, cell_ids)

  list(
    active_node = active_node,
    merge_matrix = merge_matrix,
    n_vertices = as.integer(n_vertices),
    n_merges = as.integer(n_merges),
    initial_steps = steps
  )
}

.SCChooseHierarchyRepairTarget <- function(
    source, candidates, group_sizes, gamma) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates <- unique(candidates)
  if (!length(candidates)) return(NA_character_)
  resulting_size <- group_sizes[[source]] + unname(group_sizes[candidates])
  score <- data.frame(
    candidate = candidates,
    distance_to_gamma = abs(resulting_size - gamma),
    resulting_size = resulting_size,
    stringsAsFactors = FALSE
  )
  score <- score[order(
    score$distance_to_gamma,
    score$resulting_size,
    score$candidate
  ), , drop = FALSE]
  score$candidate[[1L]]
}

.SCRepairSmallConditionMetacells <- function(
    hierarchy, cell_ids, condition, condition_value,
    initial_membership, initial_k, gamma,
    min.metacell.size = 1L, min.metacells.per.condition = 1L) {
  selected <- cell_ids[condition[cell_ids] == condition_value]
  initial_labels <- unname(initial_membership[selected])
  if (!length(selected) || anyNA(initial_labels) || any(!nzchar(initial_labels))) {
    stop("Initial condition hierarchy membership is incomplete.", call. = FALSE)
  }

  topology <- .SCConditionHierarchyInitialNodes(
    hierarchy = hierarchy,
    cell_ids = cell_ids,
    initial_k = initial_k
  )
  initial_groups <- split(selected, initial_labels)
  group_labels <- sort(names(initial_groups))
  initial_groups <- initial_groups[group_labels]
  group_ids <- paste0("G", seq_along(initial_groups))
  names(initial_groups) <- group_ids

  group_members <- initial_groups
  group_origins <- lapply(group_labels, function(label) label)
  names(group_origins) <- group_ids
  group_sizes <- stats::setNames(
    vapply(group_members, length, integer(1)), group_ids
  )
  group_active <- stats::setNames(rep(TRUE, length(group_ids)), group_ids)

  node_groups <- vector(
    "list", topology$n_vertices + topology$n_merges
  )
  for (j in seq_along(group_ids)) {
    gid <- group_ids[[j]]
    cells_j <- group_members[[gid]]
    nodes_j <- unique(unname(topology$active_node[cells_j]))
    if (length(nodes_j) != 1L) {
      stop("An initial condition metacell spans multiple hierarchy nodes.",
           call. = FALSE)
    }
    node <- nodes_j[[1L]]
    node_groups[[node]] <- c(node_groups[[node]], gid)
  }

  repair_records <- list()
  repair_iteration <- 0L
  if (min.metacell.size > 1L &&
      topology$initial_steps < topology$n_merges) {
    rows <- seq.int(topology$initial_steps + 1L, topology$n_merges)
    for (r in rows) {
      children <- topology$merge_matrix[r, ]
      parent_node <- topology$n_vertices + r
      left <- unique(node_groups[[children[[1L]]]])
      right <- unique(node_groups[[children[[2L]]]])
      left <- left[group_active[left] %in% TRUE]
      right <- right[group_active[right] %in% TRUE]

      small_left <- left[group_sizes[left] < min.metacell.size]
      small_right <- right[group_sizes[right] < min.metacell.size]
      edge_from <- edge_to <- character()

      if (length(small_left) && length(right)) {
        for (source in sort(small_left)) {
          target <- .SCChooseHierarchyRepairTarget(
            source, right, group_sizes, gamma
          )
          if (!is.na(target)) {
            edge_from <- c(edge_from, source)
            edge_to <- c(edge_to, target)
          }
        }
      }
      if (length(small_right) && length(left)) {
        for (source in sort(small_right)) {
          target <- .SCChooseHierarchyRepairTarget(
            source, left, group_sizes, gamma
          )
          if (!is.na(target)) {
            edge_from <- c(edge_from, source)
            edge_to <- c(edge_to, target)
          }
        }
      }

      if (length(edge_from)) {
        vertices <- sort(unique(c(edge_from, edge_to)))
        adjacency <- stats::setNames(vector("list", length(vertices)), vertices)
        for (e in seq_along(edge_from)) {
          a <- edge_from[[e]]
          b <- edge_to[[e]]
          adjacency[[a]] <- unique(c(adjacency[[a]], b))
          adjacency[[b]] <- unique(c(adjacency[[b]], a))
        }
        seen <- stats::setNames(rep(FALSE, length(vertices)), vertices)
        components <- list()
        for (start in vertices) {
          if (seen[[start]]) next
          queue <- start
          component <- character()
          seen[[start]] <- TRUE
          while (length(queue)) {
            current <- queue[[1L]]
            queue <- queue[-1L]
            component <- c(component, current)
            neighbours <- adjacency[[current]]
            unseen <- neighbours[!seen[neighbours]]
            if (length(unseen)) {
              seen[unseen] <- TRUE
              queue <- c(queue, unseen)
            }
          }
          if (length(component) > 1L) {
            components[[length(components) + 1L]] <- sort(component)
          }
        }

        for (component in components) {
          component <- component[group_active[component] %in% TRUE]
          if (length(component) < 2L) next
          valid_before <- component[
            group_sizes[component] >= min.metacell.size
          ]
          if (length(valid_before) > 1L) {
            stop(
              "Local hierarchy repair would merge multiple already-valid metacells; ",
              "this violates the minimum-change partition contract.",
              call. = FALSE
            )
          }
          survivor <- if (length(valid_before) == 1L) {
            valid_before[[1L]]
          } else {
            sort(component)[[1L]]
          }
          absorbed <- setdiff(component, survivor)
          if (!length(absorbed)) next
          active_count <- sum(group_active)
          if (active_count - length(absorbed) < min.metacells.per.condition) {
            stop(
              "Condition `", condition_value,
              "` is hierarchy-infeasible: local minimum-size repair would ",
              "reduce the final metacell count below min.metacells.per.condition=",
              min.metacells.per.condition, ".", call. = FALSE
            )
          }
          sizes_before <- unname(group_sizes[component])
          names(sizes_before) <- component
          group_members[[survivor]] <- unlist(
            group_members[component], use.names = FALSE
          )
          group_origins[[survivor]] <- sort(unique(unlist(
            group_origins[component], use.names = FALSE
          )))
          group_sizes[[survivor]] <- sum(group_sizes[component])
          group_active[absorbed] <- FALSE
          repair_iteration <- repair_iteration + 1L
          repair_records[[length(repair_records) + 1L]] <- data.frame(
            condition = condition_value,
            repair_iteration = repair_iteration,
            hierarchy_merge_k = as.integer(topology$n_vertices - r),
            survivor_group = survivor,
            absorbed_groups = paste(absorbed, collapse = ","),
            component_groups = paste(component, collapse = ","),
            component_sizes_before = paste(
              paste(names(sizes_before), sizes_before, sep = ":"),
              collapse = ","
            ),
            size_after = as.integer(group_sizes[[survivor]]),
            stringsAsFactors = FALSE
          )
        }
      }

      parent_groups <- unique(c(left, right))
      parent_groups <- parent_groups[group_active[parent_groups] %in% TRUE]
      node_groups[[parent_node]] <- parent_groups
    }
  }

  active_groups <- names(group_active)[group_active]
  final_sizes <- unname(group_sizes[active_groups])
  if (any(final_sizes < min.metacell.size)) {
    unresolved <- active_groups[final_sizes < min.metacell.size]
    stop(
      "Condition `", condition_value,
      "` is hierarchy-infeasible: ", length(unresolved),
      " local metacell(s) remain below min.metacell.size=",
      min.metacell.size,
      " after exhausting same-condition merges in the shared Walktrap hierarchy.",
      call. = FALSE
    )
  }
  if (length(active_groups) < min.metacells.per.condition) {
    stop(
      "Condition `", condition_value,
      "` is hierarchy-infeasible after local repair: ", length(active_groups),
      " metacells remain, below min.metacells.per.condition=",
      min.metacells.per.condition, ".", call. = FALSE
    )
  }

  final_key <- stats::setNames(character(length(selected)), selected)
  final_group_label <- stats::setNames(character(length(active_groups)), active_groups)
  for (gid in active_groups) {
    origins <- sort(unique(group_origins[[gid]]))
    label <- paste(condition_value, paste(origins, collapse = "+"), sep = "\001")
    final_group_label[[gid]] <- label
    final_key[group_members[[gid]]] <- label
  }
  if (any(!nzchar(final_key))) {
    stop("Local hierarchy repair did not assign every condition cell.",
         call. = FALSE)
  }

  repairs <- if (length(repair_records)) {
    do.call(rbind, repair_records)
  } else {
    data.frame(
      condition = character(),
      repair_iteration = integer(),
      hierarchy_merge_k = integer(),
      survivor_group = character(),
      absorbed_groups = character(),
      component_groups = character(),
      component_sizes_before = character(),
      size_after = integer(),
      stringsAsFactors = FALSE
    )
  }
  list(
    final_key = final_key,
    repairs = repairs,
    initial_metacells = length(initial_groups),
    repaired_initial_small_metacells = sum(
      vapply(initial_groups, length, integer(1)) < min.metacell.size
    ),
    realized_metacells = length(active_groups),
    sizes = stats::setNames(final_sizes, final_group_label[active_groups])
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
  repairs <- vector("list", length(condition_levels))
  names(diagnostics) <- names(repairs) <- condition_levels

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
    repaired <- .SCRepairSmallConditionMetacells(
      hierarchy = hierarchy,
      cell_ids = cell_ids,
      condition = condition,
      condition_value = condition_value,
      initial_membership = selected$membership,
      initial_k = selected$k,
      gamma = gamma,
      min.metacell.size = min.metacell.size,
      min.metacells.per.condition = min.metacells.per.condition
    )
    final_key[cells_here] <- unname(repaired$final_key[cells_here])
    parent_key[cells_here] <- paste0(
      "K", selected$k, "::C", community
    )
    shared_community[cells_here] <- community
    shared_cut_k[cells_here] <- selected$k

    final_sizes <- unname(repaired$sizes)
    diagnostic <- selected$diagnostic
    diagnostic$repaired_initial_small_metacells <-
      repaired$repaired_initial_small_metacells
    diagnostic$repair_merge_events <- nrow(repaired$repairs)
    diagnostic$repair_metacell_reduction <-
      repaired$initial_metacells - repaired$realized_metacells
    diagnostic$realized_metacells <- repaired$realized_metacells
    diagnostic$min_realized_metacell_size <- as.integer(min(final_sizes))
    diagnostic$max_realized_metacell_size <- as.integer(max(final_sizes))
    diagnostic$mean_realized_metacell_size <- mean(final_sizes)
    diagnostic$feasibility_status <- "ok"
    diagnostic$repair_status <- if (nrow(repaired$repairs)) {
      "local_hierarchy_repair_applied"
    } else {
      "no_repair_needed"
    }
    diagnostics[[i]] <- diagnostic
    repairs[[i]] <- repaired$repairs
  }
  if (any(!nzchar(final_key)) || any(!nzchar(parent_key)) ||
      any(!nzchar(shared_community)) || any(shared_cut_k < 1L)) {
    stop("Condition-constrained hierarchy partition did not assign every cell.",
         call. = FALSE)
  }
  repair_table <- repairs[vapply(repairs, nrow, integer(1)) > 0L]
  list(
    final_key = final_key,
    parent_key = parent_key,
    shared_cut_k = shared_cut_k,
    shared_community_id = shared_community,
    diagnostics = do.call(rbind, diagnostics),
    repairs = if (length(repair_table)) {
      do.call(rbind, repair_table)
    } else {
      data.frame()
    }
  )
}
