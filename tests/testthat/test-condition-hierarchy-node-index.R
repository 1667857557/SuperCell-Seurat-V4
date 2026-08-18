test_that("initial hierarchy nodes preserve cell IDs for character indexing", {
  graph <- igraph::make_ring(24)
  hierarchy <- igraph::cluster_walktrap(graph)
  cells <- paste0("cell_", seq_len(24))

  topology <- SuperCell:::.SCConditionHierarchyInitialNodes(
    hierarchy = hierarchy,
    cell_ids = cells,
    initial_k = 6L
  )

  expect_identical(names(topology$active_node), cells)
  expect_false(anyNA(topology$active_node[cells]))
  expect_true(all(topology$active_node[cells] >= 1L))
  expect_true(all(
    topology$active_node[cells] <= topology$n_vertices + topology$n_merges
  ))
})

test_that("local hierarchy repair accepts barcode-indexed initial groups", {
  graph <- igraph::make_ring(40)
  hierarchy <- igraph::cluster_walktrap(graph)
  cells <- paste0("cell_", seq_len(40))
  initial_k <- 8L
  initial <- SuperCell:::.SCConditionHierarchyCut(
    hierarchy = hierarchy,
    cell_ids = cells,
    k = initial_k
  )

  groups <- split(cells, initial)
  ordered <- names(sort(vapply(groups, length, integer(1)), decreasing = TRUE))
  expect_gte(length(ordered), 4L)

  condition <- stats::setNames(rep("B", length(cells)), cells)
  for (label in ordered[1:3]) {
    condition[groups[[label]]] <- "A"
  }
  small_cell <- groups[[ordered[[4L]]]][[1L]]
  condition[[small_cell]] <- "A"

  expect_error(
    SuperCell:::.SCRepairSmallConditionMetacells(
      hierarchy = hierarchy,
      cell_ids = cells,
      condition = condition,
      condition_value = "A",
      initial_membership = initial,
      initial_k = initial_k,
      gamma = 5,
      min.metacell.size = 3L,
      min.metacells.per.condition = 1L
    ),
    NA
  )
})
