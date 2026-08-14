test_that("small metacells merge only by original shared-WNN affinity", {
  ids <- paste0("c", 1:6)
  W <- matrix(0, 6, 6, dimnames = list(ids, ids))
  W["c1", "c2"] <- W["c2", "c1"] <- 4
  W["c1", "c4"] <- W["c4", "c1"] <- 0.1
  W["c2", "c3"] <- W["c3", "c2"] <- 2
  W["c4", "c5"] <- W["c5", "c4"] <- 2
  W["c5", "c6"] <- W["c6", "c5"] <- 2
  graph <- igraph::graph_from_adjacency_matrix(
    W, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  provisional <- stats::setNames(
    c("M1", "M2", "M2", "M3", "M3", "M3"), ids
  )
  stratum <- stats::setNames(rep("A", 6), ids)
  repaired <- SuperCell:::.SCRepairSplitMembership(
    graph = graph,
    provisional_membership = provisional,
    stratum = stratum,
    min_metacell_size = 2L,
    min_metacells_per_stratum = 2L,
    min_merge_affinity = 0.05,
    unresolved_small_policy = "error"
  )
  expect_identical(repaired$membership[["c1"]], repaired$membership[["c2"]])
  expect_equal(length(unique(repaired$membership)), 2L)
  expect_true(all(table(repaired$membership) >= 2L))
  expect_equal(nrow(repaired$unresolved), 0L)
  expect_equal(nrow(repaired$merge_diagnostics), 1L)
})

test_that("repair never crosses a condition stratum", {
  ids <- paste0("c", 1:4)
  W <- matrix(0, 4, 4, dimnames = list(ids, ids))
  W["c1", "c2"] <- W["c2", "c1"] <- 10
  W["c1", "c3"] <- W["c3", "c1"] <- 0.1
  W["c3", "c4"] <- W["c4", "c3"] <- 1
  graph <- igraph::graph_from_adjacency_matrix(
    W, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  provisional <- stats::setNames(c("A1", "B1", "A2", "A2"), ids)
  stratum <- stats::setNames(c("A", "B", "A", "A"), ids)
  repaired <- SuperCell:::.SCRepairSplitMembership(
    graph = graph,
    provisional_membership = provisional,
    stratum = stratum,
    min_metacell_size = 2L,
    min_metacells_per_stratum = 1L,
    min_merge_affinity = 0,
    unresolved_small_policy = "keep"
  )
  expect_identical(repaired$membership[["c1"]], "A2")
  expect_identical(repaired$membership[["c2"]], "B1")
})

test_that("repair requires an explicit affinity threshold when enabled", {
  ids <- c("a", "b")
  graph <- igraph::make_full_graph(2)
  igraph::V(graph)$name <- ids
  expect_error(
    SuperCell:::.SCRepairSplitMembership(
      graph,
      stats::setNames(c("M1", "M2"), ids),
      stats::setNames(rep("A", 2), ids),
      min_metacell_size = 2L,
      min_metacells_per_stratum = 1L,
      min_merge_affinity = NULL
    ),
    "min_merge_affinity",
    fixed = TRUE
  )
})
