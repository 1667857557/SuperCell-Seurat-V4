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
  ids <- paste0("c", 1:5)
  W <- matrix(0, 5, 5, dimnames = list(ids, ids))
  W["c1", "c2"] <- W["c2", "c1"] <- 10
  W["c1", "c3"] <- W["c3", "c1"] <- 0.1
  W["c3", "c4"] <- W["c4", "c3"] <- 1
  W["c2", "c5"] <- W["c5", "c2"] <- 1
  graph <- igraph::graph_from_adjacency_matrix(
    W, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  provisional <- stats::setNames(c("A1", "B1", "A2", "A2", "B1"), ids)
  stratum <- stats::setNames(c("A", "B", "A", "A", "B"), ids)
  repaired <- SuperCell:::.SCRepairSplitMembership(
    graph = graph,
    provisional_membership = provisional,
    stratum = stratum,
    min_metacell_size = 2L,
    min_metacells_per_stratum = 1L,
    min_merge_affinity = 0,
    unresolved_small_policy = "error"
  )
  expect_identical(repaired$membership[["c1"]], "A2")
  expect_identical(repaired$membership[["c2"]], "B1")
  expect_identical(repaired$membership[["c5"]], "B1")
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

test_that("SCimplify_for_Seurat exposes the original graph only for transient compact return", {
  expect_true("return.graph" %in% names(formals(SCimplify_for_Seurat)))
  expect_identical(eval(formals(SCimplify_for_Seurat)$return.graph), FALSE)
  body_text <- paste(deparse(body(SCimplify_for_Seurat)), collapse = "\n")
  expect_match(body_text, "seurat.mc$graph <- graph", fixed = TRUE)
  expect_match(body_text, "requires `return.seurat = FALSE`", fixed = TRUE)
})

test_that("grouped repair obtains the graph from SCimplify_for_Seurat and discards it", {
  body_text <- paste(deparse(body(SCimplify_by_graph_group)), collapse = "\n")
  expect_match(body_text, "result_i <- SCimplify_for_Seurat", fixed = TRUE)
  expect_match(body_text, "return.graph = min_metacell_size > 1L", fixed = TRUE)
  expect_match(body_text, "graph = graph_i", fixed = TRUE)
  expect_match(body_text, "result_i$graph <- NULL", fixed = TRUE)
  expect_false(grepl(".SCBuildGraphMembership", body_text, fixed = TRUE))
})
