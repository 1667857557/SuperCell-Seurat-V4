test_that("graph groups are independent while conditions share each graph", {
  cell_id <- c(
    "A_C_1", "A_D_1", "A_C_2", "A_D_2", "A_C_3", "A_D_3",
    "B_C_1", "B_D_1", "B_C_2", "B_D_2", "B_C_3", "B_D_3"
  )
  graph_group <- stats::setNames(rep(c("A", "B"), each = 6L), cell_id)
  condition <- stats::setNames(rep(c("C", "D"), 6L), cell_id)
  X <- rbind(
    c(0, 0), c(0.01, 0), c(1, 0), c(1.01, 0), c(2, 0), c(2.01, 0),
    c(0, 0.001), c(0.01, 0.001), c(1, 0.001), c(1.01, 0.001),
    c(2, 0.001), c(2.01, 0.001)
  )
  rownames(X) <- cell_id
  colnames(X) <- c("dim1", "dim2")

  result <- SCimplify_by_graph_group_from_embedding(
    X = X,
    cell.graph.group = graph_group,
    cell.split.condition = condition,
    gamma = 2,
    k.knn = 1,
    n.pc = 1:2,
    seed = 11,
    return.singlecell.NW = TRUE,
    return.hierarchical.structure = TRUE,
    return.group.results = TRUE
  )

  expect_identical(result$graph_scope, "independent_by_cell.graph.group")
  expect_identical(
    result$condition_scope,
    "joint_within_graph_group_then_membership_split"
  )
  expect_setequal(names(result$group.results), c("A", "B"))

  combined_edges <- igraph::as_edgelist(result$graph.singlecell, names = TRUE)
  expect_true(nrow(combined_edges) > 0L)
  expect_true(all(
    graph_group[combined_edges[, 1L]] == graph_group[combined_edges[, 2L]]
  ))
  expect_true(any(
    condition[combined_edges[, 1L]] != condition[combined_edges[, 2L]]
  ))

  membership_groups <- split(names(result$membership), result$membership)
  expect_true(all(vapply(membership_groups, function(ids) {
    length(unique(graph_group[ids])) == 1L
  }, logical(1))))
  expect_true(all(vapply(membership_groups, function(ids) {
    length(unique(condition[ids])) == 1L
  }, logical(1))))
})

test_that("graph grouping vectors are aligned by cell ID", {
  X <- matrix(seq_len(16), nrow = 8L)
  rownames(X) <- paste0("cell", seq_len(nrow(X)))
  graph_group <- stats::setNames(
    rep(c("type1", "type2"), each = 4L), rev(rownames(X))
  )
  condition <- stats::setNames(
    rep(c("control", "treated"), 4L), rev(rownames(X))
  )
  expect_silent(
    SCimplify_by_graph_group_from_embedding(
      X,
      cell.graph.group = graph_group,
      cell.split.condition = condition,
      gamma = 2,
      k.knn = 1,
      n.pc = 1:2,
      return.singlecell.NW = FALSE,
      return.hierarchical.structure = FALSE
    )
  )
})
