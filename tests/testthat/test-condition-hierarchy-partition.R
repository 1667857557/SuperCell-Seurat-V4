test_that("hierarchy-constrained API exposes hard feasibility controls", {
  args <- names(formals(SCimplify_by_graph_group))
  expect_true(all(c(
    "condition.partition", "min.metacell.size",
    "min.metacells.per.condition"
  ) %in% args))
  expect_identical(
    eval(formals(SCimplify_by_graph_group)$condition.partition)[[1L]],
    "legacy_post_split"
  )
  expect_identical(eval(formals(SCimplify_by_graph_group)$min.metacell.size), 1L)
  expect_identical(
    eval(formals(SCimplify_by_graph_group)$min.metacells.per.condition), 1L
  )
})

test_that("condition hierarchy cut enforces necessary feasibility before search", {
  graph <- igraph::make_ring(6)
  hierarchy <- igraph::cluster_walktrap(graph)
  cells <- paste0("c", seq_len(6))
  condition <- stats::setNames(rep("A", 6), cells)
  expect_error(
    SuperCell:::.SCFindConditionHierarchyCut(
      hierarchy = hierarchy,
      cell_ids = cells,
      condition = condition,
      condition_value = "A",
      gamma = 3,
      min.metacell.size = 4L,
      min.metacells.per.condition = 2L
    ),
    "infeasible"
  )
})

test_that("condition hierarchy cut uses one shared Walktrap hierarchy", {
  graph <- igraph::make_ring(12)
  hierarchy <- igraph::cluster_walktrap(graph)
  cells <- paste0("c", seq_len(12))
  condition <- stats::setNames(rep(c("A", "B"), each = 6), cells)
  partition <- SuperCell:::.SCPartitionSharedHierarchyByCondition(
    hierarchy = hierarchy,
    cell_ids = cells,
    condition = condition,
    gamma = 6,
    min.metacell.size = 2L,
    min.metacells.per.condition = 1L
  )
  expect_setequal(names(partition$final_key), cells)
  expect_true(all(partition$shared_cut_k >= 1L))
  groups <- split(cells, partition$final_key)
  expect_true(all(vapply(groups, function(ids) {
    length(unique(condition[ids])) == 1L && length(ids) >= 2L
  }, logical(1))))
  expect_setequal(partition$diagnostics$condition, c("A", "B"))
  expect_true(all(partition$diagnostics$feasibility_status == "ok"))
})
