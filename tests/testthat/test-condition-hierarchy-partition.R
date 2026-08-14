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

test_that("disconnected Walktrap hierarchy starts at its coarsest valid cut", {
  graph <- igraph::disjoint_union(igraph::make_ring(4), igraph::make_ring(4))
  hierarchy <- igraph::cluster_walktrap(graph)
  cells <- paste0("c", seq_len(8))
  names(igraph::V(graph)) <- cells
  range <- SuperCell:::.SCConditionHierarchyRange(hierarchy, length(cells))
  expect_equal(unname(range[["min"]]), 2L)
  condition <- stats::setNames(rep(c("A", "B"), each = 4), cells)
  partition <- SuperCell:::.SCPartitionSharedHierarchyByCondition(
    hierarchy = hierarchy,
    cell_ids = cells,
    condition = condition,
    gamma = 4,
    min.metacell.size = 2L,
    min.metacells.per.condition = 1L
  )
  expect_true(all(partition$diagnostics$hierarchy_min_cut_k == 2L))
  expect_true(all(partition$diagnostics$min_realized_metacell_size >= 2L))
})

test_that("hierarchy mode builds one shared hierarchy per graph group", {
  cells <- unlist(lapply(c("T", "B"), function(group) {
    unlist(lapply(c("A", "B", "C"), function(condition) {
      paste(group, condition, seq_len(4), sep = "_")
    }), use.names = FALSE)
  }), use.names = FALSE)
  counts <- matrix(
    seq_len(4L * length(cells)), nrow = 4L,
    dimnames = list(paste0("g", seq_len(4)), cells)
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(length(cells) * 3L), nrow = length(cells),
      dimnames = list(cells, paste0("PC_", seq_len(3)))
    ), key = "PC_", assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      rev(seq_len(length(cells) * 3L)), nrow = length(cells),
      dimnames = list(cells, paste0("LSI_", seq_len(3)))
    ), key = "LSI_", assay = "ATAC"
  )
  graph_group <- stats::setNames(sub("_.*$", "", cells), cells)
  condition <- stats::setNames(sub("^[^_]+_([^_]+)_.*$", "\\1", cells), cells)
  calls <- list()

  testthat::local_mocked_bindings(
    SCimplify_for_Seurat = function(seurat, ...) {
      ids <- sort(colnames(seurat))
      graph <- igraph::make_ring(length(ids))
      igraph::V(graph)$name <- ids
      hierarchy <- igraph::cluster_walktrap(graph)
      local <- igraph::cut_at(hierarchy, no = max(1L, length(ids) %/% 4L))
      names(local) <- ids
      calls[[length(calls) + 1L]] <<- ids
      list(membership = local, h_membership = hierarchy)
    },
    .package = "SuperCell"
  )

  result <- SCimplify_by_graph_group(
    seurat = object,
    cell.graph.group = graph_group,
    cell.split.condition = condition,
    condition.partition = "hierarchy_constrained",
    min.metacell.size = 2L,
    min.metacells.per.condition = 1L,
    assay = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims = list(1:3, 1:3),
    gamma = 4
  )

  expect_length(calls, 2L)
  expect_true(all(vapply(calls, function(ids) {
    length(unique(graph_group[ids])) == 1L &&
      setequal(unique(condition[ids]), c("A", "B", "C"))
  }, logical(1))))
  expect_identical(result$partition_policy, "hierarchy_constrained")
  expect_identical(
    result$partition_schema_version,
    "shared_walktrap_condition_cut_v1"
  )
  groups <- split(result$membership_table, result$membership_table$metacell_id)
  expect_true(all(vapply(groups, function(tab) {
    length(unique(tab$graph_group)) == 1L &&
      length(unique(tab$condition)) == 1L &&
      nrow(tab) >= 2L
  }, logical(1))))
})

test_that("hierarchy mode canonical IDs are invariant to cell order", {
  cells <- unlist(lapply(c("T", "B"), function(group) {
    unlist(lapply(c("A", "B"), function(condition) {
      paste(group, condition, seq_len(6), sep = "_")
    }), use.names = FALSE)
  }), use.names = FALSE)
  counts <- matrix(
    seq_len(3L * length(cells)), nrow = 3L,
    dimnames = list(paste0("g", seq_len(3)), cells)
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  pca_embedding <- matrix(
    seq_len(length(cells) * 3L), nrow = length(cells),
    dimnames = list(cells, paste0("PC_", seq_len(3)))
  )
  lsi_embedding <- pca_embedding[, 3:1, drop = FALSE]
  colnames(lsi_embedding) <- paste0("LSI_", seq_len(3))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca_embedding, key = "PC_", assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = lsi_embedding, key = "LSI_", assay = "ATAC"
  )
  graph_group <- stats::setNames(sub("_.*$", "", cells), cells)
  condition <- stats::setNames(sub("^[^_]+_([^_]+)_.*$", "\\1", cells), cells)

  testthat::local_mocked_bindings(
    SCimplify_for_Seurat = function(seurat, ...) {
      ids <- sort(colnames(seurat))
      graph <- igraph::make_ring(length(ids))
      igraph::V(graph)$name <- ids
      hierarchy <- igraph::cluster_walktrap(graph)
      local <- igraph::cut_at(hierarchy, no = 2L)
      names(local) <- ids
      list(membership = local, h_membership = hierarchy)
    },
    .package = "SuperCell"
  )

  run <- function(x) {
    SCimplify_by_graph_group(
      seurat = x,
      cell.graph.group = graph_group[colnames(x)],
      cell.split.condition = condition[colnames(x)],
      condition.partition = "hierarchy_constrained",
      min.metacell.size = 2L,
      min.metacells.per.condition = 1L,
      assay = c("RNA", "ATAC"),
      reduction = list("pca", "lsi"),
      dims = list(1:3, 1:3),
      gamma = 6
    )$membership
  }

  first <- run(object)
  within_group_reverse <- unlist(lapply(c("T", "B"), function(group) {
    rev(cells[graph_group[cells] == group])
  }), use.names = FALSE)
  second <- run(object[, within_group_reverse])
  expect_identical(first[sort(names(first))], second[sort(names(second))])
})
