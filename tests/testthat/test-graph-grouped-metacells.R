test_that("graph groups use one joint-condition multimodal call each", {
  counts <- matrix(
    seq_len(32), nrow = 4,
    dimnames = list(paste0("gene", 1:4), paste0("cell", 1:8))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(24), nrow = 8,
      dimnames = list(colnames(object), paste0("PC_", 1:3))
    ), key = "PC_", assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      rev(seq_len(24)), nrow = 8,
      dimnames = list(colnames(object), paste0("LSI_", 1:3))
    ), key = "LSI_", assay = "ATAC"
  )
  graph_group <- stats::setNames(rep(c("T", "B"), each = 4), colnames(object))
  condition <- stats::setNames(rep(c("control", "treated"), 4), colnames(object))
  calls <- list()

  testthat::local_mocked_bindings(
    SCimplify_for_Seurat = function(seurat, assay, reduction, dims, ...) {
      calls[[length(calls) + 1L]] <<- list(
        cells = colnames(seurat), assay = assay,
        reduction = reduction, dims = dims
      )
      local <- rep(c("1", "2"), each = 2, length.out = ncol(seurat))
      list(
        membership = stats::setNames(local, colnames(seurat)),
        h_membership = list(membership = local)
      )
    },
    .package = "SuperCell"
  )

  result <- SCimplify_by_graph_group(
    seurat = object,
    cell.graph.group = graph_group,
    cell.split.condition = condition,
    assay = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims = list(1:3, 1:3),
    gamma = 2,
    return.group.results = TRUE
  )

  expect_length(calls, 2L)
  expect_true(all(vapply(calls, function(x) {
    length(unique(graph_group[x$cells])) == 1L &&
      length(unique(condition[x$cells])) == 2L
  }, logical(1))))
  expect_true(all(vapply(calls, function(x) {
    identical(x$assay, c("RNA", "ATAC")) && length(x$reduction) == 2L
  }, logical(1))))
  expect_identical(result$graph_scope,
                   "independent_WNN_by_cell.graph.group")
  expect_identical(result$condition_scope,
                   "joint_within_graph_group_then_membership_split")
  expect_identical(result$modality_weighting,
                   "adaptive_WNN_within_graph_group")
  membership_groups <- split(names(result$membership), result$membership)
  expect_true(all(vapply(membership_groups, function(ids) {
    length(unique(graph_group[ids])) == 1L
  }, logical(1))))
  expect_true(all(vapply(membership_groups, function(ids) {
    length(unique(condition[ids])) == 1L
  }, logical(1))))
  expect_setequal(result$membership_table$cell_id, colnames(object))
  expect_setequal(names(result$supercell_size), unique(result$membership))
})

test_that("graph grouping vectors align by cell name", {
  counts <- matrix(
    seq_len(16), nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", 1:8))
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(seq_len(16), nrow = 8,
      dimnames = list(colnames(object), paste0("PC_", 1:2))),
    key = "PC_", assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(rev(seq_len(16)), nrow = 8,
      dimnames = list(colnames(object), paste0("LSI_", 1:2))),
    key = "LSI_", assay = "ATAC"
  )
  group <- stats::setNames(rep(c("A", "B"), each = 4), rev(colnames(object)))
  condition <- stats::setNames(rep(c("C", "D"), 4), rev(colnames(object)))

  testthat::local_mocked_bindings(
    SCimplify_for_Seurat = function(seurat, ...) {
      list(
        membership = stats::setNames(
          rep("1", ncol(seurat)), colnames(seurat)
        ),
        h_membership = list()
      )
    },
    .package = "SuperCell"
  )

  result <- SCimplify_by_graph_group(
    object, group, condition,
    assay = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims = list(1:2, 1:2)
  )
  expect_identical(names(result$membership), colnames(object))
})

test_that("obsolete graph-group embedding API is removed", {
  expect_false(
    "SCimplify_by_graph_group_from_embedding" %in%
      getNamespaceExports("SuperCell")
  )
  expect_true("SCimplify_by_graph_group" %in% getNamespaceExports("SuperCell"))
})

test_that("invalid graph-group multimodal inputs fail early", {
  counts <- matrix(1, nrow = 2, ncol = 4,
                   dimnames = list(c("g1", "g2"), paste0("c", 1:4)))
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  expect_error(
    SCimplify_by_graph_group(
      object,
      cell.graph.group = c("A", "A", "B", "B")
    ),
    "Missing assay"
  )
  expect_error(
    .SCAlignGraphGrouping(c("A", "B"), colnames(object), "group"),
    "one value per input cell"
  )
})
