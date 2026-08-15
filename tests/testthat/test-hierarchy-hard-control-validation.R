test_that("hierarchy feasibility controls reject fractional values", {
  cells <- paste0("c", seq_len(4))
  counts <- matrix(
    seq_len(8), nrow = 2,
    dimnames = list(c("g1", "g2"), cells)
  )
  object <- SeuratObject::CreateSeuratObject(counts = counts)
  object[["ATAC"]] <- SeuratObject::CreateAssayObject(counts = counts)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      seq_len(8), nrow = 4,
      dimnames = list(cells, c("PC_1", "PC_2"))
    ),
    key = "PC_", assay = "RNA"
  )
  object[["lsi"]] <- SeuratObject::CreateDimReducObject(
    embeddings = matrix(
      rev(seq_len(8)), nrow = 4,
      dimnames = list(cells, c("LSI_1", "LSI_2"))
    ),
    key = "LSI_", assay = "ATAC"
  )

  common <- list(
    seurat = object,
    cell.graph.group = stats::setNames(rep("T", 4), cells),
    assay = c("RNA", "ATAC"),
    reduction = list("pca", "lsi"),
    dims = list(1:2, 1:2)
  )

  expect_error(
    do.call(
      SCimplify_by_graph_group,
      c(common, list(min.metacell.size = 1.5))
    ),
    "min.metacell.size.*positive integer"
  )
  expect_error(
    do.call(
      SCimplify_by_graph_group,
      c(common, list(min.metacells.per.condition = 1.5))
    ),
    "min.metacells.per.condition.*positive integer"
  )
})
