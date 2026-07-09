test_that("MetacellExpression preserves grouping names", {
  skip_if_not_installed("Seurat")
  mat <- Matrix::Matrix(c(1, 0, 0, 1, 2, 0), nrow = 2, sparse = TRUE)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("cell1", "cell2", "cell3")
  obj <- Seurat::CreateSeuratObject(counts = mat)
  obj$metacell_id <- c("MC_B", "MC_A", "MC_B")

  agg <- SuperCell:::MetacellExpression(
    obj,
    assays = "RNA",
    group.by = "metacell_id",
    return.seurat = FALSE
  )[["RNA"]]
  expect_identical(colnames(agg), sort(unique(obj$metacell_id)))

  agg_named <- SuperCell:::MetacellExpression(
    obj,
    assays = "RNA",
    group.by = "metacell_id",
    metacell.names = c("A", "B"),
    return.seurat = FALSE
  )[["RNA"]]
  expect_identical(colnames(agg_named), c("A", "B"))
})
