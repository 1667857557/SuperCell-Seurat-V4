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

test_that("MetacellExpression maps named metacell.names by grouping level", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(c(10, 1), nrow = 1, sparse = TRUE)
  rownames(counts) <- "gene1"
  colnames(counts) <- c("cell1", "cell2")
  object <- Seurat::CreateSeuratObject(counts = counts)
  object$metacell <- c("MC2", "MC1")

  out <- MetacellExpression(
    object,
    assays = "RNA",
    group.by = "metacell",
    metacell.names = stats::setNames(c("MC2", "MC1"), c("MC2", "MC1")),
    return.seurat = FALSE
  )[["RNA"]]

  expect_equal(as.numeric(out["gene1", "MC1"]), 1)
  expect_equal(as.numeric(out["gene1", "MC2"]), 10)
})

test_that("MetacellExpression preserves a single grouping level", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      1:12,
      nrow = 3,
      dimnames = list(
        paste0("g", 1:3),
        paste0("c", 1:4)
      )
    ),
    sparse = TRUE
  )

  object <- Seurat::CreateSeuratObject(counts = counts)
  object$metacell_g20 <- "MC1"

  out <- MetacellExpression(
    object,
    assays = "RNA",
    group.by = "metacell_g20",
    metacell.names = c(MC1 = "MC1"),
    return.seurat = FALSE
  )

  expect_identical(colnames(out$RNA), "MC1")
})
