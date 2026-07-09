test_that("assay filtering ignores global .FilterObjects shims", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  old <- if (exists(".FilterObjects", envir = .GlobalEnv, inherits = FALSE)) {
    get(".FilterObjects", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  assign(".FilterObjects", function(...) stop("global shim should not be used"), envir = .GlobalEnv)
  on.exit({
    if (is.null(old)) rm(".FilterObjects", envir = .GlobalEnv) else assign(".FilterObjects", old, envir = .GlobalEnv)
  }, add = TRUE)

  mat <- Matrix::Matrix(c(1, 0, 0, 1), nrow = 2, sparse = TRUE)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("cell1", "cell2")
  obj <- Seurat::CreateSeuratObject(counts = mat)
  expect_true("RNA" %in% SuperCell:::.SCFilterAssays(obj))
})
