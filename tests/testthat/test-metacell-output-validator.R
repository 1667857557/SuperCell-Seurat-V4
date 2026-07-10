test_that("validate_metacell_output validates membership contract", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(c(1, 0, 2, 0, 3, 1), nrow = 2, sparse = TRUE)
  rownames(counts) <- c("gene1", "gene2")
  colnames(counts) <- c("MC1", "MC2", "MC3")
  object <- Seurat::CreateSeuratObject(counts = counts)
  object@misc$membership_table <- data.frame(
    cell_id = c("cell1", "cell2", "cell3"),
    metacell_id = c("MC1", "MC2", "MC3"),
    stringsAsFactors = FALSE
  )
  object@misc$fragment_manifest <- data.frame(
    assay = character(),
    input_file = character(),
    fragment_file = character(),
    index_file = character(),
    n_input_membership_cells = integer(),
    n_metacells = integer(),
    status = character(),
    stringsAsFactors = FALSE
  )

  expect_true(isTRUE(SuperCell::validate_metacell_output(object)))
})
