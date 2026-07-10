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


test_that("validate_metacell_output rejects malformed fragment manifests", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(c(1, 2), nrow = 1, sparse = TRUE)
  rownames(counts) <- "gene1"
  colnames(counts) <- c("MC1", "MC2")
  object <- Seurat::CreateSeuratObject(counts = counts)
  object@misc$membership_table <- data.frame(
    cell_id = c("cell1", "cell2"),
    metacell_id = c("MC1", "MC2"),
    stringsAsFactors = FALSE
  )
  object@misc$fragment_manifest <- data.frame(
    assay = "ATAC",
    status = "ok",
    stringsAsFactors = FALSE
  )

  expect_error(
    SuperCell::validate_metacell_output(object, require_fragments = TRUE),
    "required columns"
  )
})

test_that("validate_metacell_output validates fragment cell maps and manifest counts", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(c(1, 0, 2, 0), nrow = 2, sparse = TRUE)
  rownames(counts) <- c("gene1", "gene2")
  colnames(counts) <- c("MC1", "MC2")
  object <- Seurat::CreateSeuratObject(counts = counts)
  object@misc$membership_table <- data.frame(
    cell_id = c("cell1", "cell2"),
    metacell_id = c("MC1", "MC2"),
    stringsAsFactors = FALSE
  )
  td <- tempfile("fragment_manifest_")
  dir.create(td)
  fragment_file <- file.path(td, "fragments.tsv.gz")
  index_file <- paste0(fragment_file, ".tbi")
  writeBin(as.raw(1:3), fragment_file)
  writeBin(as.raw(1:3), index_file)
  object@misc$fragment_manifest <- data.frame(
    assay = "ATAC",
    input_file = "input.tsv.gz",
    fragment_file = fragment_file,
    index_file = index_file,
    n_input_membership_cells = 2L,
    n_metacells = 2L,
    n_fragment_rows = 2L,
    status = "ok",
    stringsAsFactors = FALSE
  )
  object@misc$fragment_cell_map <- data.frame(
    assay = "ATAC",
    fragment_file = fragment_file,
    object_cell = c("MC1", "MC2"),
    fragment_barcode = c("MC1", "MC2"),
    stringsAsFactors = FALSE
  )

  expect_true(isTRUE(SuperCell::validate_metacell_output(object, require_fragments = TRUE,
                                                         require_complete_fragment_coverage = TRUE)))
})

test_that("validate_metacell_output rejects duplicated fragment object cells", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(c(1, 2), nrow = 1, sparse = TRUE)
  rownames(counts) <- "gene1"
  colnames(counts) <- c("MC1", "MC2")
  object <- Seurat::CreateSeuratObject(counts = counts)
  object@misc$membership_table <- data.frame(
    cell_id = c("cell1", "cell2"),
    metacell_id = c("MC1", "MC2"),
    stringsAsFactors = FALSE
  )
  td <- tempfile("fragment_manifest_dup_")
  dir.create(td)
  fragment_files <- file.path(td, c("frag1.tsv.gz", "frag2.tsv.gz"))
  index_files <- paste0(fragment_files, ".tbi")
  for (path in c(fragment_files, index_files)) writeBin(as.raw(1:3), path)
  object@misc$fragment_manifest <- data.frame(
    fragment_file = fragment_files,
    index_file = index_files,
    n_metacells = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  object@misc$fragment_cell_map <- data.frame(
    fragment_file = fragment_files,
    object_cell = c("MC1", "MC1"),
    fragment_barcode = c("MC1", "MC1"),
    stringsAsFactors = FALSE
  )

  expect_error(
    SuperCell::validate_metacell_output(object, require_fragments = TRUE),
    "assigned to multiple fragment files"
  )
})
