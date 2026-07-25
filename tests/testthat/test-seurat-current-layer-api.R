test_that("compatibility helpers read and write v3 assays", {
  old_option <- getOption("Seurat.object.assay.version")
  on.exit(options(Seurat.object.assay.version = old_option), add = TRUE)
  options(Seurat.object.assay.version = "v3")

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 0, 1, 4, 0), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("g1", "g2"), paste0("cell", 1:4))
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  expect_false(inherits(object[["RNA"]], "Assay5"))
  expect_equal(
    as.matrix(SuperCell:::.sc_get_assay_data(
      object, assay = "RNA", slot = "counts"
    )),
    as.matrix(counts)
  )

  normalized <- log1p(counts)
  object <- SuperCell:::.sc_set_assay_data(
    object,
    assay = "RNA",
    slot = "data",
    new.data = normalized
  )
  expect_equal(
    as.matrix(SuperCell:::.sc_get_assay_data(
      object, assay = "RNA", slot = "data"
    )),
    as.matrix(normalized)
  )
  fetched <- SuperCell:::.sc_fetch_data(
    object,
    vars = "g1",
    slot = "data"
  )
  expect_equal(as.numeric(fetched$g1), as.numeric(normalized["g1", ]))
})

test_that("compatibility helpers read Signac ChromatinAssay", {
  skip_if_not_installed("Signac")
  old_option <- getOption("Seurat.object.assay.version")
  on.exit(options(Seurat.object.assay.version = old_option), add = TRUE)
  options(Seurat.object.assay.version = "v3")

  rna <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 0, 1, 4, 0), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("g1", "g2"), paste0("cell", 1:4))
  )
  peaks <- Matrix::Matrix(
    matrix(c(1, 0, 0, 2, 0, 3, 1, 0), nrow = 2),
    sparse = TRUE,
    dimnames = list(c("chr1-1-10", "chr1-20-30"), colnames(rna))
  )
  object <- Seurat::CreateSeuratObject(counts = rna)
  object[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = peaks,
    sep = c("-", "-")
  )
  expect_true(inherits(object[["ATAC"]], "ChromatinAssay"))
  expect_equal(
    as.matrix(SuperCell:::.sc_get_assay_data(
      object, assay = "ATAC", slot = "counts"
    )),
    as.matrix(peaks)
  )
  expect_true(SuperCell:::.sc_assay_has_slot(
    object, assay = "ATAC", slot = "counts"
  ))
})
