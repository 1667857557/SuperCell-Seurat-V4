test_that("named vector membership preserves cell IDs", {
  x <- c(cell1 = "MC1", cell2 = "MC1", cell3 = "MC2")
  y <- SuperCell:::.SCNormalizeMembership(x)

  expect_identical(names(y), names(x))
  expect_identical(unname(y), unname(x))
})

test_that("membership data.frame is accepted", {
  x <- data.frame(
    cell_id = c("cell1", "cell2"),
    metacell_id = c("MC1", "MC1")
  )

  y <- SuperCell:::.SCNormalizeMembership(x)

  expect_identical(names(y), x$cell_id)
  expect_identical(unname(y), x$metacell_id)
})

test_that("duplicated cell IDs are rejected", {
  x <- c(cell1 = "MC1", cell1 = "MC2")
  expect_error(SuperCell:::.SCNormalizeMembership(x), "Duplicated")
})

test_that("fragment barcode map supports unique sample-prefixed cell IDs", {
  membership <- c(Pool1_cell1 = "MC1", Pool1_cell2 = "MC2")
  barcode_map <- SuperCell:::.SCBuildFragmentBarcodeMap(SuperCell:::.SCNormalizeMembership(membership))

  expect_identical(barcode_map["cell1", "super_cell_names"], "MC1")
  expect_identical(barcode_map["cell2", "super_cell_names"], "MC2")
})

test_that("fragment barcode map does not add ambiguous stripped IDs", {
  membership <- c(Pool1_cell1 = "MC1", Pool2_cell1 = "MC2")
  barcode_map <- SuperCell:::.SCBuildFragmentBarcodeMap(SuperCell:::.SCNormalizeMembership(membership))

  expect_false("cell1" %in% rownames(barcode_map))
})
