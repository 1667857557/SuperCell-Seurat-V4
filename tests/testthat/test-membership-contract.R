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
