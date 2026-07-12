test_that("target metacell count is clamped to at least one", {
  expect_warning(
    target <- SuperCell:::.SCResolveTargetMetacells(n_cells = 4, gamma = 20),
    "one metacell"
  )
  expect_identical(target, 1L)
})

test_that("target metacell count validates gamma and cell count", {
  expect_error(
    SuperCell:::.SCResolveTargetMetacells(n_cells = 4, gamma = 0),
    "`gamma` must be a positive finite number"
  )
  expect_error(
    SuperCell:::.SCResolveTargetMetacells(n_cells = 1, gamma = 20),
    "At least two cells are required"
  )
})
