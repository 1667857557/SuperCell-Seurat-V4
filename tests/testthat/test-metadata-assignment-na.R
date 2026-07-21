test_that("metadata assignment retains metacells with only missing labels", {
  membership <- c("MC1", "MC1", "MC2", "MC2")
  labels <- c(NA, NA, "B", "B")

  assignment <- supercell_assign(
    clusters = labels,
    supercell_membership = membership,
    method = "absolute"
  )
  purity <- supercell_purity(
    clusters = labels,
    supercell_membership = membership,
    method = "max_proportion"
  )

  expect_identical(names(assignment), c("MC1", "MC2"))
  expect_true(is.na(assignment[["MC1"]]))
  expect_identical(assignment[["MC2"]], "B")
  expect_identical(names(purity), c("MC1", "MC2"))
  expect_true(is.na(purity[["MC1"]]))
  expect_equal(purity[["MC2"]], 1)
})

test_that("metadata assignment uses observed cells when labels are partly missing", {
  membership <- c("MC1", "MC1", "MC1", "MC2", "MC2")
  labels <- c(NA, "A", "A", "B", NA)

  assignment <- supercell_assign(
    clusters = labels,
    supercell_membership = membership,
    method = "absolute"
  )
  purity <- supercell_purity(
    clusters = labels,
    supercell_membership = membership,
    method = "max_proportion"
  )

  expect_identical(unname(assignment), c("A", "B"))
  expect_equal(unname(purity), c(1, 1))
})

test_that("metadata helpers reject unequal input lengths", {
  expect_error(
    supercell_assign(c("A", "B"), "MC1"),
    "same length"
  )
  expect_error(
    supercell_purity(c("A", "B"), "MC1"),
    "same length"
  )
})
