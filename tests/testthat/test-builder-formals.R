test_that("SuperCell2 builders do not expose sample_col", {
  builders <- c(
    "SCimplify",
    "SCimplify_from_embedding",
    "SCimplify_for_Seurat",
    "SCimplify_for_velocity"
  )

  builder_formals <- lapply(builders, function(builder) {
    names(formals(getExportedValue("SuperCell", builder)))
  })
  names(builder_formals) <- builders

  expect_false(any(vapply(
    builder_formals,
    function(arguments) "sample_col" %in% arguments,
    logical(1)
  )))
})
