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

test_that("condition and cell-type separation retain their native inputs", {
  matrix_builders <- list(
    SCimplify = SCimplify,
    SCimplify_from_embedding = SCimplify_from_embedding
  )

  for (builder in matrix_builders) {
    builder_arguments <- names(formals(builder))
    expect_true("cell.split.condition" %in% builder_arguments)
    expect_true("cell.annotation" %in% builder_arguments)
  }

  seurat_builder_arguments <- names(formals(SCimplify_for_Seurat))
  expect_true("label" %in% seurat_builder_arguments)
  expect_false("condition_col" %in% seurat_builder_arguments)
  expect_false("sample_col" %in% seurat_builder_arguments)
})
