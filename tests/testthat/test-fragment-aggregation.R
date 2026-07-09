test_that("AggregateFragmentFile filters unmatched barcodes", {
  skip_if(Sys.which("bgzip") == "")
  skip_if(Sys.which("tabix") == "")

  td <- tempdir()
  fragment_tsv <- file.path(td, "fragments.tsv")
  writeLines(c(
    "chr1\t100\t200\tAAAC-1\t1",
    "chr1\t300\t400\tAAAG-1\t1",
    "chr1\t500\t600\tUNMATCHED-1\t1"
  ), fragment_tsv)
  system2("bgzip", c("-f", fragment_tsv))
  input <- paste0(fragment_tsv, ".gz")

  out <- SuperCell:::AggregateFragmentFile(
    input_file = input,
    membership = c("AAAC-1" = "Metacell_1", "AAAG-1" = "Metacell_2"),
    output_name = "mc_fragments.tsv.gz",
    output_path = td,
    tmp_path = file.path(td, "tmp"),
    nb_cl = 1L,
    returnOutputFileName = TRUE
  )
  expect_true(file.exists(out))
  expect_true(file.exists(paste0(out, ".tbi")))
  got <- data.table::fread(cmd = paste("bgzip -dc", shQuote(out)))
  expect_setequal(unique(got[[4]]), c("Metacell_1", "Metacell_2"))
  expect_false(any(is.na(got[[4]])))
  expect_false("UNMATCHED-1" %in% got[[4]])
})
