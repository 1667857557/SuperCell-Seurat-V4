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


test_that("AggregateFragmentFile accepts data.frame membership", {
  skip_if(Sys.which("bgzip") == "")
  skip_if(Sys.which("tabix") == "")

  td <- tempdir()
  fragment_tsv <- file.path(td, "fragments_df.tsv")
  writeLines(c(
    "chr1\t100\t150\tcell1\t1",
    "chr1\t200\t240\tcell2\t1",
    "chr2\t100\t120\tcell3\t1"
  ), fragment_tsv)
  system2("bgzip", c("-f", fragment_tsv))
  input <- paste0(fragment_tsv, ".gz")

  membership_df <- data.frame(
    cell_id = c("cell1", "cell2", "cell3"),
    metacell_id = c("MC1", "MC1", "MC2"),
    stringsAsFactors = FALSE
  )
  out <- SuperCell:::AggregateFragmentFile(
    input_file = input,
    membership = membership_df,
    output_name = "mc_fragments_df.tsv.gz",
    output_path = td,
    nb_cl = 1L,
    returnOutputFileName = TRUE
  )

  expect_true(file.exists(out))
  expect_true(file.exists(paste0(out, ".tbi")))
  got <- data.table::fread(cmd = paste("bgzip -dc", shQuote(out)))
  expect_setequal(unique(got[[4]]), c("MC1", "MC2"))
  expect_false(any(is.na(got[[4]])))
})
