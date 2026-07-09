#' Transform fragments file
#'
#' This function transforms a fragment file associating each fragment to a single cell
#' to a fragment file associating fragments to metacells.
#'
#' @param input_file path to original fragments file
#' @param output_name name of the supercell fragments file
#' @param output_path path where the new fragments file will be saved
#' @param membership named vector mapping single-cell barcodes to output metacell barcodes
#' @param n_skip ignored; retained for backward compatibility
#' @return Invisibly returns the output fragment filename when generated.
#' @examples
#' \dontrun{
#' transform_fragment_file("fragments.tsv.gz", "MC_fragments.tsv.gz", "out", membership)
#' }
#' @export
transform_fragment_file <- function(input_file, output_name, output_path, membership, n_skip = 0) {
  AggregateFragmentFile(
    input_file = input_file,
    membership = membership,
    output_name = output_name,
    output_path = output_path,
    returnOutputFileName = TRUE
  )
}

#' Transform fragments file parallel
#'
#' Backward-compatible wrapper around \code{AggregateFragmentFile}.
#'
#' @inheritParams AggregateFragmentFile
#' @param prefixMC prefix to prepend to membership values for legacy callers. New code
#' should pass final metacell barcodes in \code{membership} and leave this empty.
#' @export
transform_fragment_file_parallel <- function(input_file,
                                             membership,
                                             output_name = NULL,
                                             output_path = NULL,
                                             tmp_path = "./tmp/",
                                             nb_cl = NULL,
                                             nb_row_split = 10000000L,
                                             prefixMC = "SC_",
                                             bgzip_path = NULL,
                                             tabix_path = NULL,
                                             returnOutputFileName = TRUE) {
  if (!is.null(prefixMC) && nzchar(prefixMC)) {
    membership_names <- names(membership)
    membership <- paste0(prefixMC, membership)
    names(membership) <- membership_names
  }
  AggregateFragmentFile(
    input_file = input_file,
    membership = membership,
    output_name = output_name,
    output_path = output_path,
    tmp_path = tmp_path,
    nb_cl = nb_cl,
    nb_row_split = nb_row_split,
    bgzip_path = bgzip_path,
    tabix_path = tabix_path,
    returnOutputFileName = returnOutputFileName
  )
}

.SCFragmentCommands <- function(bgzip_path = NULL, tabix_path = NULL) {
  bgzip_command <- bgzip_path %||% Sys.which("bgzip")
  tabix_command <- tabix_path %||% Sys.which("tabix")
  if (!nzchar(bgzip_command)) {
    stop("`bgzip` was not found. Provide `bgzip_path`.", call. = FALSE)
  }
  if (!nzchar(tabix_command)) {
    stop("`tabix` was not found. Provide `tabix_path`.", call. = FALSE)
  }
  list(bgzip = bgzip_command, tabix = tabix_command)
}

.SCFragmentOutputNames <- function(input_file, output_name = NULL, output_path = NULL) {
  if (is.null(output_name)) {
    output_name <- paste0("MC_", fs::path_file(input_file))
  }
  if (is.null(output_path)) {
    output_path <- fs::path_dir(input_file)
  }
  output_path <- fs::path_abs(output_path)
  if (grepl("\\.gz$", output_name)) {
    gz <- file.path(output_path, output_name)
    tsv <- sub("\\.gz$", "", gz)
  } else if (grepl("\\.tsv$", output_name)) {
    tsv <- file.path(output_path, output_name)
    gz <- paste0(tsv, ".gz")
  } else {
    tsv <- file.path(output_path, paste0(output_name, ".tsv"))
    gz <- paste0(tsv, ".gz")
  }
  list(output_path = output_path, gz = gz, tsv = tsv)
}

.SCCopyFilesBinary <- function(files, output_file) {
  if (length(files) == 0L) {
    stop("No updated fragment chunks were produced.", call. = FALSE)
  }
  con_out <- file(output_file, open = "wb")
  close_out <- TRUE
  on.exit(if (close_out) close(con_out), add = TRUE)
  for (f in files) {
    con_in <- file(f, open = "rb")
    on.exit(close(con_in), add = TRUE)
    repeat {
      buf <- readBin(con_in, what = "raw", n = 1024 * 1024)
      if (!length(buf)) break
      writeBin(buf, con_out)
    }
    close(con_in)
  }
  close(con_out)
  close_out <- FALSE
  if (!file.exists(output_file) || file.info(output_file)$size == 0) {
    stop("Fragment concatenation produced an empty output file.", call. = FALSE)
  }
  invisible(output_file)
}

.SCValidateFragmentTable <- function(fragment_table, expected_cells) {
  preview <- data.table::fread(fragment_table, nrows = 10000)
  if (nrow(preview) == 0L || ncol(preview) < 4L) {
    stop("Aggregated fragment table is empty or malformed.", call. = FALSE)
  }
  found <- unique(as.character(preview[[4]]))
  if (any(is.na(found)) || any(found == "NA")) {
    stop("Aggregated fragment table contains NA barcodes.", call. = FALSE)
  }
  expected_cells <- unique(as.character(expected_cells))
  if (length(intersect(expected_cells, found)) == 0L) {
    stop("None of the expected metacell barcodes were found in aggregated fragment table.", call. = FALSE)
  }
  unexpected <- setdiff(found, expected_cells)
  if (length(unexpected) > 0L) {
    stop("Aggregated fragment table contains unexpected barcodes: ",
         paste(utils::head(unexpected, 10), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.SCValidateAggregatedFragments <- function(fragment_file, expected_cells, tabix_path = "tabix") {
  if (!file.exists(fragment_file)) {
    stop("Missing fragment file: ", fragment_file, call. = FALSE)
  }
  if (file.info(fragment_file)$size == 0) {
    stop("Aggregated fragment file is empty: ", fragment_file, call. = FALSE)
  }
  index_file <- paste0(fragment_file, ".tbi")
  if (!file.exists(index_file)) {
    stop("Missing tabix index: ", index_file, call. = FALSE)
  }
  expected_cells <- unique(as.character(expected_cells))
  preview <- tryCatch(
    data.table::fread(cmd = paste(shQuote(tabix_path), "-H", shQuote(fragment_file)), nrows = 1000),
    error = function(e) NULL
  )
  if (!is.null(preview) && ncol(preview) >= 4L && nrow(preview) > 0L) {
    found <- unique(as.character(preview[[4]]))
    if (any(is.na(found)) || any(found == "NA")) {
      stop("Aggregated fragment file contains NA barcodes.", call. = FALSE)
    }
    if (length(intersect(expected_cells, found)) == 0L) {
      stop("None of the expected metacell barcodes were found in aggregated fragment file.", call. = FALSE)
    }
    unexpected <- setdiff(found, expected_cells)
    if (length(unexpected) > 0L) {
      stop("Aggregated fragment file contains unexpected barcodes: ",
           paste(utils::head(unexpected, 10), collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' AggregateFragmentFile parallel
#'
#' Transform fragments from single-cell barcodes to metacell barcodes.
#'
#' @param input_file path to original fragments file
#' @param membership named vector mapping single-cell barcodes to final metacell barcodes
#' @param output_name name of the metacell fragments file
#' @param output_path directory where the new fragments file will be saved
#' @param tmp_path directory path for temporary files
#' @param nb_cl number of workers to use
#' @param nb_row_split row number for split files
#' @param bgzip_path optional path to bgzip
#' @param tabix_path optional path to tabix
#' @param returnOutputFileName whether to return the new fragment file name
#' @param keep_unmatched whether to keep fragments whose barcode is absent from membership
#' @importFrom foreach %dopar%
#' @export
AggregateFragmentFile <- function(input_file,
                                  membership,
                                  output_name = NULL,
                                  output_path = NULL,
                                  tmp_path = "./tmp/",
                                  nb_cl = NULL,
                                  nb_row_split = 10000000L,
                                  bgzip_path = NULL,
                                  tabix_path = NULL,
                                  returnOutputFileName = TRUE,
                                  keep_unmatched = FALSE) {
  if (is.null(names(membership)) || any(!nzchar(names(membership)))) {
    stop("`membership` must be a named vector with single-cell barcodes as names.", call. = FALSE)
  }
  membership <- as.character(membership)
  names(membership) <- names(membership)
  commands <- .SCFragmentCommands(bgzip_path = bgzip_path, tabix_path = tabix_path)
  outputs <- .SCFragmentOutputNames(input_file, output_name, output_path)
  output_path <- outputs$output_path
  full_output_name <- outputs$gz
  full_output_name_tsv <- outputs$tsv
  tmp_path <- fs::path_abs(tmp_path)

  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(tmp_path, showWarnings = FALSE, recursive = TRUE)

  unique_prefix <- paste0("frags_subset_", Sys.getpid(), "_", as.integer(Sys.time()), "_", sample.int(1e6, 1))
  split_prefix <- file.path(tmp_path, unique_prefix)
  decompressed <- file.path(tmp_path, paste0(unique_prefix, "_fragments.tsv"))
  cleanup_pattern <- paste0("^", basename(unique_prefix))
  on.exit(unlink(list.files(path = tmp_path, pattern = cleanup_pattern, full.names = TRUE), force = TRUE), add = TRUE)

  message("Start fragment decompression")
  decompress_status <- system2(commands$bgzip, args = c("-dc", input_file), stdout = decompressed)
  if (!identical(decompress_status, 0L) || !file.exists(decompressed) || file.info(decompressed)$size == 0) {
    stop("Fragment decompression produced empty file.", call. = FALSE)
  }

  message("Start fragments file split")
  split_status <- system2("split", args = c("-l", as.character(as.integer(nb_row_split)), decompressed, split_prefix))
  if (!identical(split_status, 0L)) {
    stop("Fragment file split failed.", call. = FALSE)
  }
  list_fragments <- list.files(path = tmp_path, pattern = paste0("^", basename(unique_prefix), "[a-z]+$"), full.names = TRUE)
  if (length(list_fragments) == 0L) {
    stop("Fragment file split produced no chunks.", call. = FALSE)
  }

  matching_names <- data.frame(
    single_cell_names = names(membership),
    super_cell_names = membership,
    stringsAsFactors = FALSE
  )
  rownames(matching_names) <- matching_names$single_cell_names

  nb_cl <- .SCResolveNbCl(nb_cl)
  update_one <- function(fragment) {
    tmp <- data.table::fread(fragment)
    if (!keep_unmatched) {
      tmp <- tmp[tmp$V4 %in% rownames(matching_names), ]
    }
    tmp$V4 <- matching_names[tmp$V4, "super_cell_names"]
    tmp <- tmp[!is.na(tmp$V4) & tmp$V4 != "NA", ]
    out <- paste0(fragment, "_up")
    data.table::fwrite(tmp, out, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    out
  }

  message("Start fragment barcode update")
  if (nb_cl > 1L && length(list_fragments) > 1L) {
    cl <- parallel::makeCluster(nb_cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    up_files <- foreach::foreach(fragment = list_fragments, .packages = c("data.table")) %dopar% {
      update_one(fragment)
    }
  } else {
    up_files <- lapply(list_fragments, update_one)
  }
  up_files <- unlist(up_files, use.names = FALSE)

  message("Start fragment concatenation")
  .SCCopyFilesBinary(up_files, full_output_name_tsv)

  .SCValidateFragmentTable(full_output_name_tsv, expected_cells = unique(membership))

  message("Start bgzip")
  bgzip_status <- system2(commands$bgzip, args = c("-f", full_output_name_tsv))
  if (!identical(bgzip_status, 0L) || !file.exists(full_output_name) || file.info(full_output_name)$size == 0) {
    stop("bgzip failed or produced an empty output file.", call. = FALSE)
  }

  message("Start tabix")
  tabix_status <- system2(commands$tabix, args = c("-f", "-p", "bed", full_output_name))
  if (!identical(tabix_status, 0L)) {
    stop("tabix indexing failed.", call. = FALSE)
  }
  .SCValidateAggregatedFragments(full_output_name, expected_cells = unique(membership), tabix_path = commands$tabix)

  if (returnOutputFileName) {
    return(full_output_name)
  }
  invisible(full_output_name)
}
