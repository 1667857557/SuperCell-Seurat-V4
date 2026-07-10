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
                                             tmp_path = NULL,
                                             nb_cl = NULL,
                                             nb_row_split = 10000000L,
                                             prefixMC = "SC_",
                                             bgzip_path = NULL,
                                             tabix_path = NULL,
                                             returnOutputFileName = TRUE) {
  membership <- .SCNormalizeMembership(membership)
  if (!is.null(prefixMC) && nzchar(prefixMC)) {
    cell_ids <- names(membership)
    membership <- stats::setNames(paste0(prefixMC, unname(membership)), cell_ids)
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

.SCNormalizeMembership <- function(membership) {
  if (is.data.frame(membership)) {
    if (all(c("cell_id", "metacell_id") %in% colnames(membership))) {
      cell_id <- membership$cell_id
      metacell_id <- membership$metacell_id
    } else if (ncol(membership) >= 2L) {
      cell_id <- membership[[1L]]
      metacell_id <- membership[[2L]]
    } else {
      stop("Membership data.frame must contain cell_id and metacell_id.", call. = FALSE)
    }
    membership <- stats::setNames(as.character(metacell_id), as.character(cell_id))
  } else {
    cell_id <- names(membership)
    if (is.null(cell_id) || length(cell_id) != length(membership) || anyNA(cell_id) || any(!nzchar(cell_id))) {
      stop("`membership` must be a named vector or a data.frame with cell_id and metacell_id.", call. = FALSE)
    }
    membership <- stats::setNames(as.character(unname(membership)), as.character(cell_id))
  }
  if (!length(membership)) stop("Membership is empty.", call. = FALSE)
  if (anyDuplicated(names(membership))) {
    duplicated_ids <- unique(names(membership)[duplicated(names(membership))])
    stop("Duplicated single-cell barcodes in membership: ", paste(utils::head(duplicated_ids, 10L), collapse = ", "), call. = FALSE)
  }
  if (anyNA(membership) || any(!nzchar(membership))) {
    stop("Membership contains missing or empty metacell IDs.", call. = FALSE)
  }
  membership
}

.SCBuildFragmentBarcodeMap <- function(membership) {
  direct <- data.frame(
    fragment_barcode = names(membership),
    super_cell_names = unname(membership),
    stringsAsFactors = FALSE
  )

  stripped <- sub("^[^_]+_", "", names(membership))
  can_strip <- stripped != names(membership) & nzchar(stripped)
  unique_stripped <- can_strip & !duplicated(stripped) & !duplicated(stripped, fromLast = TRUE)
  stripped <- stripped[unique_stripped]

  stripped_map <- data.frame(
    fragment_barcode = stripped,
    super_cell_names = unname(membership)[unique_stripped],
    stringsAsFactors = FALSE
  )
  stripped_map <- stripped_map[!(stripped_map$fragment_barcode %in% direct$fragment_barcode), , drop = FALSE]

  matching_names <- rbind(direct, stripped_map)
  rownames(matching_names) <- matching_names$fragment_barcode
  matching_names
}

.SCFragmentCommands <- function(bgzip_path = NULL, tabix_path = NULL, split_path = NULL) {
  bgzip_command <- bgzip_path %||% Sys.which("bgzip")
  tabix_command <- tabix_path %||% Sys.which("tabix")
  split_command <- split_path %||% Sys.which("split")
  if (!nzchar(bgzip_command)) {
    stop("`bgzip` was not found. Provide `bgzip_path`.", call. = FALSE)
  }
  if (!nzchar(tabix_command)) {
    stop("`tabix` was not found. Provide `tabix_path`.", call. = FALSE)
  }
  if (!nzchar(split_command)) {
    stop("GNU/coreutils `split` was not found.", call. = FALSE)
  }
  list(bgzip = bgzip_command, tabix = tabix_command, split = split_command)
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
  preview <- data.table::fread(fragment_table, nrows = 10000, header = FALSE)
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

.SCValidateAggregatedFragments <- function(fragment_file, expected_cells, bgzip_path = "bgzip") {
  if (!file.exists(fragment_file)) {
    stop("Missing fragment file: ", fragment_file, call. = FALSE)
  }
  if (file.info(fragment_file)$size == 0) {
    stop("Aggregated fragment file is empty: ", fragment_file, call. = FALSE)
  }
  index_file <- paste0(fragment_file, ".tbi")
  if (!file.exists(index_file) || file.info(index_file)$size == 0) {
    stop("Missing or empty tabix index: ", index_file, call. = FALSE)
  }
  cmd <- paste(shQuote(bgzip_path), "-dc", shQuote(fragment_file), "| head -n 10000")
  preview <- tryCatch(data.table::fread(cmd = cmd, header = FALSE, showProgress = FALSE), error = function(e) e)
  if (inherits(preview, "error")) {
    stop("Failed to read aggregated fragment preview: ", conditionMessage(preview), call. = FALSE)
  }
  if (nrow(preview) == 0L || ncol(preview) < 4L) {
    stop("Aggregated fragment file is empty or malformed.", call. = FALSE)
  }
  found <- unique(as.character(preview[[4L]]))
  expected_cells <- unique(as.character(expected_cells))
  if (anyNA(found) || any(!nzchar(found)) || any(found == "NA")) {
    stop("Aggregated fragment file contains missing barcodes.", call. = FALSE)
  }
  unexpected <- setdiff(found, expected_cells)
  if (length(unexpected)) {
    stop("Aggregated fragment file contains unexpected metacell barcodes: ", paste(utils::head(unexpected, 10L), collapse = ", "), call. = FALSE)
  }
  if (!length(intersect(found, expected_cells))) {
    stop("None of the expected metacell barcodes were found.", call. = FALSE)
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
                                  tmp_path = NULL,
                                  nb_cl = NULL,
                                  nb_row_split = 10000000L,
                                  bgzip_path = NULL,
                                  tabix_path = NULL,
                                  returnOutputFileName = TRUE,
                                  keep_unmatched = FALSE) {
  membership <- .SCNormalizeMembership(membership)
  commands <- .SCFragmentCommands(bgzip_path = bgzip_path, tabix_path = tabix_path)
  if (!file.exists(input_file)) stop("Fragment file not found: ", input_file, call. = FALSE)
  if (file.access(input_file, mode = 4L) != 0L) stop("Fragment file is not readable: ", input_file, call. = FALSE)
  outputs <- .SCFragmentOutputNames(input_file, output_name, output_path)
  output_path <- outputs$output_path
  full_output_name <- outputs$gz
  full_output_name_tsv <- outputs$tsv
  if (is.null(tmp_path)) {
    tmp_path <- tempfile(pattern = paste0("SuperCell_fragments_", Sys.getpid(), "_"))
  }
  tmp_path <- fs::path_abs(tmp_path)

  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  if (file.access(output_path, mode = 2L) != 0L) stop("Output directory is not writable: ", output_path, call. = FALSE)
  dir.create(tmp_path, showWarnings = FALSE, recursive = TRUE)

  unique_prefix <- paste0("frags_subset_", Sys.getpid(), "_", as.integer(Sys.time()), "_", sample.int(1e6, 1))
  split_prefix <- file.path(tmp_path, unique_prefix)
  decompressed <- file.path(tmp_path, paste0(unique_prefix, "_fragments.tsv"))
  cleanup_pattern <- paste0("^", basename(unique_prefix))
  on.exit(unlink(list.files(path = tmp_path, pattern = cleanup_pattern, full.names = TRUE), force = TRUE), add = TRUE)

  preview_cmd <- paste(shQuote(commands$bgzip), "-dc", shQuote(input_file), "| head -n 100000")
  preview <- data.table::fread(cmd = preview_cmd, header = FALSE, showProgress = FALSE)
  if (nrow(preview) == 0L || ncol(preview) < 4L) stop("Input fragment file is empty or malformed.", call. = FALSE)
  preview_barcodes <- unique(as.character(preview[[4L]]))
  preview_match <- mean(preview_barcodes %in% names(membership))
  if (!is.finite(preview_match) || preview_match == 0) {
    stop("No fragment barcodes matched the supplied membership. Check barcode prefixes and sample mapping.", call. = FALSE)
  }

  message("Start fragment decompression")
  decompress_status <- system2(commands$bgzip, args = c("-dc", input_file), stdout = decompressed)
  if (!identical(decompress_status, 0L) || !file.exists(decompressed) || file.info(decompressed)$size == 0) {
    stop("Fragment decompression produced empty file.", call. = FALSE)
  }

  message("Start fragments file split")
  split_status <- system2(commands$split, args = c("-l", as.character(as.integer(nb_row_split)), decompressed, split_prefix))
  if (!identical(split_status, 0L)) {
    stop("Fragment file split failed.", call. = FALSE)
  }
  list_fragments <- list.files(path = tmp_path, pattern = paste0("^", basename(unique_prefix), "[a-z]+$"), full.names = TRUE)
  if (length(list_fragments) == 0L) {
    stop("Fragment file split produced no chunks.", call. = FALSE)
  }

  matching_names <- .SCBuildFragmentBarcodeMap(membership)

  nb_cl <- .SCResolveNbCl(nb_cl)
  update_one <- function(fragment) {
    tmp <- data.table::fread(fragment, header = FALSE)
    if (nrow(tmp) == 0L || ncol(tmp) < 4L) {
      stop("Input fragment file is empty or malformed.", call. = FALSE)
    }
    if (!keep_unmatched) {
      tmp <- tmp[tmp$V4 %in% rownames(matching_names), ]
    }
    tmp$V4 <- matching_names[tmp$V4, "super_cell_names"]
    tmp <- tmp[!is.na(tmp$V4) & tmp$V4 != "NA", ]
    out <- paste0(fragment, "_up")
    data.table::fwrite(tmp, out, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    list(file = out, n = nrow(tmp))
  }

  message("Start fragment barcode update")
  if (nb_cl > 1L && length(list_fragments) > 1L) {
    cl <- parallel::makeCluster(nb_cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    update_results <- foreach::foreach(fragment = list_fragments, .packages = c("data.table")) %dopar% {
      update_one(fragment)
    }
  } else {
    update_results <- lapply(list_fragments, update_one)
  }
  matched_rows <- sum(vapply(update_results, function(x) x$n, integer(1)))
  if (matched_rows == 0L) {
    stop("No fragment barcodes matched the supplied membership. Check barcode prefixes and sample mapping.", call. = FALSE)
  }
  up_files <- vapply(update_results, function(x) x$file, character(1))

  message("Start fragment concatenation")
  .SCCopyFilesBinary(up_files, full_output_name_tsv)

  .SCValidateFragmentTable(full_output_name_tsv, expected_cells = unique(membership))

  message("Start bgzip")
  bgzip_status <- system2(commands$bgzip, args = c("-f", full_output_name_tsv))
  if (!identical(bgzip_status, 0L) || !file.exists(full_output_name) || file.info(full_output_name)$size == 0) {
    stop("bgzip failed or produced an empty output file.", call. = FALSE)
  }

  message("Start tabix")
  sort_check <- system2(commands$tabix, args = c("-f", "-p", "bed", full_output_name), stdout = TRUE, stderr = TRUE)
  tabix_status <- attr(sort_check, "status") %||% 0L
  if (!identical(tabix_status, 0L)) {
    stop("tabix indexing failed. The aggregated fragment file may be unsorted or malformed: ", paste(sort_check, collapse = "\n"), call. = FALSE)
  }
  .SCValidateAggregatedFragments(full_output_name, expected_cells = unique(membership), bgzip_path = commands$bgzip)

  if (returnOutputFileName) {
    return(full_output_name)
  }
  invisible(full_output_name)
}
