.SCResolveTargetMetacells <- function(n_cells, gamma) {
  if (n_cells < 2L) {
    stop("At least two cells are required to identify metacells.", call. = FALSE)
  }
  gamma <- suppressWarnings(as.numeric(gamma[1L]))
  if (!is.finite(gamma) || gamma <= 0) {
    stop("`gamma` must be a positive finite number.", call. = FALSE)
  }
  n_target <- floor(n_cells / gamma)
  n_target <- max(1L, min(as.integer(n_target), as.integer(n_cells)))
  if (n_target == 1L) {
    warning("The requested gamma produces one metacell for ", n_cells,
            " cells. This is valid for aggregation but provides no within-stratum variability.",
            call. = FALSE)
  } else if (n_target < 5L) {
    warning("Very few metacells will be generated; downstream correlation or differential analysis may be unstable.",
            call. = FALSE)
  }
  n_target
}

#' SCimplify_for_Seurat
#'
#' \code{SCimplify_for_Seurat}
#' Build metacells from a Seurat single-cell object.
#' @param seurat A Seurat single-cell object. It has to be preprocessed (eg. latent space computed) for the assay(s) used to identify metacells
#' @param sobj.mc A metacell seurat object that will be rescaled (optional) it requires a metacell_hierarchy object in the slot misc.
#' @param gamma graining level.
#' @param assay a list of one or two assays to use to build the knn graph on which metacell are identified.
#' @param reduction a list of corresponding reduction name in the seurat single cell object.
#' @param dims a list of corresponding dimensions to use.
#' @param membership a vector of metacell membership as in SuperCell v1 to use to directly aggregate the data (optionnal).
#' @param label optional metadata column used to keep labeled groups separate
#' during metacell construction. This may contain condition, sample, cell-type,
#' or another categorical annotation; `NA` values enable partial annotation.
#' @details `sample_col` and `condition_col` are not formal arguments of this
#' Seurat builder. Pass a Seurat metadata column name through `label` when
#' metacells must not mix known conditions, samples, cell types, or other
#' labeled groups. The matrix builders [SCimplify()] and
#' [SCimplify_from_embedding()] instead accept condition and cell-type vectors
#' through `cell.split.condition` and `cell.annotation`, respectively.
#' @return  A Seurat metacell object with metacell_hierarchy and memberships in the slot misc.
#' @examples
#' sobj.mc <- SCimplify_for_Seurat(seurat = pbmc,
#'                          gamma = 30)
#' @import Seurat
#' @export

SCimplify_for_Seurat <- function(seurat,
                                 seurat.mc = NULL,
                                 k.knn = 30,
                                 kith = NULL,
                                 kernel = T,
                                 gamma = 20,
                                 graph.name = NULL,
                                 assay = c("RNA"),
                                 reduction = list("pca"),
                                 dims = list(c(1:30)),
                                 membership = NULL,
                                 metacellNormalization = F,
                                 avg.in.data = F,
                                 fragmentFiles = NULL,
                                 tmpPath = NULL,
                                 outputDirMcFragment = NULL,
                                 bgzip_path = NULL,
                                 tabix_path = NULL,
                                 prefixMC = "",
                                 seed = 12345L,
                                 return_membership_table = TRUE,
                                 return_fragment_manifest = TRUE,
                                 peakSep = c("-", "-"),
                                 label = NULL,
                                 return.seurat = T,
                                 nb_cl = NULL,
                                 verbose = FALSE)
{
  library(Signac)
  library(Seurat)
  seed <- as.integer(seed)[1L]
  if (!is.finite(seed)) stop("`seed` must be a finite integer.", call. = FALSE)
  old_random_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_random_seed <- if (old_random_seed_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (old_random_seed_exists) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  fragment_manifest_rows <- list()
  fragment_cell_map_rows <- list()

  if (!is.null(label)) {
    seurat[[paste0(label,"_with_unknown")]] <- seurat[[label]]
    seurat[[paste0(label,"_with_unknown")]][is.na(seurat[[paste0(label,"_with_unknown")]])] <- "unknown"
  }
  if (is.null(seurat.mc) & is.null(membership)) {
    if (length(assay) == 1) {
      if (is.null(graph.name)) {graph.name = "nn"}
      if (is.null(label)) {
        graph <- ComputeUnimodalKnn(seurat = seurat,
                                    k.knn = k.knn, kith = kith, kernel = kernel,
                                    graph.name = graph.name, assay = assay, reduction = reduction,
                                    dims = dims,
                                    verbose = verbose)
      }
      else {
        if (length(which(is.na(seurat[[label]][, 1]))) >
            0) {
          message("using partial annotation")
          unknowns <- colnames(seurat)[which(is.na(seurat[[label]][,
                                                                   1]))]
          graphUnknown <- ComputeUnimodalKnn(seurat = seurat,
                                             k.knn = k.knn, kith = kith, kernel = kernel,
                                             graph.name = graph.name, assay = assay, reduction = reduction,
                                             dims = dims,
                                             verbose = verbose)

          igraph::V(graphUnknown)$label <- seurat[[label]][,
                                                           1]
          unknownAndNeighbors <- unique(c(unknowns, unlist(sapply(which(is.na(seurat[[label]][,
                                                                                              1])), FUN = function(X) {
                                                                                                names(igraph::neighbors(graphUnknown, X))
                                                                                              }))))
          graphUnknown <- igraph::subgraph(graphUnknown,
                                           unknownAndNeighbors)
          edgeList <- igraph::as_edgelist(graphUnknown)
          igraph::E(graphUnknown)$kept <- is.na(igraph::V(graphUnknown)[edgeList[,
                                                                                 1]]$label) | is.na(igraph::V(graphUnknown)[edgeList[,
                                                                                                                                     2]]$label)
          graphUnknown = igraph::subgraph.edges(graphUnknown,
                                                which(igraph::E(graphUnknown)$kept), delete.vertices = FALSE)
        }
        graphList <- lapply(X = as.vector(na.exclude(unique(seurat[[label]][,
                                                                            1]))), FUN = function(X) {
                                                                              ComputeUnimodalKnn(seurat = seurat, k.knn = k.knn,
                                                                                                 kith = kith, kernel = kernel, graph.name = graph.name,
                                                                                                 assay = assay, reduction = reduction, dims = dims,
                                                                                                 label = label, subsetLabel = X,
                                                                                                 verbose = verbose)
                                                                            })
        if (length(which(is.na(seurat[[label]][, 1]))) >
            0) {
          graphList[[length(graphList) + 1]] <- graphUnknown
        }
        if (length(graphList) > 1) {
          graph <- do.call(igraph::union, graphList)
          allWeigths <- lapply(X = 1:length(graphList),
                               FUN = function(X) {
                                 igraph::get.edge.attribute(graph, name = paste0("weight_",
                                                                                 X))
                               })
          dfWeights <- do.call("cbind", allWeigths)
          igraph::E(graph)$weight <- rowSums(dfWeights,
                                             na.rm = T)
          graph = igraph::permute(graph, match(igraph::V(graph)$name,
                                               colnames(seurat)))
        } else {
          graph <- graphList[[1]]
        }
      }
    }
    else {
      if (is.null(graph.name)) {graph.name = "knn"}
      if (is.null(label)) {
        graph <- ComputeMultimodalKnn(seurat = seurat,
                                      k.knn = k.knn, kith = kith, kernel = kernel,
                                      graph.name = graph.name, assay = assay, reduction = reduction,
                                      dims = dims,
                                      verbose = verbose)
      }
      else {
        if (length(which(is.na(seurat[[label]][, 1]))) >
            0) {
          message("using partial annotation")
          unknowns <- colnames(seurat)[which(is.na(seurat[[label]][,
                                                                   1]))]
          graphUnknown <- ComputeMultimodalKnn(seurat = seurat,
                                               k.knn = k.knn, kith = kith, kernel = kernel,
                                               graph.name = graph.name, assay = assay, reduction = reduction,
                                               dims = dims,
                                               verbose = verbose)
          igraph::V(graphUnknown)$label <- seurat[[label]][,
                                                           1]
          unknownAndNeighbors <- unique(c(unknowns, unlist(sapply(which(is.na(seurat[[label]][,
                                                                                              1])), FUN = function(X) {
                                                                                                names(igraph::neighbors(graphUnknown, X))
                                                                                              }))))
          graphUnknown <- igraph::subgraph(graphUnknown,
                                           unknownAndNeighbors)
          edgeList <- igraph::as_edgelist(graphUnknown)
          igraph::E(graphUnknown)$kept <- is.na(igraph::V(graphUnknown)[edgeList[,
                                                                                 1]]$label) | is.na(igraph::V(graphUnknown)[edgeList[,
                                                                                                                                     2]]$label)
          graphUnknown = igraph::subgraph.edges(graphUnknown,
                                                which(igraph::E(graphUnknown)$kept), delete.vertices = FALSE)
        }
        graphList <- lapply(X = na.exclude(unique(seurat[[label]][,
                                                                  1])), FUN = function(X) {
                                                                    ComputeMultimodalKnn(seurat = seurat,
                                                                                         k.knn = k.knn,
                                                                                         kith = kith,
                                                                                         kernel = kernel,
                                                                                         graph.name = graph.name,
                                                                                         assay = assay,
                                                                                         reduction = reduction,
                                                                                         dims = dims,
                                                                                         label = label,
                                                                                         subsetLabel = X,
                                                                                         verbose = verbose)
                                                                  })
        if (length(which(is.na(seurat[[label]][, 1]))) >
            0) {
          graphList[[length(graphList) + 1]] <- graphUnknown
        }
        if (length(graphList) > 1) {
          graph <- do.call(igraph::union, graphList)
          allWeigths <- lapply(X = 1:length(graphList),
                               FUN = function(X) {
                                 igraph::edge_attr(graph, name = paste0("weight_",
                                                                        X))
                               })
          dfWeights <- do.call("cbind", allWeigths)
          if (any(dfWeights == -1, na.rm = T)) {
            medWeights <- median(na.exclude(dfWeights[dfWeights >
                                                        0]))
            dfWeights[which(dfWeights == -1)] <- medWeights
          }
          igraph::E(graph)$weight <- rowSums(dfWeights,
                                             na.rm = T)
          graph = igraph::permute(graph, match(igraph::V(graph)$name,
                                               colnames(seurat)))
        }
        else {
          graph <- graphList[[1]]
        }
      }
    }
    walktrap <- igraph::cluster_walktrap(graph)
    seurat[[paste0("walktrap_clusters_", assay[[1]])]] <- walktrap$membership
    n_target <- .SCResolveTargetMetacells(n_cells = ncol(seurat), gamma = gamma)
    membership <- igraph::cut_at(walktrap, no = n_target)
    names(membership) <- colnames(seurat)
    message("metacells identified")
  }
  else {
    if (is.null(membership) & !is.null(seurat.mc)) {
      walktrap <- seurat.mc@misc$metacells_hierarchy
      n_target <- .SCResolveTargetMetacells(n_cells = ncol(seurat), gamma = gamma)
      membership <- igraph::cut_at(walktrap, no = n_target)
      names(membership) <- colnames(seurat)
    }
    else {
      if (!is.null(membership)) {
        membership <- .SCNormalizeMembership(membership)
        walktrap <- list(membership = membership)
        gamma = floor(length(membership)/length(unique(membership)))
      }
    }
  }
  membership_names <- .SCFormatMetacellNames(membership, prefixMC = prefixMC)
  membership_names <- stats::setNames(as.character(unname(membership_names)), as.character(names(membership)))
  if (anyDuplicated(names(membership_names))) stop("Duplicated single-cell IDs in membership.", call. = FALSE)
  seurat[[paste0("metacell_g", gamma)]] <- membership_names
  metacell_ids <- sort(unique(unname(membership_names)))
  metacell_name_map <- stats::setNames(metacell_ids, metacell_ids)
  membership_table <- data.frame(cell_id = names(membership_names), metacell_id = unname(membership_names), stringsAsFactors = FALSE)
  fragment_membership <- membership_names
  if (return.seurat) {
    assaysToAgg <- Assays(seurat)[sapply(X = Assays(seurat),
                                         FUN = function(X) {
                                           .sc_assay_has_slot(seurat, X, "counts")
                                         })]

    isChromAssay <- sapply(X = assaysToAgg, FUN = function(X) {
      is(GetAssay(seurat,assay = X))[1] == "ChromatinAssay"
    })

    chrom.assay.list <- list()

    for (chromAssay in assaysToAgg[isChromAssay]) {
      if (avg.in.data) {
        chrom.assay.list[[chromAssay]] <- Signac::CreateChromatinAssay(counts =  MetacellExpression(seurat,
                                                                                            assays = chromAssay,
                                                                                            group.by = paste0("metacell_g", gamma),
                                                                                            metacell.names = metacell_name_map,
                                                                                            return.seurat = F)[[chromAssay]],
                                                               genome = genome(seurat[[chromAssay]]),
                                                               ranges = Signac::StringToGRanges(rownames(seurat[[chromAssay]]),
                                                                                                sep = peakSep),
                                                               annotation = Signac::Annotation(seurat[[chromAssay]]))
        chrom.assay.list[[chromAssay]] <- .sc_set_assay_data(
          object = chrom.assay.list[[chromAssay]],
          new.data = MetacellExpression(seurat,
                                        assays = chromAssay,
                                        pb.method = "average",
                                        group.by = paste0("metacell_g", gamma),
                                        metacell.names = metacell_name_map,
                                        layer = "data",
                                        return.seurat = F)[[chromAssay]],
          slot = "data"
        )
      } else {
        chrom.assay.list[[chromAssay]] <- Signac::CreateChromatinAssay(counts =  MetacellExpression(seurat,
                                                                                            assays = chromAssay,
                                                                                            group.by = paste0("metacell_g", gamma),
                                                                                            metacell.names = metacell_name_map,
                                                                                            return.seurat = F)[[chromAssay]],
                                                               genome = genome(seurat[[chromAssay]]),
                                                               ranges = Signac::StringToGRanges(rownames(seurat[[chromAssay]]),
                                                                                                sep = peakSep),
                                                               annotation = Signac::Annotation(seurat[[chromAssay]]))
      }


      if (!is.null(fragmentFiles[[chromAssay]])) {
        if (is.null(tmpPath)) {
          tmpPath <- file.path(tempdir(), paste0("SuperCell_fragments_", Sys.getpid()))
        }
        frag_paths <- as.character(unlist(fragmentFiles[[chromAssay]], use.names = FALSE))
        fragment_results <- lapply(seq_along(frag_paths), function(i) {
          AggregateFragmentFile(input_file = frag_paths[[i]],
                                tmp_path = file.path(tmpPath, paste0("fragment_", chromAssay, "_", sprintf("%03d", i))),
                                output_name = paste0("MC_", i, "_", fs::path_file(frag_paths[[i]])),
                                output_path = outputDirMcFragment,
                                membership = fragment_membership,
                                returnOutputFileName = TRUE,
                                return_details = TRUE,
                                bgzip_path = bgzip_path,
                                tabix_path = tabix_path,
                                nb_cl = nb_cl)
        })
        fragment_cell_maps <- lapply(fragment_results, function(result) {
          ids <- intersect(as.character(result$metacell_ids), colnames(chrom.assay.list[[chromAssay]]))
          stats::setNames(ids, ids)
        })
        empty_maps <- which(!lengths(fragment_cell_maps))
        if (length(empty_maps)) {
          stop("No object metacells from aggregated fragment file(s) were present in the ChromatinAssay: ",
               paste(utils::head(empty_maps, 10L), collapse = ", "), call. = FALSE)
        }
        all_registered_cells <- unlist(lapply(fragment_cell_maps, names), use.names = FALSE)
        if (anyDuplicated(all_registered_cells)) {
          duplicated_cells <- unique(all_registered_cells[duplicated(all_registered_cells)])
          stop("Metacells are assigned to more than one fragment file: ",
               paste(utils::head(duplicated_cells, 10L), collapse = ", "), call. = FALSE)
        }
        if (isTRUE(return_fragment_manifest)) {
          for (i in seq_along(frag_paths)) {
            fragment_manifest_rows[[length(fragment_manifest_rows) + 1L]] <- data.frame(
              assay = chromAssay,
              input_file = normalizePath(frag_paths[[i]], mustWork = FALSE),
              fragment_file = fragment_results[[i]]$fragment_file,
              index_file = fragment_results[[i]]$index_file,
              n_input_membership_cells = length(fragment_membership),
              n_metacells = fragment_results[[i]]$n_metacells,
              n_fragment_rows = fragment_results[[i]]$n_fragment_rows,
              status = "ok",
              stringsAsFactors = FALSE
            )
            fragment_cell_map_rows[[length(fragment_cell_map_rows) + 1L]] <- data.frame(
              assay = chromAssay,
              fragment_file = fragment_results[[i]]$fragment_file,
              object_cell = names(fragment_cell_maps[[i]]),
              fragment_barcode = unname(fragment_cell_maps[[i]]),
              stringsAsFactors = FALSE
            )
          }
        }
        message("Fragment file aggregated")
        mcFragments <- Map(function(result, cell_map) {
          Signac::CreateFragmentObject(path = result$fragment_file,
                                       cells = cell_map,
                                       validate.fragments = TRUE)
        }, result = fragment_results, cell_map = fragment_cell_maps)
        Signac::Fragments(chrom.assay.list[[chromAssay]]) <- mcFragments
      }

    }

    if (length(which(!isChromAssay)) >0 ) {
      if (avg.in.data) {
        std.assay.list <- list()

        for (assay_name in assaysToAgg[!isChromAssay]) {
          std.assay.list[[assay_name]] <- .sc_create_assay_object(counts = MetacellExpression(seurat,
                                                                                    assays = assay_name,
                                                                                    group.by = paste0("metacell_g", gamma),
                                                                                    metacell.names = metacell_name_map,
                                                                                    return.seurat = F)[[assay_name]],
                                                        data =  MetacellExpression(seurat,
                                                                                   assays = assay_name, pb.method = "average",
                                                                                   group.by = paste0("metacell_g", gamma),
                                                                                   metacell.names = metacell_name_map,
                                                                                   layer = "data",
                                                                                   return.seurat = F)[[assay_name]]
          )
        }

        seurat.mc <- .sc_create_seurat_object(std.assay.list[[1]], assay = names(std.assay.list)[1])
        for (std.a in names(std.assay.list)[-1]) {
          seurat.mc[[std.a]] <- std.assay.list[[std.a]]
        }
        for (a in names(chrom.assay.list)) {
          seurat.mc[[a]] <- chrom.assay.list[[a]]
        }
      } else {
        seurat.mc <- MetacellExpression(seurat,
                                        assays = assaysToAgg[!isChromAssay],
                                        group.by = paste0("metacell_g", gamma),
                                        metacell.names = metacell_name_map,
                                        return.seurat = T)
      }
      for (a in names(chrom.assay.list)) {
        seurat.mc[[a]] <- chrom.assay.list[[a]]
      }


    } else {

      seurat.mc <- .sc_create_seurat_object(chrom.assay.list[[1]], assay = names(chrom.assay.list)[1])
      for (a in names(chrom.assay.list)[-1]) {
        seurat.mc[[a]] <- chrom.assay.list[[a]]
      }

    }
    if (metacellNormalization) {
      normalizations <- names(seurat@commands)[startsWith(names(seurat@commands),
                                                          prefix = "NormalizeData")]
      for (n in normalizations) {
        DefaultAssay(seurat.mc) <- seurat@commands[[n]]$assay
        seurat.mc <- NormalizeData(object = seurat.mc,
                                   normalization.method = seurat@commands[[n]]$normalization.method,
                                   scale.factor = seurat@commands[[n]]$scale.factor,
                                   margin = seurat@commands[[n]]$margin, verbose = TRUE)
      }
    }
  }
  else {
    seurat.mc <- list(membership = membership_names, supercell_size = as.numeric(table(factor(membership_names, levels = metacell_ids))),
                      h_membership = walktrap)
  }
  fields <- sapply(X = colnames(seurat@meta.data), FUN = function(X) {
    is.character(seurat[[X]][, 1]) | is.factor(seurat[[X]][,
                                                           1])
  })
  message("metadata assignement")
  for (f in colnames(seurat@meta.data)[fields]) {
    assign_res <- supercell_assign(clusters = seurat[[f]][,1],
                                   supercell_membership = membership_names, method = "absolute")
    purity_res <- supercell_purity(clusters = seurat[[f]][, 1], supercell_membership = membership_names)
    if (!all(colnames(seurat.mc) %in% names(assign_res)) || !all(colnames(seurat.mc) %in% names(purity_res))) {
      stop("Metacell metadata assignment failed: names do not match metacell object colnames.", call. = FALSE)
    }
    seurat.mc[[f]] <- assign_res[colnames(seurat.mc)]
    seurat.mc[[paste0(f, "_purity")]] <- purity_res[colnames(seurat.mc)]
  }
  if (!is.null(label)) {
    # annotate metacell containing only unknown cell as unknown
    seurat.mc[[label]][seurat.mc[[paste0(label,"_purity")]]==0] <- "unknown"
  }
  if (return.seurat) {
    mc_ids <- as.character(colnames(seurat.mc))
    if (!setequal(unique(membership_table$metacell_id), mc_ids)) {
      stop("Membership metacell IDs do not match metacell object colnames.", call. = FALSE)
    }
    membership_table$metacell_id <- factor(membership_table$metacell_id, levels = mc_ids)
    membership_table <- membership_table[order(membership_table$metacell_id), , drop = FALSE]
    membership_table$metacell_id <- as.character(membership_table$metacell_id)
    membership_names <- stats::setNames(membership_table$metacell_id, membership_table$cell_id)
    fragment_manifest <- if (length(fragment_manifest_rows) && isTRUE(return_fragment_manifest)) {
      do.call(rbind, fragment_manifest_rows)
    } else {
      data.frame(assay = character(), input_file = character(), fragment_file = character(), index_file = character(), n_input_membership_cells = integer(), n_metacells = integer(), n_fragment_rows = integer(), status = character(), stringsAsFactors = FALSE)
    }
    fragment_cell_map <- if (length(fragment_cell_map_rows) && isTRUE(return_fragment_manifest)) {
      do.call(rbind, fragment_cell_map_rows)
    } else {
      data.frame(assay = character(), fragment_file = character(), object_cell = character(), fragment_barcode = character(), stringsAsFactors = FALSE)
    }
    seurat.mc$metacell_id <- mc_ids
    seurat.mc$size <- as.integer(table(factor(membership_table$metacell_id, levels = mc_ids)))
    seurat.mc@misc$schema_version <- "supercell2_metacell_v2"
    seurat.mc@misc$metacells_hierarchy <- walktrap
    seurat.mc@misc$walktrap_clusters <- walktrap$membership
    seurat.mc@misc$gamma <- gamma
    seurat.mc@misc$seed <- seed
    seurat.mc@misc$membership <- membership_names
    if (isTRUE(return_membership_table)) {
      seurat.mc@misc$membership_table <- membership_table
      stopifnot(identical(sort(unique(seurat.mc@misc$membership_table$metacell_id)), sort(colnames(seurat.mc))))
      stopifnot(identical(names(seurat.mc@misc$membership), seurat.mc@misc$membership_table$cell_id))
    } else {
      seurat.mc@misc$membership_table <- NULL
    }
    seurat.mc@misc$fragment_manifest <- fragment_manifest
    seurat.mc@misc$fragment_cell_map <- fragment_cell_map
  }
  return(seurat.mc)
}
