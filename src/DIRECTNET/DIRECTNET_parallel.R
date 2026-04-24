Run_DIRECT_NET_parallel <- function(
  object,
  peakcalling = FALSE,
  macs2.path = NULL,
  fragments = NULL,
  k_neigh = 50,
  atacbinary = TRUE,
  max_overlap = 0.8,
  reduction.name = NULL,
  size_factor_normalize = FALSE,
  genome.info,
  focus_markers,
  params = NULL,
  nthread = 1,
  early_stop = FALSE,
  HC_cutoff = NULL,
  LC_cutoff = NULL,
  rescued = FALSE,
  seed = 123,
  verbose = TRUE,
  n_workers = 1,
  parallel_backend = c("mclapply", "lapply"),
  progress_every = 1,
  coordinate_format = c("auto", "directnet_original", "cicero_underscore"),
  debug_coordinate_overlap = FALSE
) {
  parallel_backend <- match.arg(parallel_backend)
  coordinate_format <- match.arg(coordinate_format)

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required but is not available.")
  }
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required but is not available.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but is not available.")
  }
  if (!requireNamespace("cicero", quietly = TRUE)) {
    stop("Package 'cicero' is required but is not available.")
  }

  n_workers <- as.integer(n_workers)
  if (is.na(n_workers) || n_workers < 1) {
    n_workers <- 1L
  }

  nthread <- as.integer(nthread)
  if (is.na(nthread) || nthread < 1) {
    nthread <- 1L
  }

  progress_every <- as.integer(progress_every)
  if (is.na(progress_every) || progress_every < 1) {
    progress_every <- 100L
  }

  if (peakcalling) {
    if (verbose) {
      message("Calling Peak")
    }

    object$cluster <- Idents(object)

    if (is.null(macs2.path)) {
      stop("Please give the path to macs2!")
    }
    if (is.null(fragments)) {
      stop("Please input fragments!")
    }

    peaks_called <- CallPeaks(
      object = object,
      group.by = "cluster",
      macs2.path = macs2.path
    )

    new_atac_data <- FeatureMatrix(
      fragments = fragments,
      features = peaks_called
    )

    object@assays$ATAC@counts <- new_atac_data

    if (verbose) {
      message("Peak calling finished")
    }
  }

  if (verbose) {
    message("Generating aggregated data")
  }

  if ("aggregated.data" %in% names(Misc(object))) {
    agg.data <- Misc(object, slot = "aggregated.data")
  } else {
    agg.data <- Aggregate_data(
      object,
      k_neigh = k_neigh,
      atacbinary = atacbinary,
      max_overlap = max_overlap,
      reduction.name = reduction.name,
      size_factor_normalize = size_factor_normalize,
      seed = seed,
      verbose = verbose
    )
    Misc(object, slot = "aggregated.data") <- agg.data
  }

  options(stringsAsFactors = FALSE)

  if (is.null(params)) {
    params <- list(
      eta = 0.3,
      max_depth = 6,
      min_child_weight = 1,
      subsample = 1,
      colsample_bytree = 1,
      lambda = 1
    )
  }

  has_rna <- "rna" %in% names(agg.data)

  if (has_rna) {
    data_rna <- as.matrix(agg.data$rna)
    rna <- rownames(data_rna)
    rna <- unlist(lapply(rna, function(x) strsplit(x, "[.]")[[1]][1]))
    rownames(data_rna) <- rna
    data_rna <- data_rna[!duplicated(rna), , drop = FALSE]
  } else {
    data_rna <- NULL
  }

  data_atac_raw <- as.matrix(agg.data$atac)
  raw_peak_names <- rownames(data_atac_raw)

  genome.info <- as.data.frame(genome.info, stringsAsFactors = FALSE)
  genome.info$genes <- unlist(lapply(genome.info$genes, function(x) strsplit(x, "[|]")[[1]][1]))
  genome.info$genes <- unlist(lapply(genome.info$genes, function(x) strsplit(x, "[.]")[[1]][1]))
  genome.info <- genome.info[!duplicated(genome.info$genes), , drop = FALSE]
  genome.info$Starts <- as.numeric(genome.info$Starts)
  genome.info$Ends <- as.numeric(genome.info$Ends)
  genome.info <- genome.info[!is.na(genome.info$Chrom) & !is.na(genome.info$Starts) & !is.na(genome.info$Ends), , drop = FALSE]

  focus_markers <- unique(unlist(lapply(focus_markers, function(x) strsplit(x, "[.]")[[1]][1])))
  focus_markers <- genome.info$genes[genome.info$genes %in% focus_markers]
  genome.info.used <- genome.info[match(focus_markers, genome.info$genes), , drop = FALSE]

  Chr <- genome.info.used$Chrom
  Starts <- genome.info.used$Starts
  Ends <- genome.info.used$Ends

  if (length(focus_markers) == 0) {
    stop("No focus_markers overlap genome.info$genes after normalization.")
  }

  build_peak_names <- function(format) {
    if (format == "cicero_underscore") {
      gsub("[:-]", "_", raw_peak_names)
    } else {
      # convert chr1_100_200 → chr1:100_200
      sub("^([^_]+)_", "\\1:", raw_peak_names)
    }
  }

  make_region <- function(chr, start, end, format) {
    if (format == "cicero_underscore") {
      paste(chr, start, end, sep = "_")
    } else {
      paste(chr, ":", start, "-", end, sep = "")
    }
  }

  score_coordinate_format <- function(format, n_test = 100L) {
    test_indices <- seq_len(min(length(focus_markers), n_test))
    test_peaks <- build_peak_names(format)
    sum(vapply(test_indices, function(j) {
      chr_j <- genome.info.used$Chrom[j]
      start_j <- as.numeric(genome.info.used$Starts[j])
      p2 <- make_region(chr_j, start_j - 250000, start_j + 250000, format)
      length(cicero::find_overlapping_coordinates(test_peaks, p2)) > 1
    }, logical(1)))
  }

  if (coordinate_format == "auto") {
    directnet_score <- score_coordinate_format("directnet_original")
    cicero_score <- score_coordinate_format("cicero_underscore")
    if (cicero_score > directnet_score) {
      coordinate_format <- "cicero_underscore"
    } else if (directnet_score > cicero_score) {
      coordinate_format <- "directnet_original"
    } else if (any(grepl("^chr[^_]+_[0-9]+_[0-9]+$", raw_peak_names))) {
      coordinate_format <- "cicero_underscore"
    } else {
      coordinate_format <- "directnet_original"
    }
    if (verbose) {
      message("Selected coordinate_format=", coordinate_format, " using directnet_original score=", directnet_score, ", cicero_underscore score=", cicero_score)
    }
  }

  data_atac <- data_atac_raw
  rownames(data_atac) <- build_peak_names(coordinate_format)
  peaks <- rownames(data_atac)

  if (isTRUE(debug_coordinate_overlap)) {
    valid_debug_idx <- which(!is.na(genome.info.used$Starts) & !is.na(genome.info.used$Chrom))
    if (length(valid_debug_idx) > 0) {
      test_i <- valid_debug_idx[1]
      chr_test <- genome.info.used$Chrom[test_i]
      start_test <- as.numeric(genome.info.used$Starts[test_i])
      p1_test <- make_region(chr_test, start_test - 500, start_test, coordinate_format)
      p2_test <- make_region(chr_test, start_test - 250000, start_test + 250000, coordinate_format)
      message("Coordinate format: ", coordinate_format)
      message("Example ATAC peaks:")
      message(paste(utils::head(peaks, 10), collapse = "
"))
      message("Example promoter query: ", p1_test)
      message("Example enhancer query: ", p2_test)
      message("Promoter overlap count: ", length(cicero::find_overlapping_coordinates(peaks, p1_test)))
      message("Enhancer overlap count: ", length(cicero::find_overlapping_coordinates(peaks, p2_test)))
    }
  }

  if (verbose) {
    message(
      "Inferring links for ", length(focus_markers),
      " genes using n_workers=", n_workers,
      ", xgboost nthread=", nthread,
      ", backend=", parallel_backend,
      ", coordinate_format=", coordinate_format
    )
  }

  infer_one_gene <- function(i) {
    if (verbose && (i == 1L || i %% progress_every == 0L)) {
      message("Inferring links for gene ", i, "/", length(focus_markers), ": ", focus_markers[i])
    }

    set.seed(seed + i)

    p1 <- make_region(Chr[i], Starts[i] - 500, Starts[i], coordinate_format)
    p2 <- make_region(Chr[i], Starts[i] - 250000, Starts[i] + 250000, coordinate_format)

    promoters <- cicero::find_overlapping_coordinates(peaks, p1)
    enhancers <- cicero::find_overlapping_coordinates(peaks, p2)
    enhancers <- setdiff(enhancers, promoters)

    if (has_rna) {
      idx <- which(rownames(data_rna) == focus_markers[i])
    } else {
      idx <- 1L
    }

    if (!((length(promoters) > 0 && length(enhancers) > 1) && length(idx) != 0)) {
      # if (verbose) {
      #   message("There are less than two peaks detected within 500 kb for ", focus_markers[i])
      # }
      return(list(i = i, conns = NULL, X = NULL, Y = NULL, error = NULL, reason = "peak_filter"))
    }

    id1 <- match(promoters, peaks)
    id1 <- id1[!is.na(id1)]
    id2 <- match(enhancers, peaks)
    id2 <- id2[!is.na(id2)]
    id2_new <- setdiff(id2, id1)

    if (length(id1) == 0 || length(id2_new) <= 1) {
      # if (verbose) {
      #   message("There are less than two usable peaks detected within 500 kb for ", focus_markers[i])
      # }
      return(list(i = i, conns = NULL, X = NULL, Y = NULL, error = NULL, reason = "usable_peak_filter"))
    }

    X <- data_atac[id2_new, , drop = FALSE]
    if (nrow(X) < 2 || ncol(X) < 2) {
      if (verbose) {
        message("Skipping ", focus_markers[i], ": insufficient X dimensions after filtering; dim(X)=", paste(dim(X), collapse = " x "))
      }
      return(list(i = i, conns = NULL, X = NULL, Y = NULL, error = NULL, reason = "x_dim_filter"))
    }

    Y_mat <- data_atac[id1, , drop = FALSE]
    Y_vec <- if (nrow(Y_mat) > 1) Matrix::colSums(Y_mat) else as.numeric(Y_mat[1, ])
    Y <- matrix(Y_vec, nrow = 1)
    rownames(Y) <- peaks[id1[1]]
    colnames(Y) <- colnames(data_atac)

    if (has_rna) {
      Z_mat <- data_rna[idx, , drop = FALSE]
      Z_vec <- if (nrow(Z_mat) > 1) Matrix::colSums(Z_mat) else as.numeric(Z_mat[1, ])
      Z <- matrix(Z_vec, nrow = 1)
      rownames(Z) <- focus_markers[i]
      colnames(Z) <- colnames(data_rna)
    } else {
      Z <- Y
    }

    X <- as.matrix(X)
    Y <- as.matrix(Y)
    Z <- as.matrix(Z)

    z_vec <- as.numeric(Z[1, ])
    if (length(z_vec) != ncol(X)) {
      if (verbose) {
        message("Skipping ", focus_markers[i], ": response length does not match X columns; length(Z)=", length(z_vec), ", ncol(X)=", ncol(X))
      }
      return(list(i = i, conns = NULL, X = NULL, Y = NULL, error = NULL, reason = "response_length_mismatch"))
    }

    if (any(!is.finite(z_vec))) {
      if (verbose) {
        message("Skipping ", focus_markers[i], ": non-finite response values in Z")
      }
      return(list(i = i, conns = NULL, X = NULL, Y = NULL, error = NULL, reason = "nonfinite_response"))
    }


    if (early_stop) {
      n_obs <- ncol(X)
      cv_idx <- sample(1:5, size = n_obs, replace = TRUE)
      test_idx <- which(cv_idx == 1)
      validation_idx <- which(cv_idx == 2)
      holdout_idx <- c(test_idx, validation_idx)

      if (length(test_idx) == 0 || length(validation_idx) == 0 || length(holdout_idx) >= n_obs) {
        nrounds_final <- 100
      } else {
        x_train <- as.matrix(X[, -holdout_idx, drop = FALSE])
        x_validation <- as.matrix(X[, validation_idx, drop = FALSE])
        y_vec <- as.numeric(Y[1, ])
        y_train <- y_vec[-holdout_idx]
        y_validation <- y_vec[validation_idx]

        dtrain <- xgboost::xgb.DMatrix(data = t(x_train), label = y_train)
        dvalidation <- xgboost::xgb.DMatrix(data = t(x_validation), label = y_validation)

        xgb_v <- xgboost::xgb.train(
          params = params,
          data = dtrain,
          watchlist = list(train = dtrain, test = dvalidation),
          nrounds = 100,
          nthread = nthread,
          objective = "reg:squarederror",
          verbose = 0
        )

        cv1 <- xgb_v$evaluation_log
        rmse_d <- cv1$test_rmse - cv1$train_rmse
        if (length(rmse_d) >= 2 && is.finite(rmse_d[1]) && rmse_d[1] != 0) {
          rmse_dd <- abs(diff(rmse_d)) / abs(rmse_d[1])
          if (all(is.na(rmse_dd))) {
            nrounds_final <- 100
          } else {
            stop_index <- which(rmse_dd == min(rmse_dd, na.rm = TRUE))
            nrounds_final <- if (length(stop_index) > 0 && is.finite(stop_index[1])) stop_index[1] else 100
          }
        } else {
          nrounds_final <- 100
        }
      }
    } else {
      nrounds_final <- 100
    }

    xgb.fit.final <- xgboost::xgboost(
      params = params,
      data = t(X),
      label = z_vec,
      nrounds = nrounds_final,
      nthread = nthread,
      objective = "reg:squarederror",
      verbose = 0
    )

    importance_matrix <- tryCatch(
      xgboost::xgb.importance(model = xgb.fit.final),
      error = function(e) NULL
    )

    if (is.null(importance_matrix) || nrow(importance_matrix) == 0) {
      return(list(i = i, conns = NULL, X = if (rescued) X else NULL, Y = if (rescued) Y else NULL, error = NULL, reason = "empty_importance"))
    }

    conns_h <- data.frame(
      Peak1 = as.character(rownames(Y)),
      Peak2 = as.character(importance_matrix$Feature),
      Importance = importance_matrix$Gain,
      stringsAsFactors = FALSE
    )

    list(i = i, conns = conns_h, X = if (rescued) X else NULL, Y = if (rescued) Y else NULL, error = NULL, reason = "success")
  }

  safe_infer_one_gene <- function(i) {
    tryCatch(
      infer_one_gene(i),
      error = function(e) {
        msg <- paste0("DIRECT-NET failed for gene ", i, " / ", length(focus_markers), " (", focus_markers[i], "): ", conditionMessage(e))
        if (verbose) {
          message(msg)
        }
        list(i = i, conns = NULL, X = NULL, Y = NULL, error = msg, reason = "error")
      }
    )
  }

  gene_indices <- seq_along(focus_markers)

  if (n_workers == 1L || parallel_backend == "lapply") {
    gene_results <- lapply(gene_indices, safe_infer_one_gene)
  } else {
    gene_results <- parallel::mclapply(
      gene_indices,
      safe_infer_one_gene,
      mc.cores = n_workers,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE
    )
  }

  DIRECT_NET_Result <- vector("list", length(focus_markers))
  TXs <- vector("list", length(focus_markers))
  TYs <- vector("list", length(focus_markers))
  failed_results <- 0L

  for (res in gene_results) {
    if (!is.list(res) || is.null(res$i)) {
      failed_results <- failed_results + 1L
      next
    }
    i <- res$i
    DIRECT_NET_Result[[i]] <- res$conns
    if (rescued) {
      TXs[[i]] <- res$X
      TYs[[i]] <- res$Y
    }
  }

  failed_gene_messages <- unlist(lapply(gene_results, function(x) {
    if (is.list(x) && !is.null(x$error)) x$error else NULL
  }))

  if (verbose && length(failed_gene_messages) > 0) {
    message("DIRECT-NET skipped ", length(failed_gene_messages), " genes due to per-gene errors.")
    message("First few per-gene errors:")
    message(paste(utils::head(failed_gene_messages, 10), collapse = "\n"))
  }

  if (verbose && failed_results > 0) {
    message("Skipped ", failed_results, " malformed parallel results.")
  }

  if (verbose) {
    skip_reasons <- unlist(lapply(gene_results, function(x) {
      if (is.list(x) && !is.null(x$reason)) x$reason else "malformed_result"
    }))
    message("DIRECT-NET per-gene outcome summary:")
    message(paste(capture.output(print(sort(table(skip_reasons), decreasing = TRUE))), collapse = "
"))
  }

  nonnull <- vapply(DIRECT_NET_Result, function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, logical(1))

  if (!any(nonnull)) {
    warning("DIRECT-NET produced no links. Saving empty direct.net result.")
    DIRECT_NET_Result_all <- data.frame(
      gene = character(),
      Chr = character(),
      Starts = numeric(),
      Ends = numeric(),
      Peak1 = character(),
      Peak2 = character(),
      Importance = numeric(),
      function_type = character(),
      stringsAsFactors = FALSE
    )
    Misc(object, slot = "direct.net") <- DIRECT_NET_Result_all
    return(object)
  }

  conns <- do.call(rbind, DIRECT_NET_Result[nonnull])

  if (is.null(HC_cutoff)) {
    HC_cutoff <- max(stats::quantile(conns$Importance, 0.50, na.rm = TRUE), 0.001)
  }
  if (is.null(LC_cutoff)) {
    LC_cutoff <- min(0.001, stats::quantile(conns$Importance, 0.25, na.rm = TRUE))
  }

  if (verbose) {
    message("HC_cutoff=", HC_cutoff, "; LC_cutoff=", LC_cutoff)
  }

  for (i in seq_along(DIRECT_NET_Result)) {
    if (!is.null(DIRECT_NET_Result[[i]]) && is.data.frame(DIRECT_NET_Result[[i]]) && nrow(DIRECT_NET_Result[[i]]) > 0) {
      conns_h <- DIRECT_NET_Result[[i]]
      Imp_value <- conns_h$Importance
      index1 <- which(Imp_value > HC_cutoff)
      index2 <- intersect(which(Imp_value > LC_cutoff), which(Imp_value <= HC_cutoff))
      index3 <- which(Imp_value <= LC_cutoff)
      function_type <- rep(NA_character_, length(Imp_value))
      function_type[index1] <- "HC"
      function_type[index2] <- "MC"
      function_type[index3] <- "LC"

      if (rescued && i <= length(TXs) && !is.null(TXs[[i]]) && !is.null(TYs[[i]])) {
        X <- TXs[[i]]
        if (nrow(X) > 1) {
          CPi <- suppressWarnings(abs(stats::cor(t(X))))
          CPi[!is.finite(CPi)] <- 0
          diag(CPi) <- 0
          hic_index <- which(rownames(X) %in% conns_h$Peak2[index1])
          other_index <- which(rownames(X) %in% conns_h$Peak2[-index1])
          if (length(hic_index) > 0 && length(other_index) > 0) {
            CPi_sub <- CPi[hic_index, other_index, drop = FALSE]
            flag_matrix <- matrix(0, nrow = nrow(CPi_sub), ncol = ncol(CPi_sub))
            flag_matrix[which(CPi_sub > 0.25)] <- 1
            correlated_index <- which(colSums(flag_matrix) > 0)
            if (length(correlated_index) > 0) {
              function_type[conns_h$Peak2 %in% rownames(X)[other_index[correlated_index]]] <- "HC"
            }
          }
        }
      }

      DIRECT_NET_Result[[i]] <- cbind(
        data.frame(
          gene = focus_markers[i],
          Chr = Chr[i],
          Starts = Starts[i],
          Ends = Ends[i],
          stringsAsFactors = FALSE
        ),
        cbind(conns_h, function_type = function_type)
      )
    }
  }

  nonnull_final <- vapply(DIRECT_NET_Result, function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0 && "function_type" %in% colnames(x), logical(1))

  if (!any(nonnull_final)) {
    warning("DIRECT-NET produced no labeled links. Saving empty direct.net result.")
    DIRECT_NET_Result_all <- data.frame(
      gene = character(),
      Chr = character(),
      Starts = numeric(),
      Ends = numeric(),
      Peak1 = character(),
      Peak2 = character(),
      Importance = numeric(),
      function_type = character(),
      stringsAsFactors = FALSE
    )
  } else {
    DIRECT_NET_Result_all <- do.call(rbind, DIRECT_NET_Result[nonnull_final])
    rownames(DIRECT_NET_Result_all) <- NULL
    DIRECT_NET_Result_all$Starts <- as.numeric(DIRECT_NET_Result_all$Starts)
    DIRECT_NET_Result_all$Ends <- as.numeric(DIRECT_NET_Result_all$Ends)
    DIRECT_NET_Result_all$Importance <- as.numeric(DIRECT_NET_Result_all$Importance)
  }

  Misc(object, slot = "direct.net") <- DIRECT_NET_Result_all
  return(object)
}
