#' Default configuration for syncfmri pipeline
#'
#' @return A named list with default pipeline settings.
#' @export
syncfmri_default_config <- function() {
  list(
    timecourse_derivative = "gimmefMRI",
    timecourse_file_regex = "\\.tsv$",
    roi_x = "LR_postJHipp_200",
    roi_y = "LR_mPFC_200",
    run_boundaries = data.frame(
      run_id = c("run-1", "run-2", "run-3"),
      start_tp = c(1L, 360L, 688L),
      end_tp = c(359L, 687L, 1012L)
    ),
    tr_seconds = 2,
    window_size_tp = 12L,
    step_size_tp = 1L,
    min_points_per_window = 8L,
    regress_pair_mean_signal = FALSE,
    output_derivative_name = "syncfmri"
  )
}

#' Compute run-aware sliding-window connectivity for one ROI pair
#'
#' Windows are computed within each run boundary and never cross run edges.
#' Pearson correlations are transformed using Fisher z.
#'
#' @param x Numeric vector for ROI X.
#' @param y Numeric vector for ROI Y.
#' @param run_boundaries Data frame with columns `run_id`, `start_tp`, `end_tp`.
#' @param window_size_tp Integer window length in timepoints.
#' @param step_size_tp Integer step size in timepoints.
#' @param min_points_per_window Minimum non-missing paired timepoints required.
#' @param regress_pair_mean_signal Logical; if TRUE, regress each ROI on pair mean.
#' @param diagnostic_context Optional label included in boundary validation errors.
#'
#' @return A data frame of sliding-window connectivity values.
#' @export
compute_sliding_connectivity <- function(
    x,
    y,
    run_boundaries,
    window_size_tp,
    step_size_tp = 1L,
    min_points_per_window = 8L,
    regress_pair_mean_signal = FALSE,
    diagnostic_context = NULL) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  .validate_run_boundaries(
    run_boundaries,
    n_tp = length(x),
    diagnostic_context = diagnostic_context
  )

  x_work <- as.numeric(x)
  y_work <- as.numeric(y)

  if (isTRUE(regress_pair_mean_signal)) {
    pair_mean <- rowMeans(cbind(x_work, y_work), na.rm = TRUE)
    pair_mean[is.nan(pair_mean)] <- NA_real_
    x_work <- .residualize_on_signal(x_work, pair_mean)
    y_work <- .residualize_on_signal(y_work, pair_mean)
  }

  rows <- vector("list", nrow(run_boundaries))

  for (i in seq_len(nrow(run_boundaries))) {
    run_id <- as.character(run_boundaries$run_id[[i]])
    run_start <- as.integer(run_boundaries$start_tp[[i]])
    run_end <- as.integer(run_boundaries$end_tp[[i]])

    starts <- seq.int(
      from = run_start,
      to = max(run_start, run_end - as.integer(window_size_tp) + 1L),
      by = as.integer(step_size_tp)
    )
    starts <- starts[starts + as.integer(window_size_tp) - 1L <= run_end]

    if (length(starts) == 0L) {
      rows[[i]] <- data.frame(
        run_id = character(0),
        window_index = integer(0),
        window_start_tp = integer(0),
        window_end_tp = integer(0),
        window_center_tp = numeric(0),
        n_valid_tp = integer(0),
        pearson_r = numeric(0),
        fisher_z = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    run_rows <- lapply(seq_along(starts), function(j) {
      start_tp <- starts[[j]]
      end_tp <- start_tp + as.integer(window_size_tp) - 1L
      idx <- start_tp:end_tp

      valid <- stats::complete.cases(x_work[idx], y_work[idx])
      n_valid <- sum(valid)

      r_val <- NA_real_
      z_val <- NA_real_

      if (n_valid >= as.integer(min_points_per_window)) {
        x_valid <- x_work[idx][valid]
        y_valid <- y_work[idx][valid]

        if (stats::sd(x_valid) > 0 && stats::sd(y_valid) > 0) {
          r_val <- stats::cor(x_valid, y_valid, method = "pearson")
          r_val <- min(max(r_val, -0.999999), 0.999999)
          z_val <- atanh(r_val)
        }
      }

      data.frame(
        run_id = run_id,
        window_index = as.integer(j),
        window_start_tp = as.integer(start_tp),
        window_end_tp = as.integer(end_tp),
        window_center_tp = (start_tp + end_tp) / 2,
        n_valid_tp = as.integer(n_valid),
        pearson_r = r_val,
        fisher_z = z_val,
        stringsAsFactors = FALSE
      )
    })

    rows[[i]] <- do.call(rbind, run_rows)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Run syncfmri pipeline on a BIDS root
#'
#' This function discovers ROI timecourse TSV files in a derivative, computes
#' run-aware sliding-window connectivity for a single ROI pair, and writes
#' subject-level and group-level outputs.
#'
#' @param bids_root Path to BIDS root directory.
#' @param config Named list of pipeline settings.
#'
#' @return A list containing subject and group output paths.
#' @export
syncfmri_run_pipeline <- function(bids_root, config = syncfmri_default_config()) {
  .validate_pipeline_config(config)

  bids_root <- normalizePath(bids_root, winslash = "/", mustWork = TRUE)
  input_root <- file.path(bids_root, "derivatives", config$timecourse_derivative)
  output_root <- file.path(bids_root, "derivatives", config$output_derivative_name)

  if (!dir.exists(input_root)) {
    stop(
      sprintf("Input derivative directory does not exist: %s", input_root),
      call. = FALSE
    )
  }

  files <- .discover_timecourse_files(input_root, config$timecourse_file_regex)
  if (length(files) == 0L) {
    stop("No timecourse files matched the configured pattern.", call. = FALSE)
  }

  subject_tables <- lapply(files, function(one_file) {
    .process_single_timecourse(
      file_path = one_file,
      output_root = output_root,
      config = config
    )
  })

  group_table <- do.call(rbind, subject_tables)
  rownames(group_table) <- NULL

  group_tsv <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.tsv")
  group_json <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.json")
  group_rds <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.rds")
  plot_png <- file.path(output_root, "group", "figures", "group_desc-slidingconn_heatmap.png")
  plot_json <- file.path(output_root, "group", "figures", "group_desc-slidingconn_heatmap.json")

  dir.create(dirname(group_tsv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_png), recursive = TRUE, showWarnings = FALSE)

  readr::write_tsv(group_table, group_tsv, na = "n/a")
  saveRDS(group_table, group_rds)

  group_meta <- list(
    Description = "Group-level sliding window ROI connectivity values.",
    Columns = list(
      subject = "BIDS subject label",
      session = "BIDS session label",
      source_file = "Input ROI timecourse TSV",
      run_id = "Run identifier",
      window_index = "Window number within run",
      window_start_tp = "Window start timepoint",
      window_end_tp = "Window end timepoint",
      window_center_tp = "Window midpoint in timepoints",
      time_seconds = "Window midpoint in seconds",
      n_valid_tp = "Count of valid paired samples used in window",
      pearson_r = "Pearson correlation within window",
      fisher_z = "Fisher z transform of Pearson correlation"
    ),
    Parameters = config
  )
  jsonlite::write_json(group_meta, group_json, pretty = TRUE, auto_unbox = TRUE)

  .plot_group_heatmap(group_table, plot_png)
  jsonlite::write_json(
    list(
      Description = "Heatmap of Fisher-z sliding-window connectivity values.",
      ColorMap = "viridis",
      XAxis = "Window midpoint (seconds)",
      YAxis = "Subject",
      Fill = "Fisher z"
    ),
    plot_json,
    pretty = TRUE,
    auto_unbox = TRUE
  )

  list(
    input_root = input_root,
    output_root = output_root,
    n_files_processed = length(files),
    group_tsv = group_tsv,
    group_json = group_json,
    group_rds = group_rds,
    group_plot = plot_png,
    group_plot_json = plot_json
  )
}

.process_single_timecourse <- function(file_path, output_root, config) {
  tbl <- readr::read_tsv(
    file_path,
    show_col_types = FALSE,
    progress = FALSE,
    skip_empty_rows = FALSE,
    na = c("NA", "NaN", "n/a")
  )

  roi_x <- config$roi_x
  roi_y <- config$roi_y

  if (!(roi_x %in% names(tbl))) {
    stop(sprintf("ROI column not found: %s in %s", roi_x, file_path), call. = FALSE)
  }
  if (!(roi_y %in% names(tbl))) {
    stop(sprintf("ROI column not found: %s in %s", roi_y, file_path), call. = FALSE)
  }

  x <- as.numeric(tbl[[roi_x]])
  y <- as.numeric(tbl[[roi_y]])
  x[is.na(x)] <- NaN
  y[is.na(y)] <- NaN
  entities <- .extract_entities_from_path(file_path)

  sw <- compute_sliding_connectivity(
    x = x,
    y = y,
    run_boundaries = config$run_boundaries,
    window_size_tp = config$window_size_tp,
    step_size_tp = config$step_size_tp,
    min_points_per_window = config$min_points_per_window,
    regress_pair_mean_signal = config$regress_pair_mean_signal,
    diagnostic_context = sprintf(
      "subject=%s; session=%s; source_file=%s",
      entities$subject,
      entities$session,
      basename(file_path)
    )
  )

  sw$subject <- entities$subject
  sw$session <- entities$session
  sw$source_file <- basename(file_path)
  sw$time_seconds <- (sw$window_center_tp - 1) * config$tr_seconds

  out_dir <- file.path(output_root, entities$subject, entities$session, "func")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix <- paste(
    c(entities$subject, entities$session, "desc-slidingconn_timeseries"),
    collapse = "_"
  )

  tsv_path <- file.path(out_dir, paste0(out_prefix, ".tsv"))
  json_path <- file.path(out_dir, paste0(out_prefix, ".json"))
  rds_path <- file.path(out_dir, paste0(out_prefix, ".rds"))

  readr::write_tsv(sw, tsv_path, na = "n/a")
  saveRDS(sw, rds_path)

  sidecar <- list(
    Description = "Sliding-window ROI connectivity values for one subject/session.",
    SourceFile = basename(file_path),
    ROIPair = list(roi_x = roi_x, roi_y = roi_y),
    Parameters = config
  )
  jsonlite::write_json(sidecar, json_path, pretty = TRUE, auto_unbox = TRUE)

  sw
}

.discover_timecourse_files <- function(input_root, file_regex) {
  candidate_dirs <- Sys.glob(file.path(input_root, "sub-*", "ses-*", "func"))
  if (length(candidate_dirs) == 0L) {
    candidate_dirs <- Sys.glob(file.path(input_root, "sub-*", "func"))
  }

  files <- unlist(lapply(candidate_dirs, function(d) {
    list.files(d, pattern = file_regex, full.names = TRUE)
  }), use.names = FALSE)

  sort(unique(files))
}

.extract_entities_from_path <- function(file_path) {
  normalized <- gsub("\\\\", "/", file_path)

  subject <- .first_or_default(
    regmatches(normalized, gregexpr("sub-[^/]+", normalized, perl = TRUE))[[1]],
    "sub-unknown"
  )
  session <- .first_or_default(
    regmatches(normalized, gregexpr("ses-[^/]+", normalized, perl = TRUE))[[1]],
    "ses-01"
  )

  list(subject = subject, session = session)
}

.first_or_default <- function(x, default) {
  if (length(x) == 0L || is.na(x[[1]]) || x[[1]] == "") {
    return(default)
  }
  x[[1]]
}

.residualize_on_signal <- function(y, x) {
  keep <- stats::complete.cases(y, x)
  out <- rep(NA_real_, length(y))

  if (sum(keep) < 3L) {
    return(out)
  }

  fit <- stats::lm(y[keep] ~ x[keep])
  out[keep] <- stats::residuals(fit)
  out
}

.validate_run_boundaries <- function(run_boundaries, n_tp, diagnostic_context = NULL) {
  required <- c("run_id", "start_tp", "end_tp")
  if (!all(required %in% names(run_boundaries))) {
    stop("run_boundaries must include run_id, start_tp, end_tp.", call. = FALSE)
  }

  starts <- as.integer(run_boundaries$start_tp)
  ends <- as.integer(run_boundaries$end_tp)

  if (any(is.na(starts)) || any(is.na(ends))) {
    stop("run_boundaries start_tp/end_tp must be integers.", call. = FALSE)
  }
  if (any(starts < 1L) || any(ends < starts) || any(ends > n_tp)) {
    run_tbl <- data.frame(
      run_id = as.character(run_boundaries$run_id),
      start_tp = starts,
      end_tp = ends,
      stringsAsFactors = FALSE
    )
    run_tbl_text <- paste(capture.output(print(run_tbl, row.names = FALSE)), collapse = "\n")

    ctx <- if (is.null(diagnostic_context) || diagnostic_context == "") {
      "context=unknown"
    } else {
      diagnostic_context
    }

    msg <- paste(
      "run_boundaries are outside the available timeseries length.",
      paste0("Context: ", ctx),
      paste0("Available timeseries length (n_tp): ", n_tp),
      paste0("Configured min start_tp: ", min(starts), "; max end_tp: ", max(ends)),
      "Configured run_boundaries:",
      run_tbl_text,
      sep = "\n"
    )

    stop(msg, call. = FALSE)
  }
}

.validate_pipeline_config <- function(config) {
  required <- c(
    "timecourse_derivative",
    "timecourse_file_regex",
    "roi_x",
    "roi_y",
    "run_boundaries",
    "tr_seconds",
    "window_size_tp",
    "step_size_tp",
    "min_points_per_window",
    "regress_pair_mean_signal",
    "output_derivative_name"
  )

  missing <- setdiff(required, names(config))
  if (length(missing) > 0L) {
    stop(
      sprintf("Missing required config fields: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.data.frame(config$run_boundaries)) {
    stop("config$run_boundaries must be a data.frame.", call. = FALSE)
  }

  if (as.integer(config$window_size_tp) < 2L) {
    stop("config$window_size_tp must be >= 2.", call. = FALSE)
  }

  if (as.integer(config$step_size_tp) < 1L) {
    stop("config$step_size_tp must be >= 1.", call. = FALSE)
  }

  if (as.integer(config$min_points_per_window) < 3L) {
    stop("config$min_points_per_window must be >= 3.", call. = FALSE)
  }
}

.plot_group_heatmap <- function(group_table, png_path) {
  plot_tbl <- group_table
  plot_tbl$subject_session <- paste(plot_tbl$subject, plot_tbl$session, sep = "_")
  plot_tbl <- plot_tbl[stats::complete.cases(plot_tbl$fisher_z), , drop = FALSE]

  if (nrow(plot_tbl) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No valid windows to plot") +
      ggplot2::theme_void()

    ggplot2::ggsave(
      filename = png_path,
      plot = p,
      width = 12,
      height = 8,
      dpi = 150
    )
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(
      x = rlang::.data$time_seconds,
      y = rlang::.data$subject_session,
      fill = rlang::.data$fisher_z
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(ggplot2::vars(rlang::.data$run_id), scales = "free_x", ncol = 1) +
    viridis::scale_fill_viridis(option = "viridis", na.value = "grey80") +
    ggplot2::labs(
      title = "Sliding-window Fisher-z Connectivity",
      x = "Time (seconds)",
      y = "Subject",
      fill = "Fisher z"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  ggplot2::ggsave(
    filename = png_path,
    plot = p,
    width = 12,
    height = 8,
    dpi = 150
  )
}
