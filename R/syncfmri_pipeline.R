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
    event_seconds = numeric(0),
    peak_prominence_quantile = 0.75,
    peak_min_distance_seconds = NA_real_,
    peak_min_height = NA_real_,
    peak_positive_only = TRUE,
    peak_event_max_distance_seconds = Inf,
    subgroup_column = "group",
    subgroup_levels = c("YA", "OA"),
    peak_density_bin_seconds = 4,
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
#' @param nuisance_signal Optional numeric nuisance vector used when
#'   `regress_pair_mean_signal = TRUE`.
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
    nuisance_signal = NULL,
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
    if (is.null(nuisance_signal)) {
      stop(
        "regress_pair_mean_signal = TRUE requires a nuisance_signal vector.",
        call. = FALSE
      )
    }

    nuisance_signal <- as.numeric(nuisance_signal)
    if (length(nuisance_signal) != length(x_work)) {
      stop("nuisance_signal must be the same length as x and y.", call. = FALSE)
    }

    x_work <- .residualize_on_signal(x_work, nuisance_signal)
    y_work <- .residualize_on_signal(y_work, nuisance_signal)
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
#' @param plot_only Logical; if TRUE, skip recomputation and render outputs
#'   from previously saved group windows.
#' @param overwrite Logical; if TRUE, regenerate outputs even when cached
#'   section outputs are available on disk.
#'
#' @return A list containing subject and group output paths.
#' @export
syncfmri_run_pipeline <- function(
    bids_root,
    config = syncfmri_default_config(),
  plot_only = FALSE,
  overwrite = FALSE) {
  .validate_pipeline_config(config)

  bids_root <- normalizePath(bids_root, winslash = "/", mustWork = TRUE)
  input_root <- file.path(bids_root, "derivatives", config$timecourse_derivative)
  output_root <- file.path(bids_root, "derivatives", config$output_derivative_name)

  group_tsv <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.tsv")
  group_json <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.json")
  group_rds <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.rds")
  group_peaks_tsv <- file.path(output_root, "group", "func", "group_desc-peaks_timeseries.tsv")
  group_peaks_json <- file.path(output_root, "group", "func", "group_desc-peaks_timeseries.json")
  group_peaks_rds <- file.path(output_root, "group", "func", "group_desc-peaks_timeseries.rds")
  density_tsv <- file.path(output_root, "group", "func", "group_desc-peakdensity_timeseries.tsv")
  density_json <- file.path(output_root, "group", "func", "group_desc-peakdensity_timeseries.json")
  density_rds <- file.path(output_root, "group", "func", "group_desc-peakdensity_timeseries.rds")
  plot_png <- file.path(output_root, "group", "figures", "group_desc-slidingconn_heatmap.png")
  plot_json <- file.path(output_root, "group", "figures", "group_desc-slidingconn_heatmap.json")
  density_all_png <- file.path(output_root, "group", "figures", "group_desc-peakdensity_all.png")
  density_sub_png <- file.path(output_root, "group", "figures", "group_desc-peakdensity_subgroups.png")

  dir.create(dirname(group_tsv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_png), recursive = TRUE, showWarnings = FALSE)

  group_cached <- file.exists(group_rds) && file.exists(group_peaks_rds)
  density_cached <- file.exists(density_rds)
  heatmap_cached <- file.exists(plot_png)
  density_all_cached <- file.exists(density_all_png)
  density_sub_cached <- file.exists(density_sub_png)

  if (isTRUE(plot_only)) {
    files <- character(0)
    group_table <- .load_existing_group_windows(group_rds = group_rds, group_tsv = group_tsv)
    peak_table <- .load_existing_group_peaks(group_peaks_rds = group_peaks_rds, group_peaks_tsv = group_peaks_tsv)
  } else {
    if (!isTRUE(overwrite) && group_cached) {
      files <- character(0)
      group_table <- readRDS(group_rds)
      peak_table <- readRDS(group_peaks_rds)
    } else {
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

      group_map <- .load_participants_groups(
        bids_root = bids_root,
        subgroup_column = config$subgroup_column
      )

      subject_results <- lapply(files, function(one_file) {
        .process_single_timecourse(
          file_path = one_file,
          output_root = output_root,
          config = config,
          group_map = group_map,
          overwrite = isTRUE(overwrite)
        )
      })

      group_table <- do.call(rbind, lapply(subject_results, function(x) x$sw))
      rownames(group_table) <- NULL
      peak_table <- do.call(rbind, lapply(subject_results, function(x) x$peaks))
      rownames(peak_table) <- NULL

      readr::write_tsv(group_table, group_tsv, na = "n/a")
      saveRDS(group_table, group_rds)

      readr::write_tsv(peak_table, group_peaks_tsv, na = "n/a")
      saveRDS(peak_table, group_peaks_rds)
    }
  }

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

  peaks_meta <- list(
    Description = "Group-level detected connectivity peaks.",
    Columns = list(
      subject = "BIDS subject label",
      session = "BIDS session label",
      group = "Subgroup label from participants.tsv",
      source_file = "Input ROI timecourse TSV",
      run_id = "Run identifier",
      peak_time_seconds = "Peak time in seconds",
      peak_window_center_tp = "Peak center timepoint",
      peak_fisher_z = "Peak Fisher z value",
      prominence = "Peak prominence",
      nearest_event_seconds = "Nearest event marker time in seconds",
      event_distance_seconds = "Absolute distance to nearest event marker"
    ),
    Parameters = config
  )
  jsonlite::write_json(peaks_meta, group_peaks_json, pretty = TRUE, auto_unbox = TRUE)

  if (isTRUE(overwrite) || !heatmap_cached) {
    .plot_group_heatmap(
      group_table = group_table,
      png_path = plot_png,
      event_seconds = config$event_seconds
    )
  }

  if (!isTRUE(overwrite) && density_cached) {
    density_table <- readRDS(density_rds)
  } else {
    density_table <- .compute_peak_density(
      peak_table = peak_table,
      group_table = group_table,
      bin_seconds = as.numeric(config$peak_density_bin_seconds),
      subgroup_levels = config$subgroup_levels
    )
    readr::write_tsv(density_table, density_tsv, na = "n/a")
    saveRDS(density_table, density_rds)
  }

  density_meta <- list(
    Description = "Peak density over time for all subjects and subgroups.",
    Parameters = config
  )
  jsonlite::write_json(density_meta, density_json, pretty = TRUE, auto_unbox = TRUE)

  if (isTRUE(overwrite) || !density_all_cached) {
    .plot_peak_density(density_table = density_table, png_path = density_all_png, subgroup = "ALL")
  }
  if (isTRUE(overwrite) || !density_sub_cached) {
    .plot_peak_density_subgroups(
      density_table = density_table,
      png_path = density_sub_png,
      subgroup_levels = config$subgroup_levels
    )
  }

  jsonlite::write_json(
    list(
      Description = "Heatmap of Fisher-z sliding-window connectivity values.",
      ColorMap = "viridis",
      XAxis = "Window midpoint (seconds)",
      YAxis = "Subject",
      Fill = "Fisher z",
      EventMarkerSeconds = config$event_seconds
    ),
    plot_json,
    pretty = TRUE,
    auto_unbox = TRUE
  )

  list(
    input_root = input_root,
    output_root = output_root,
    plot_only = isTRUE(plot_only),
    overwrite = isTRUE(overwrite),
    n_files_processed = length(files),
    group_tsv = group_tsv,
    group_json = group_json,
    group_rds = group_rds,
    group_peaks_tsv = group_peaks_tsv,
    group_peaks_json = group_peaks_json,
    group_peaks_rds = group_peaks_rds,
    peak_density_tsv = density_tsv,
    peak_density_json = density_json,
    peak_density_rds = density_rds,
    group_plot = plot_png,
    group_plot_json = plot_json,
    peak_density_all_plot = density_all_png,
    peak_density_subgroups_plot = density_sub_png
  )
}

.load_existing_group_windows <- function(group_rds, group_tsv) {
  if (file.exists(group_rds)) {
    return(readRDS(group_rds))
  }

  if (file.exists(group_tsv)) {
    return(readr::read_tsv(group_tsv, show_col_types = FALSE, progress = FALSE))
  }

  stop(
    paste(
      "plot_only = TRUE requires previously saved windows.",
      sprintf("Missing files: %s and %s", group_rds, group_tsv),
      sep = "\n"
    ),
    call. = FALSE
  )
}

.load_existing_group_peaks <- function(group_peaks_rds, group_peaks_tsv) {
  if (file.exists(group_peaks_rds)) {
    return(readRDS(group_peaks_rds))
  }

  if (file.exists(group_peaks_tsv)) {
    return(readr::read_tsv(group_peaks_tsv, show_col_types = FALSE, progress = FALSE))
  }

  stop(
    paste(
      "plot_only = TRUE requires previously saved peak outputs.",
      sprintf("Missing files: %s and %s", group_peaks_rds, group_peaks_tsv),
      sep = "\n"
    ),
    call. = FALSE
  )
}

.process_single_timecourse <- function(file_path, output_root, config, group_map, overwrite = FALSE) {
  entities <- .extract_entities_from_path(file_path)
  roi_x <- config$roi_x
  roi_y <- config$roi_y

  out_dir <- file.path(output_root, entities$subject, entities$session, "func")
  fig_dir <- file.path(output_root, entities$subject, entities$session, "figures")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  out_prefix <- paste(
    c(entities$subject, entities$session, "desc-slidingconn_timeseries"),
    collapse = "_"
  )

  tsv_path <- file.path(out_dir, paste0(out_prefix, ".tsv"))
  json_path <- file.path(out_dir, paste0(out_prefix, ".json"))
  rds_path <- file.path(out_dir, paste0(out_prefix, ".rds"))
  peak_tsv_path <- file.path(out_dir, paste0(entities$subject, "_", entities$session, "_desc-peaks_timeseries.tsv"))
  peak_json_path <- file.path(out_dir, paste0(entities$subject, "_", entities$session, "_desc-peaks_timeseries.json"))
  peak_rds_path <- file.path(out_dir, paste0(entities$subject, "_", entities$session, "_desc-peaks_timeseries.rds"))
  subj_plot_path <- file.path(fig_dir, paste0(entities$subject, "_", entities$session, "_desc-timeseries_peaks.png"))

  sw_cached <- file.exists(rds_path)
  peaks_cached <- file.exists(peak_rds_path)

  subject_group <- .subject_group_from_map(entities$subject, group_map)
  sw <- NULL
  peaks <- NULL

  if (!isTRUE(overwrite) && sw_cached && peaks_cached) {
    sw <- readRDS(rds_path)
    peaks <- readRDS(peak_rds_path)
  }

  needs_plot <- isTRUE(overwrite) || !file.exists(subj_plot_path)
  needs_compute <- is.null(sw) || is.null(peaks)

  x <- NULL
  y <- NULL
  tbl <- NULL

  if (needs_compute || needs_plot) {
    tbl <- readr::read_tsv(
      file_path,
      show_col_types = FALSE,
      progress = FALSE,
      skip_empty_rows = FALSE,
      na = c("NA", "NaN", "n/a")
    )

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
  }

  if (needs_compute) {
    nuisance_signal <- NULL
    if (isTRUE(config$regress_pair_mean_signal)) {
      numeric_cols <- vapply(tbl, is.numeric, logical(1))
      if (!any(numeric_cols)) {
        stop(
          sprintf("No numeric columns available to build nuisance signal in %s", file_path),
          call. = FALSE
        )
      }

      nuisance_signal <- rowMeans(tbl[, numeric_cols, drop = FALSE], na.rm = TRUE)
      nuisance_signal[is.nan(nuisance_signal)] <- NA_real_
    }

    sw <- compute_sliding_connectivity(
      x = x,
      y = y,
      run_boundaries = config$run_boundaries,
      window_size_tp = config$window_size_tp,
      step_size_tp = config$step_size_tp,
      min_points_per_window = config$min_points_per_window,
      regress_pair_mean_signal = config$regress_pair_mean_signal,
      nuisance_signal = nuisance_signal,
      diagnostic_context = sprintf(
        "subject=%s; session=%s; source_file=%s",
        entities$subject,
        entities$session,
        basename(file_path)
      )
    )

    sw$subject <- entities$subject
    sw$session <- entities$session
    sw$group <- subject_group
    sw$source_file <- basename(file_path)
    sw$time_seconds <- (sw$window_center_tp - 1) * config$tr_seconds

    peaks <- .detect_subject_peaks(
      sw = sw,
      event_seconds = config$event_seconds,
      prominence_quantile = as.numeric(config$peak_prominence_quantile),
      min_distance_seconds = as.numeric(config$peak_min_distance_seconds),
      min_height = as.numeric(config$peak_min_height),
      positive_only = isTRUE(config$peak_positive_only),
      event_max_distance_seconds = as.numeric(config$peak_event_max_distance_seconds),
      default_min_distance_seconds = as.numeric(config$window_size_tp) * as.numeric(config$tr_seconds)
    )

    readr::write_tsv(sw, tsv_path, na = "n/a")
    saveRDS(sw, rds_path)
    readr::write_tsv(peaks, peak_tsv_path, na = "n/a")
    saveRDS(peaks, peak_rds_path)
  }

  if (!("group" %in% names(sw))) {
    sw$group <- subject_group
  }

  if (!("group" %in% names(peaks))) {
    peaks$group <- subject_group
  }

  if (needs_plot) {
    if (is.null(x) || is.null(y)) {
      tbl <- readr::read_tsv(
        file_path,
        show_col_types = FALSE,
        progress = FALSE,
        skip_empty_rows = FALSE,
        na = c("NA", "NaN", "n/a")
      )
      x <- as.numeric(tbl[[roi_x]])
      y <- as.numeric(tbl[[roi_y]])
      x[is.na(x)] <- NaN
      y[is.na(y)] <- NaN
    }

    .plot_subject_timeseries_with_peaks(
      x = x,
      y = y,
      sw = sw,
      peaks = peaks,
      roi_x = roi_x,
      roi_y = roi_y,
      event_seconds = config$event_seconds,
      tr_seconds = config$tr_seconds,
      png_path = subj_plot_path,
      subject_label = paste(entities$subject, entities$session)
    )
  }

  sidecar <- list(
    Description = "Sliding-window ROI connectivity values for one subject/session.",
    SourceFile = basename(file_path),
    ROIPair = list(roi_x = roi_x, roi_y = roi_y),
    Parameters = config
  )
  if (isTRUE(overwrite) || !file.exists(json_path)) {
    jsonlite::write_json(sidecar, json_path, pretty = TRUE, auto_unbox = TRUE)
  }

  peak_sidecar <- list(
    Description = "Detected connectivity peaks for one subject/session.",
    SourceFile = basename(file_path),
    ROIPair = list(roi_x = roi_x, roi_y = roi_y),
    Parameters = config
  )
  if (isTRUE(overwrite) || !file.exists(peak_json_path)) {
    jsonlite::write_json(peak_sidecar, peak_json_path, pretty = TRUE, auto_unbox = TRUE)
  }

  list(sw = sw, peaks = peaks)
}

.load_participants_groups <- function(bids_root, subgroup_column) {
  participants_path <- file.path(bids_root, "participants.tsv")
  if (!file.exists(participants_path)) {
    return(setNames(character(0), character(0)))
  }

  participants <- readr::read_tsv(participants_path, show_col_types = FALSE, progress = FALSE)
  if (!("participant_id" %in% names(participants))) {
    return(setNames(character(0), character(0)))
  }
  if (!(subgroup_column %in% names(participants))) {
    return(setNames(character(0), character(0)))
  }

  ids <- as.character(participants$participant_id)
  vals <- as.character(participants[[subgroup_column]])
  setNames(vals, ids)
}

.subject_group_from_map <- function(subject, group_map) {
  if (length(group_map) == 0L) {
    return("UNKNOWN")
  }
  if (!(subject %in% names(group_map))) {
    return("UNKNOWN")
  }
  val <- group_map[[subject]]
  if (is.na(val) || val == "") {
    return("UNKNOWN")
  }
  val
}

.detect_subject_peaks <- function(
    sw,
    event_seconds,
    prominence_quantile,
    min_distance_seconds,
    min_height,
    positive_only,
    event_max_distance_seconds,
    default_min_distance_seconds) {
  run_ids <- unique(as.character(sw$run_id))
  peak_rows <- lapply(run_ids, function(run_id) {
    run_tbl <- sw[sw$run_id == run_id, , drop = FALSE]
    run_tbl <- run_tbl[order(run_tbl$time_seconds), , drop = FALSE]

    valid <- which(is.finite(run_tbl$fisher_z))
    if (length(valid) < 3L) {
      return(data.frame(
        run_id = character(0),
        peak_time_seconds = numeric(0),
        peak_window_center_tp = numeric(0),
        peak_fisher_z = numeric(0),
        prominence = numeric(0),
        nearest_event_seconds = numeric(0),
        event_distance_seconds = numeric(0),
        stringsAsFactors = FALSE
      ))
    }

    seg_starts <- c(valid[1], valid[which(diff(valid) > 1L) + 1L])
    seg_ends <- c(valid[which(diff(valid) > 1L)], valid[length(valid)])

    peak_candidates <- vector("list", length(seg_starts))
    for (s in seq_along(seg_starts)) {
      idx <- seg_starts[[s]]:seg_ends[[s]]
      z <- run_tbl$fisher_z[idx]
      t <- run_tbl$time_seconds[idx]
      tp <- run_tbl$window_center_tp[idx]

      if (length(z) < 3L) {
        peak_candidates[[s]] <- data.frame()
        next
      }

      local_idx <- which(z[2:(length(z) - 1)] > z[1:(length(z) - 2)] &
        z[2:(length(z) - 1)] >= z[3:length(z)]) + 1L

      if (length(local_idx) == 0L) {
        peak_candidates[[s]] <- data.frame()
        next
      }

      prom <- vapply(local_idx, function(i) {
        left_min <- min(z[1:i], na.rm = TRUE)
        right_min <- min(z[i:length(z)], na.rm = TRUE)
        z[i] - max(left_min, right_min)
      }, numeric(1))

      peak_candidates[[s]] <- data.frame(
        run_id = run_id,
        peak_time_seconds = t[local_idx],
        peak_window_center_tp = tp[local_idx],
        peak_fisher_z = z[local_idx],
        prominence = prom,
        stringsAsFactors = FALSE
      )
    }

    cand <- do.call(rbind, peak_candidates)
    if (is.null(cand) || nrow(cand) == 0L) {
      return(data.frame(
        run_id = character(0),
        peak_time_seconds = numeric(0),
        peak_window_center_tp = numeric(0),
        peak_fisher_z = numeric(0),
        prominence = numeric(0),
        nearest_event_seconds = numeric(0),
        event_distance_seconds = numeric(0),
        stringsAsFactors = FALSE
      ))
    }

    prom_thr <- stats::quantile(cand$prominence, probs = prominence_quantile, na.rm = TRUE)
    keep <- cand$prominence >= prom_thr

    if (isTRUE(positive_only)) {
      keep <- keep & cand$peak_fisher_z > 0
    }
    if (is.finite(min_height)) {
      keep <- keep & cand$peak_fisher_z >= min_height
    }

    cand <- cand[keep, , drop = FALSE]
    if (nrow(cand) == 0L) {
      return(data.frame(
        run_id = character(0),
        peak_time_seconds = numeric(0),
        peak_window_center_tp = numeric(0),
        peak_fisher_z = numeric(0),
        prominence = numeric(0),
        nearest_event_seconds = numeric(0),
        event_distance_seconds = numeric(0),
        stringsAsFactors = FALSE
      ))
    }

    min_dist <- if (is.finite(min_distance_seconds)) min_distance_seconds else default_min_distance_seconds
    cand <- cand[order(cand$peak_fisher_z, decreasing = TRUE), , drop = FALSE]
    selected <- cand[0, , drop = FALSE]
    for (i in seq_len(nrow(cand))) {
      if (nrow(selected) == 0L) {
        selected <- cand[i, , drop = FALSE]
      } else {
        d <- abs(selected$peak_time_seconds - cand$peak_time_seconds[[i]])
        if (all(d >= min_dist)) {
          selected <- rbind(selected, cand[i, , drop = FALSE])
        }
      }
    }

    selected <- selected[order(selected$peak_time_seconds), , drop = FALSE]
    if (length(event_seconds) > 0L) {
      nearest_idx <- vapply(selected$peak_time_seconds, function(t0) {
        which.min(abs(event_seconds - t0))
      }, integer(1))
      nearest <- event_seconds[nearest_idx]
      d_event <- abs(selected$peak_time_seconds - nearest)
      if (is.finite(event_max_distance_seconds)) {
        nearest[d_event > event_max_distance_seconds] <- NA_real_
      }
      selected$nearest_event_seconds <- nearest
      selected$event_distance_seconds <- d_event
    } else {
      selected$nearest_event_seconds <- NA_real_
      selected$event_distance_seconds <- NA_real_
    }

    selected
  })

  peaks <- do.call(rbind, peak_rows)
  if (is.null(peaks) || nrow(peaks) == 0L) {
    peaks <- data.frame(
      run_id = character(0),
      peak_time_seconds = numeric(0),
      peak_window_center_tp = numeric(0),
      peak_fisher_z = numeric(0),
      prominence = numeric(0),
      nearest_event_seconds = numeric(0),
      event_distance_seconds = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  peaks$subject <- rep(sw$subject[[1]], nrow(peaks))
  peaks$session <- rep(sw$session[[1]], nrow(peaks))
  peaks$group <- rep(sw$group[[1]], nrow(peaks))
  peaks$source_file <- rep(sw$source_file[[1]], nrow(peaks))
  peaks
}

.plot_subject_timeseries_with_peaks <- function(
    x,
    y,
    sw,
    peaks,
    roi_x,
    roi_y,
    event_seconds,
    tr_seconds,
    png_path,
    subject_label) {
  tp_seconds <- (seq_along(x) - 1) * tr_seconds
  roi_df <- rbind(
    data.frame(
      time_seconds = tp_seconds,
      value = as.numeric(x),
      series = roi_x,
      panel = "ROI activity",
      stringsAsFactors = FALSE
    ),
    data.frame(
      time_seconds = tp_seconds,
      value = as.numeric(y),
      series = roi_y,
      panel = "ROI activity",
      stringsAsFactors = FALSE
    )
  )

  conn_df <- data.frame(
    time_seconds = sw$time_seconds,
    value = sw$fisher_z,
    series = "fisher_z",
    panel = "Sliding connectivity",
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(roi_df, conn_df)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = "time_seconds", y = "value", color = "series")) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.9) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = paste0("Subject ", subject_label, ": ROI activity and connectivity peaks"),
      x = "Time (seconds)",
      y = "Value",
      color = "Series"
    )

  if (length(event_seconds) > 0L) {
    p <- p + ggplot2::geom_vline(
      xintercept = event_seconds,
      color = "red",
      linewidth = 0.15,
      alpha = 0.7
    )
  }

  if (nrow(peaks) > 0L) {
    peak_plot <- data.frame(
      time_seconds = peaks$peak_time_seconds,
      value = peaks$peak_fisher_z,
      panel = "Sliding connectivity",
      stringsAsFactors = FALSE
    )
    p <- p + ggplot2::geom_point(
      data = peak_plot,
      mapping = ggplot2::aes_string(x = "time_seconds", y = "value"),
      color = "black",
      size = 1.5,
      inherit.aes = FALSE
    )
  }

  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 8, dpi = 150)
}

.compute_peak_density <- function(peak_table, group_table, bin_seconds, subgroup_levels) {
  if (!all(c("subject", "run_id", "group") %in% names(group_table))) {
    stop("group_table missing required columns for density computation.", call. = FALSE)
  }

  run_levels <- sort(unique(as.character(group_table$run_id)))
  run_min <- tapply(group_table$time_seconds, group_table$run_id, min, na.rm = TRUE)
  run_max <- tapply(group_table$time_seconds, group_table$run_id, max, na.rm = TRUE)
  run_ranges <- data.frame(
    run_id = run_levels,
    run_min_seconds = as.numeric(run_min[run_levels]),
    run_max_seconds = as.numeric(run_max[run_levels]),
    stringsAsFactors = FALSE
  )

  subject_runs <- unique(group_table[, c("subject", "run_id", "group")])

  out_rows <- lapply(seq_len(nrow(run_ranges)), function(i) {
    run_id <- as.character(run_ranges$run_id[[i]])
    min_t <- as.numeric(run_ranges$run_min_seconds[[i]])
    max_t <- as.numeric(run_ranges$run_max_seconds[[i]])
    breaks <- seq(min_t, max_t + bin_seconds, by = bin_seconds)
    centers <- head(breaks, -1) + bin_seconds / 2

    run_subjects <- subject_runs[subject_runs$run_id == run_id, , drop = FALSE]
    all_n <- nrow(run_subjects)

    run_peaks <- peak_table[peak_table$run_id == run_id, , drop = FALSE]
    all_counts <- hist(run_peaks$peak_time_seconds, breaks = breaks, plot = FALSE)$counts
    all_df <- data.frame(
      run_id = run_id,
      group = "ALL",
      time_bin_seconds = centers,
      peak_count = all_counts,
      subject_count = all_n,
      peak_rate = if (all_n > 0L) all_counts / all_n else NA_real_,
      stringsAsFactors = FALSE
    )

    sub_df <- lapply(subgroup_levels, function(g) {
      g_subjects <- run_subjects[run_subjects$group == g, , drop = FALSE]
      g_n <- nrow(g_subjects)
      g_peaks <- run_peaks[run_peaks$group == g, , drop = FALSE]
      g_counts <- hist(g_peaks$peak_time_seconds, breaks = breaks, plot = FALSE)$counts
      data.frame(
        run_id = run_id,
        group = g,
        time_bin_seconds = centers,
        peak_count = g_counts,
        subject_count = g_n,
        peak_rate = if (g_n > 0L) g_counts / g_n else NA_real_,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, c(list(all_df), sub_df))
  })

  out <- do.call(rbind, out_rows)
  rownames(out) <- NULL
  out
}

.plot_peak_density <- function(density_table, png_path, subgroup = "ALL") {
  tbl <- density_table[density_table$group == subgroup, , drop = FALSE]
  p <- ggplot2::ggplot(tbl, ggplot2::aes_string(x = "time_bin_seconds", y = "peak_rate")) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.4) +
    ggplot2::facet_wrap(~run_id, scales = "free_x", ncol = 1) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Group Peak Density (All Subjects)",
      x = "Time (seconds)",
      y = "Peaks per subject"
    )
  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 7, dpi = 150)
}

.plot_peak_density_subgroups <- function(density_table, png_path, subgroup_levels) {
  tbl <- density_table[density_table$group %in% subgroup_levels, , drop = FALSE]
  p <- ggplot2::ggplot(tbl, ggplot2::aes_string(x = "time_bin_seconds", y = "peak_rate", color = "group")) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::facet_grid(group ~ run_id, scales = "free_x") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Group Peak Density by Subgroup",
      x = "Time (seconds)",
      y = "Peaks per subject",
      color = "Group"
    )
  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 8, dpi = 150)
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
    "event_seconds",
    "peak_prominence_quantile",
    "peak_min_distance_seconds",
    "peak_min_height",
    "peak_positive_only",
    "peak_event_max_distance_seconds",
    "subgroup_column",
    "subgroup_levels",
    "peak_density_bin_seconds",
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

  if (!is.numeric(config$event_seconds)) {
    stop("config$event_seconds must be numeric.", call. = FALSE)
  }

  if (!is.numeric(config$peak_prominence_quantile) ||
      config$peak_prominence_quantile < 0 ||
      config$peak_prominence_quantile > 1) {
    stop("config$peak_prominence_quantile must be numeric in [0, 1].", call. = FALSE)
  }

  if (!is.numeric(config$peak_min_distance_seconds)) {
    stop("config$peak_min_distance_seconds must be numeric.", call. = FALSE)
  }

  if (!is.numeric(config$peak_min_height)) {
    stop("config$peak_min_height must be numeric.", call. = FALSE)
  }

  if (!is.logical(config$peak_positive_only) || length(config$peak_positive_only) != 1L) {
    stop("config$peak_positive_only must be TRUE/FALSE.", call. = FALSE)
  }

  if (!is.numeric(config$peak_event_max_distance_seconds)) {
    stop("config$peak_event_max_distance_seconds must be numeric.", call. = FALSE)
  }

  if (!is.character(config$subgroup_column) || length(config$subgroup_column) != 1L) {
    stop("config$subgroup_column must be a single character string.", call. = FALSE)
  }

  if (!is.character(config$subgroup_levels) || length(config$subgroup_levels) < 1L) {
    stop("config$subgroup_levels must be a character vector.", call. = FALSE)
  }

  if (!is.numeric(config$peak_density_bin_seconds) || config$peak_density_bin_seconds <= 0) {
    stop("config$peak_density_bin_seconds must be > 0.", call. = FALSE)
  }
}

.plot_group_heatmap <- function(group_table, png_path, event_seconds = numeric(0)) {
  plot_tbl <- group_table
  plot_tbl$subject_session <- paste(plot_tbl$subject, plot_tbl$session, sep = "_")
  plot_tbl <- plot_tbl[stats::complete.cases(plot_tbl$fisher_z), , drop = FALSE]
  event_seconds <- as.numeric(event_seconds)
  event_seconds <- event_seconds[is.finite(event_seconds)]

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

  run_levels <- sort(unique(as.character(plot_tbl$run_id)))
  run_min <- tapply(plot_tbl$time_seconds, plot_tbl$run_id, min, na.rm = TRUE)
  run_max <- tapply(plot_tbl$time_seconds, plot_tbl$run_id, max, na.rm = TRUE)

  run_ranges <- data.frame(
    run_id = run_levels,
    run_min_seconds = as.numeric(run_min[run_levels]),
    run_max_seconds = as.numeric(run_max[run_levels]),
    stringsAsFactors = FALSE
  )

  if (length(event_seconds) > 0L) {
    event_df <- merge(
      run_ranges,
      data.frame(event_seconds = event_seconds),
      by = NULL
    )
    event_df <- event_df[
      event_df$event_seconds >= event_df$run_min_seconds &
        event_df$event_seconds <= event_df$run_max_seconds,
      c("run_id", "event_seconds"),
      drop = FALSE
    ]
  } else {
    event_df <- data.frame(run_id = character(0), event_seconds = numeric(0))
  }

  p <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes_string(
      x = "time_seconds",
      y = "subject_session",
      fill = "fisher_z"
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(
      data = event_df,
      mapping = ggplot2::aes_string(xintercept = "event_seconds"),
      inherit.aes = FALSE,
      color = "red",
      linewidth = 0.2,
      alpha = 0.9
    ) +
    ggplot2::facet_wrap(~run_id, scales = "free_x", ncol = 1) +
    viridis::scale_fill_viridis(option = "viridis", na.value = "grey80") +
    ggplot2::labs(
      title = "Sliding-window Fisher-z Connectivity",
      x = "Time (seconds)",
      y = "",
      fill = "Fisher z"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    filename = png_path,
    plot = p,
    width = 12,
    height = 8,
    dpi = 150
  )
}
