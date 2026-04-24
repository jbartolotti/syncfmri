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
    zone_smoothing_window = 3L,
    zone_threshold_quantile = 0.75,
    zone_min_consecutive_windows = 2L,
    zone_merge_gap_windows = 1L,
    zone_event_max_distance_seconds = Inf,
    zone_density_bin_seconds = 4,
    subgroup_column = "group",
    subgroup_levels = c("YA", "OA"),
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
#' @param overwrite_windows Logical; if TRUE, force recomputation of sliding
#'   window connectivity.
#' @param overwrite_zones Logical; if TRUE, force recomputation of smoothing
#'   and zone identification (can reuse cached sliding windows).
#'
#' @return A list containing subject and group output paths.
#' @export
syncfmri_run_pipeline <- function(
    bids_root,
    config = syncfmri_default_config(),
    plot_only = FALSE,
    overwrite = FALSE,
    overwrite_windows = FALSE,
    overwrite_zones = FALSE) {
  .validate_pipeline_config(config)

  overwrite_windows <- isTRUE(overwrite) || isTRUE(overwrite_windows)
  overwrite_zones <- isTRUE(overwrite) || isTRUE(overwrite_zones) || isTRUE(overwrite_windows)

  bids_root <- normalizePath(bids_root, winslash = "/", mustWork = TRUE)
  input_root <- file.path(bids_root, "derivatives", config$timecourse_derivative)
  output_root <- file.path(bids_root, "derivatives", config$output_derivative_name)

  group_tsv <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.tsv")
  group_json <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.json")
  group_rds <- file.path(output_root, "group", "func", "group_desc-slidingconn_timeseries.rds")
  group_zones_tsv <- file.path(output_root, "group", "func", "group_desc-zones_timeseries.tsv")
  group_zones_json <- file.path(output_root, "group", "func", "group_desc-zones_timeseries.json")
  group_zones_rds <- file.path(output_root, "group", "func", "group_desc-zones_timeseries.rds")
  occupancy_tsv <- file.path(output_root, "group", "func", "group_desc-zoneoccupancy_timeseries.tsv")
  occupancy_json <- file.path(output_root, "group", "func", "group_desc-zoneoccupancy_timeseries.json")
  occupancy_rds <- file.path(output_root, "group", "func", "group_desc-zoneoccupancy_timeseries.rds")
  plot_png <- file.path(output_root, "group", "figures", "group_desc-slidingconn_heatmap.png")
  plot_json <- file.path(output_root, "group", "figures", "group_desc-slidingconn_heatmap.json")
  occupancy_sub_png <- file.path(output_root, "group", "figures", "group_desc-zoneoccupancy_subgroups.png")
  zone_heatmap_png <- file.path(output_root, "group", "figures", "group_desc-zones_heatmap.png")

  dir.create(dirname(group_tsv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(plot_png), recursive = TRUE, showWarnings = FALSE)

  group_cached <- file.exists(group_rds) && file.exists(group_zones_rds)
  occupancy_cached <- file.exists(occupancy_rds)
  heatmap_cached <- file.exists(plot_png)
  occupancy_sub_cached <- file.exists(occupancy_sub_png)
  zone_heatmap_cached <- file.exists(zone_heatmap_png)

  if (isTRUE(plot_only)) {
    files <- character(0)
    group_table <- .load_existing_group_windows(group_rds = group_rds, group_tsv = group_tsv)
    zone_table <- .load_existing_group_zones(group_zones_rds = group_zones_rds, group_zones_tsv = group_zones_tsv)
  } else {
    if (!overwrite_windows && !overwrite_zones && group_cached) {
      files <- character(0)
      group_table <- readRDS(group_rds)
      zone_table <- readRDS(group_zones_rds)
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
          overwrite_windows = overwrite_windows,
          overwrite_zones = overwrite_zones
        )
      })

      group_table <- do.call(rbind, lapply(subject_results, function(x) x$sw))
      rownames(group_table) <- NULL
      zone_table <- do.call(rbind, lapply(subject_results, function(x) x$zones))
      rownames(zone_table) <- NULL

      readr::write_tsv(group_table, group_tsv, na = "n/a")
      saveRDS(group_table, group_rds)

      readr::write_tsv(zone_table, group_zones_tsv, na = "n/a")
      saveRDS(zone_table, group_zones_rds)
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

  zones_meta <- list(
    Description = "Group-level detected sustained connectivity zones.",
    Columns = list(
      subject = "BIDS subject label",
      session = "BIDS session label",
      group = "Subgroup label from participants.tsv",
      source_file = "Input ROI timecourse TSV",
      run_id = "Run identifier",
      zone_id = "Zone identifier within subject/session",
      zone_start_seconds = "Zone start time in seconds",
      zone_end_seconds = "Zone end time in seconds",
      zone_duration_seconds = "Zone duration in seconds",
      zone_peak_seconds = "Time of maximum smoothed Fisher z within zone",
      zone_peak_fisher_z = "Peak smoothed Fisher z within zone",
      zone_mean_fisher_z = "Mean smoothed Fisher z across zone",
      zone_auc = "Area under smoothed Fisher z within zone",
      zone_n_windows = "Number of windows in zone",
      nearest_event_seconds = "Nearest event marker time in seconds",
      event_distance_seconds = "Absolute distance to nearest event marker"
    ),
    Parameters = config
  )
  jsonlite::write_json(zones_meta, group_zones_json, pretty = TRUE, auto_unbox = TRUE)

  if (overwrite_windows || !heatmap_cached) {
    .plot_group_heatmap(
      group_table = group_table,
      png_path = plot_png,
      event_seconds = config$event_seconds
    )
  }

  if (!overwrite_zones && occupancy_cached) {
    occupancy_table <- readRDS(occupancy_rds)
  } else {
    occupancy_table <- .compute_zone_occupancy(
      zone_table = zone_table,
      group_table = group_table,
      bin_seconds = as.numeric(config$zone_density_bin_seconds),
      subgroup_levels = config$subgroup_levels
    )
    readr::write_tsv(occupancy_table, occupancy_tsv, na = "n/a")
    saveRDS(occupancy_table, occupancy_rds)
  }

  occupancy_meta <- list(
    Description = "Zone occupancy over time for all subjects and subgroups.",
    Parameters = config
  )
  jsonlite::write_json(occupancy_meta, occupancy_json, pretty = TRUE, auto_unbox = TRUE)

  if (overwrite_zones || !occupancy_sub_cached) {
    .plot_zone_occupancy_subgroups(
      occupancy_table = occupancy_table,
      png_path = occupancy_sub_png,
      event_seconds = config$event_seconds,
      subgroup_levels = config$subgroup_levels
    )
  }
  if (overwrite_zones || !zone_heatmap_cached) {
    .plot_zone_heatmap(
      group_table = group_table,
      png_path = zone_heatmap_png,
      event_seconds = config$event_seconds
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
    overwrite_windows = overwrite_windows,
    overwrite_zones = overwrite_zones,
    n_files_processed = length(files),
    group_tsv = group_tsv,
    group_json = group_json,
    group_rds = group_rds,
    group_zones_tsv = group_zones_tsv,
    group_zones_json = group_zones_json,
    group_zones_rds = group_zones_rds,
    zone_occupancy_tsv = occupancy_tsv,
    zone_occupancy_json = occupancy_json,
    zone_occupancy_rds = occupancy_rds,
    group_plot = plot_png,
    group_plot_json = plot_json,
    zone_occupancy_subgroups_plot = occupancy_sub_png,
    zone_heatmap_plot = zone_heatmap_png
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

.load_existing_group_zones <- function(group_zones_rds, group_zones_tsv) {
  if (file.exists(group_zones_rds)) {
    return(readRDS(group_zones_rds))
  }

  if (file.exists(group_zones_tsv)) {
    return(readr::read_tsv(group_zones_tsv, show_col_types = FALSE, progress = FALSE))
  }

  stop(
    paste(
      "plot_only = TRUE requires previously saved zone outputs.",
      sprintf("Missing files: %s and %s", group_zones_rds, group_zones_tsv),
      sep = "\n"
    ),
    call. = FALSE
  )
}

.process_single_timecourse <- function(
  file_path,
  output_root,
  config,
  group_map,
  overwrite_windows = FALSE,
  overwrite_zones = FALSE) {
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
  zone_tsv_path <- file.path(out_dir, paste0(entities$subject, "_", entities$session, "_desc-zones_timeseries.tsv"))
  zone_json_path <- file.path(out_dir, paste0(entities$subject, "_", entities$session, "_desc-zones_timeseries.json"))
  zone_rds_path <- file.path(out_dir, paste0(entities$subject, "_", entities$session, "_desc-zones_timeseries.rds"))
  subj_plot_path <- file.path(fig_dir, paste0(entities$subject, "_", entities$session, "_desc-timeseries_zones.png"))

  sw_cached <- file.exists(rds_path)
  zones_cached <- file.exists(zone_rds_path)

  subject_group <- .subject_group_from_map(entities$subject, group_map)
  sw <- NULL
  zones <- NULL

  if (!overwrite_windows && sw_cached) {
    sw <- readRDS(rds_path)
  }
  if (!overwrite_zones && zones_cached) {
    zones <- readRDS(zone_rds_path)
  }

  if (overwrite_windows) {
    sw <- NULL
  }
  if (overwrite_zones) {
    zones <- NULL
  }

  needs_compute_windows <- is.null(sw)
  needs_compute_zones <- is.null(zones)
  if (needs_compute_windows) {
    needs_compute_zones <- TRUE
  }

  needs_plot <- overwrite_zones || overwrite_windows || !file.exists(subj_plot_path)

  x <- NULL
  y <- NULL
  tbl <- NULL

  if (needs_compute_windows || needs_plot) {
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

  if (needs_compute_windows) {
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
    readr::write_tsv(sw, tsv_path, na = "n/a")
    saveRDS(sw, rds_path)
  }

  if (needs_compute_zones) {
    smooth_result <- .smooth_and_detect_zones(
      sw = sw,
      smoothing_window = as.integer(config$zone_smoothing_window),
      threshold_quantile = as.numeric(config$zone_threshold_quantile),
      min_consecutive_windows = as.integer(config$zone_min_consecutive_windows),
      merge_gap_windows = as.integer(config$zone_merge_gap_windows),
      event_seconds = config$event_seconds,
      event_max_distance_seconds = as.numeric(config$zone_event_max_distance_seconds)
    )
    sw$smoothed_fisher_z <- smooth_result$sw$smoothed_fisher_z
    sw$zone_id <- smooth_result$sw$zone_id
    sw$in_zone <- smooth_result$sw$in_zone

    zones <- smooth_result$zones
    zones$subject <- rep(entities$subject, nrow(zones))
    zones$session <- rep(entities$session, nrow(zones))
    zones$group <- rep(subject_group, nrow(zones))
    zones$source_file <- rep(basename(file_path), nrow(zones))

    readr::write_tsv(sw, tsv_path, na = "n/a")
    saveRDS(sw, rds_path)
    readr::write_tsv(zones, zone_tsv_path, na = "n/a")
    saveRDS(zones, zone_rds_path)
  }

  if (!("group" %in% names(sw))) {
    sw$group <- subject_group
  }
  if (!("smoothed_fisher_z" %in% names(sw))) {
    stop("Cached sliding-window file is missing smoothed_fisher_z; rerun with overwrite = TRUE.", call. = FALSE)
  }
  if (!("zone_id" %in% names(sw))) {
    stop("Cached sliding-window file is missing zone_id; rerun with overwrite = TRUE.", call. = FALSE)
  }
  if (!("in_zone" %in% names(sw))) {
    stop("Cached sliding-window file is missing in_zone; rerun with overwrite = TRUE.", call. = FALSE)
  }

  if (!("group" %in% names(zones))) {
    zones$group <- rep(subject_group, nrow(zones))
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
      zones = zones,
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
  if (overwrite_windows || !file.exists(json_path)) {
    jsonlite::write_json(sidecar, json_path, pretty = TRUE, auto_unbox = TRUE)
  }

  zone_sidecar <- list(
    Description = "Detected sustained connectivity zones for one subject/session.",
    SourceFile = basename(file_path),
    ROIPair = list(roi_x = roi_x, roi_y = roi_y),
    Parameters = config
  )
  if (overwrite_zones || !file.exists(zone_json_path)) {
    jsonlite::write_json(zone_sidecar, zone_json_path, pretty = TRUE, auto_unbox = TRUE)
  }

  list(sw = sw, zones = zones)
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

.smooth_and_detect_zones <- function(
    sw,
    smoothing_window,
    threshold_quantile,
    min_consecutive_windows,
    merge_gap_windows,
    event_seconds,
    event_max_distance_seconds) {
  sw <- sw[order(sw$run_id, sw$time_seconds), , drop = FALSE]
  sw$smoothed_fisher_z <- NA_real_
  sw$zone_id <- NA_character_
  sw$in_zone <- FALSE

  zone_rows <- list()
  zone_counter <- 0L

  for (run_id in unique(as.character(sw$run_id))) {
    idx_run <- which(sw$run_id == run_id)
    run_tbl <- sw[idx_run, , drop = FALSE]
    valid <- is.finite(run_tbl$fisher_z)

    if (sum(valid) > 0L) {
      valid_idx <- which(valid)
      seg_starts <- c(valid_idx[1], valid_idx[which(diff(valid_idx) > 1L) + 1L])
      seg_ends <- c(valid_idx[which(diff(valid_idx) > 1L)], valid_idx[length(valid_idx)])

      for (s in seq_along(seg_starts)) {
        seg_idx <- seg_starts[[s]]:seg_ends[[s]]
        z <- run_tbl$fisher_z[seg_idx]
        smooth_z <- .rolling_mean_partial(z, k = smoothing_window)
        run_tbl$smoothed_fisher_z[seg_idx] <- smooth_z

        finite_seg <- is.finite(smooth_z)
        if (sum(finite_seg) < min_consecutive_windows) {
          next
        }

        thr <- stats::quantile(smooth_z[finite_seg], probs = threshold_quantile, na.rm = TRUE)
        elevated <- smooth_z >= thr
        elevated[!finite_seg] <- FALSE
        elevated <- .merge_short_false_gaps(elevated, max_gap = merge_gap_windows)
        elevated <- .drop_short_true_runs(elevated, min_len = min_consecutive_windows)

        if (any(elevated)) {
          r <- rle(elevated)
          ends <- cumsum(r$lengths)
          starts <- ends - r$lengths + 1L
          zone_runs <- which(r$values)

          for (zr in zone_runs) {
            zone_counter <- zone_counter + 1L
            local_idx <- starts[[zr]]:ends[[zr]]
            global_idx <- idx_run[seg_idx[local_idx]]
            sw$in_zone[global_idx] <- TRUE
            sw$zone_id[global_idx] <- paste0("zone-", zone_counter)

            zone_times <- sw$time_seconds[global_idx]
            zone_vals <- sw$smoothed_fisher_z[global_idx]
            peak_idx <- which.max(zone_vals)
            zone_center <- mean(range(zone_times))

            nearest_event <- NA_real_
            event_distance <- NA_real_
            if (length(event_seconds) > 0L) {
              nearest_event <- event_seconds[[which.min(abs(event_seconds - zone_center))]]
              event_distance <- abs(nearest_event - zone_center)
              if (is.finite(event_max_distance_seconds) && event_distance > event_max_distance_seconds) {
                nearest_event <- NA_real_
              }
            }

            zone_rows[[length(zone_rows) + 1L]] <- data.frame(
              run_id = run_id,
              zone_id = paste0("zone-", zone_counter),
              zone_start_seconds = min(zone_times),
              zone_end_seconds = max(zone_times),
              zone_duration_seconds = max(zone_times) - min(zone_times) + median(diff(zone_times)),
              zone_peak_seconds = zone_times[[peak_idx]],
              zone_peak_fisher_z = zone_vals[[peak_idx]],
              zone_mean_fisher_z = mean(zone_vals, na.rm = TRUE),
              zone_auc = sum(zone_vals, na.rm = TRUE) * median(diff(zone_times)),
              zone_n_windows = length(zone_vals),
              nearest_event_seconds = nearest_event,
              event_distance_seconds = event_distance,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    sw$smoothed_fisher_z[idx_run] <- run_tbl$smoothed_fisher_z
  }

  zones <- if (length(zone_rows) > 0L) do.call(rbind, zone_rows) else data.frame(
    run_id = character(0),
    zone_id = character(0),
    zone_start_seconds = numeric(0),
    zone_end_seconds = numeric(0),
    zone_duration_seconds = numeric(0),
    zone_peak_seconds = numeric(0),
    zone_peak_fisher_z = numeric(0),
    zone_mean_fisher_z = numeric(0),
    zone_auc = numeric(0),
    zone_n_windows = integer(0),
    nearest_event_seconds = numeric(0),
    event_distance_seconds = numeric(0),
    stringsAsFactors = FALSE
  )

  list(sw = sw, zones = zones)
}

.rolling_mean_partial <- function(x, k) {
  n <- length(x)
  half <- floor(k / 2)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    out[[i]] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}

.merge_short_false_gaps <- function(x, max_gap) {
  if (max_gap < 1L || length(x) == 0L) {
    return(x)
  }
  r <- rle(x)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  for (i in seq_along(r$values)) {
    if (!r$values[[i]] && r$lengths[[i]] <= max_gap && i > 1L && i < length(r$values)) {
      if (isTRUE(r$values[[i - 1L]]) && isTRUE(r$values[[i + 1L]])) {
        x[starts[[i]]:ends[[i]]] <- TRUE
      }
    }
  }
  x
}

.drop_short_true_runs <- function(x, min_len) {
  if (min_len <= 1L || length(x) == 0L) {
    return(x)
  }
  r <- rle(x)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  for (i in seq_along(r$values)) {
    if (isTRUE(r$values[[i]]) && r$lengths[[i]] < min_len) {
      x[starts[[i]]:ends[[i]]] <- FALSE
    }
  }
  x
}

.plot_subject_timeseries_with_peaks <- function(
    x,
    y,
    sw,
    zones,
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
    series = "fisher_z_raw",
    panel = "Sliding connectivity",
    stringsAsFactors = FALSE
  )
  smooth_df <- data.frame(
    time_seconds = sw$time_seconds,
    value = sw$smoothed_fisher_z,
    series = "fisher_z_smoothed",
    panel = "Sliding connectivity",
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(roi_df, conn_df, smooth_df)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = "time_seconds", y = "value", color = "series")) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.9) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      title = paste0("Subject ", subject_label, ": ROI activity and connectivity zones"),
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

  if (nrow(zones) > 0L) {
    zone_rects <- data.frame(
      xmin = zones$zone_start_seconds,
      xmax = zones$zone_end_seconds,
      ymin = -Inf,
      ymax = Inf,
      panel = "Sliding connectivity",
      stringsAsFactors = FALSE
    )
    p <- p + ggplot2::geom_rect(
      data = zone_rects,
      mapping = ggplot2::aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax"),
      fill = "gold",
      alpha = 0.15,
      inherit.aes = FALSE
    )
  }

  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 8, dpi = 150)
}

.compute_zone_occupancy <- function(zone_table, group_table, bin_seconds, subgroup_levels) {
  if (!all(c("subject", "run_id", "group") %in% names(group_table))) {
    stop("group_table missing required columns for zone occupancy computation.", call. = FALSE)
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
    run_windows <- group_table[group_table$run_id == run_id & group_table$in_zone %in% TRUE, , drop = FALSE]
    all_counts <- .count_subjects_in_bins(run_windows, breaks)
    all_df <- data.frame(
      run_id = run_id,
      group = "ALL",
      time_bin_seconds = centers,
      zone_count = all_counts,
      subject_count = all_n,
      zone_proportion = if (all_n > 0L) all_counts / all_n else NA_real_,
      stringsAsFactors = FALSE
    )

    sub_df <- lapply(subgroup_levels, function(g) {
      g_subjects <- run_subjects[run_subjects$group == g, , drop = FALSE]
      g_n <- nrow(g_subjects)
      g_windows <- group_table[group_table$run_id == run_id & group_table$group == g & group_table$in_zone %in% TRUE, , drop = FALSE]
      g_counts <- .count_subjects_in_bins(g_windows, breaks)
      data.frame(
        run_id = run_id,
        group = g,
        time_bin_seconds = centers,
        zone_count = g_counts,
        subject_count = g_n,
        zone_proportion = if (g_n > 0L) g_counts / g_n else NA_real_,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, c(list(all_df), sub_df))
  })

  out <- do.call(rbind, out_rows)
  rownames(out) <- NULL
  out
}

.count_subjects_in_bins <- function(windows_tbl, breaks) {
  if (nrow(windows_tbl) == 0L) {
    return(integer(length(breaks) - 1L))
  }

  bin_idx <- findInterval(windows_tbl$time_seconds, vec = breaks, rightmost.closed = TRUE, all.inside = FALSE)
  keep <- bin_idx >= 1L & bin_idx < length(breaks)
  if (!any(keep)) {
    return(integer(length(breaks) - 1L))
  }

  uniq <- unique(data.frame(subject = windows_tbl$subject[keep], bin = bin_idx[keep], stringsAsFactors = FALSE))
  tabulate(uniq$bin, nbins = length(breaks) - 1L)
}

.plot_zone_occupancy_subgroups <- function(occupancy_table, png_path, event_seconds, subgroup_levels) {
  tbl <- occupancy_table[occupancy_table$group %in% subgroup_levels, , drop = FALSE]
  event_df <- .build_run_specific_event_df(tbl[, c("run_id", "time_bin_seconds"), drop = FALSE], "time_bin_seconds", event_seconds)
  p <- ggplot2::ggplot(tbl, ggplot2::aes_string(x = "time_bin_seconds", y = "zone_proportion", color = "group")) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::facet_grid(group ~ run_id, scales = "free_x") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Zone Occupancy by Age Group",
      x = "Time (seconds)",
      y = "Proportion in zone",
      color = "Group"
    )
  if (nrow(event_df) > 0L) {
    p <- p + ggplot2::geom_vline(
      data = event_df,
      mapping = ggplot2::aes_string(xintercept = "event_seconds"),
      inherit.aes = FALSE,
      color = "red",
      linewidth = 0.15,
      alpha = 0.7
    )
  }
  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 8, dpi = 150)
}

.plot_zone_heatmap <- function(group_table, png_path, event_seconds = numeric(0)) {
  plot_tbl <- group_table
  plot_tbl$subject_session <- paste(plot_tbl$subject, plot_tbl$session, sep = "_")
  plot_tbl$zone_value <- ifelse(plot_tbl$in_zone %in% TRUE, plot_tbl$smoothed_fisher_z, NA_real_)
  plot_tbl <- plot_tbl[stats::complete.cases(plot_tbl$zone_value), , drop = FALSE]
  event_df <- .build_run_specific_event_df(plot_tbl[, c("run_id", "time_seconds"), drop = FALSE], "time_seconds", event_seconds)

  if (nrow(plot_tbl) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No zone windows to plot") +
      ggplot2::theme_void()
    ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 8, dpi = 150)
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes_string(x = "time_seconds", y = "subject_session", fill = "zone_value")
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~run_id, scales = "free_x", ncol = 1) +
    viridis::scale_fill_viridis(option = "viridis", na.value = "white") +
    ggplot2::labs(
      title = "Group Zone Heatmap",
      x = "Time (seconds)",
      y = "",
      fill = "Smoothed Fisher z"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  if (nrow(event_df) > 0L) {
    p <- p + ggplot2::geom_vline(
      data = event_df,
      mapping = ggplot2::aes_string(xintercept = "event_seconds"),
      inherit.aes = FALSE,
      color = "red",
      linewidth = 0.15,
      alpha = 0.7
    )
  }

  ggplot2::ggsave(filename = png_path, plot = p, width = 12, height = 8, dpi = 150)
}

.build_run_specific_event_df <- function(tbl, time_col, event_seconds) {
  event_seconds <- as.numeric(event_seconds)
  event_seconds <- event_seconds[is.finite(event_seconds)]
  if (length(event_seconds) == 0L || nrow(tbl) == 0L) {
    return(data.frame(run_id = character(0), event_seconds = numeric(0), stringsAsFactors = FALSE))
  }

  run_levels <- sort(unique(as.character(tbl$run_id)))
  time_vals <- tbl[[time_col]]
  run_min <- tapply(time_vals, tbl$run_id, min, na.rm = TRUE)
  run_max <- tapply(time_vals, tbl$run_id, max, na.rm = TRUE)
  run_ranges <- data.frame(
    run_id = run_levels,
    run_min_seconds = as.numeric(run_min[run_levels]),
    run_max_seconds = as.numeric(run_max[run_levels]),
    stringsAsFactors = FALSE
  )

  event_df <- merge(run_ranges, data.frame(event_seconds = event_seconds), by = NULL)
  event_df <- event_df[
    event_df$event_seconds >= event_df$run_min_seconds &
      event_df$event_seconds <= event_df$run_max_seconds,
    c("run_id", "event_seconds"),
    drop = FALSE
  ]
  event_df
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
    "zone_smoothing_window",
    "zone_threshold_quantile",
    "zone_min_consecutive_windows",
    "zone_merge_gap_windows",
    "zone_event_max_distance_seconds",
    "zone_density_bin_seconds",
    "subgroup_column",
    "subgroup_levels",
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

  if (as.integer(config$zone_smoothing_window) < 1L) {
    stop("config$zone_smoothing_window must be >= 1.", call. = FALSE)
  }

  if (!is.numeric(config$zone_threshold_quantile) || config$zone_threshold_quantile < 0 || config$zone_threshold_quantile > 1) {
    stop("config$zone_threshold_quantile must be numeric in [0, 1].", call. = FALSE)
  }

  if (as.integer(config$zone_min_consecutive_windows) < 1L) {
    stop("config$zone_min_consecutive_windows must be >= 1.", call. = FALSE)
  }

  if (as.integer(config$zone_merge_gap_windows) < 0L) {
    stop("config$zone_merge_gap_windows must be >= 0.", call. = FALSE)
  }

  if (!is.numeric(config$zone_event_max_distance_seconds)) {
    stop("config$zone_event_max_distance_seconds must be numeric.", call. = FALSE)
  }

  if (!is.numeric(config$zone_density_bin_seconds) || config$zone_density_bin_seconds <= 0) {
    stop("config$zone_density_bin_seconds must be > 0.", call. = FALSE)
  }

  if (!is.character(config$subgroup_column) || length(config$subgroup_column) != 1L) {
    stop("config$subgroup_column must be a single character string.", call. = FALSE)
  }

  if (!is.character(config$subgroup_levels) || length(config$subgroup_levels) < 1L) {
    stop("config$subgroup_levels must be a character vector.", call. = FALSE)
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
