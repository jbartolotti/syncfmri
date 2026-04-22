#!/usr/bin/env Rscript

# Minimal wrapper script for running syncfmri on a BIDS dataset.
# Copy this script to your execution location and update only the paths/settings.

library(syncfmri)

bids_root <- "P:/IRB_STUDY00149390_A015/MR_Data/Connectivity/natfMRI/A015_BIDS"

cfg <- syncfmri_default_config()
cfg$timecourse_derivative <- "gimmefMRI"
cfg$timecourse_file_regex <- "\\.tsv$"
cfg$roi_x <- "LR_postJHipp_200"
cfg$roi_y <- "LR_mPFC_200"
cfg$run_boundaries <- data.frame(
  run_id = c("run-1", "run-2", "run-3"),
  start_tp = c(1L, 360L, 688L),
  end_tp = c(359L, 687L, 1012L)
)
cfg$tr_seconds <- 2
cfg$window_size_tp <- 12L
cfg$step_size_tp <- 1L
cfg$min_points_per_window <- 8L
cfg$regress_pair_mean_signal <- TRUE
cfg$output_derivative_name <- "syncfmri"

result <- syncfmri_run_pipeline(bids_root = bids_root, config = cfg)
print(result)
