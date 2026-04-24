test_that("zone detection finds sustained elevated regions after smoothing", {
  sw <- data.frame(
    run_id = rep("run-1", 12),
    window_center_tp = seq_len(12),
    fisher_z = c(0.1, 0.2, 0.75, 0.8, 0.82, 0.78, 0.2, 0.1, 0.7, 0.74, 0.72, 0.1),
    subject = rep("sub-101", 12),
    session = rep("ses-01", 12),
    group = rep("YA", 12),
    source_file = rep("dummy.tsv", 12),
    time_seconds = seq(0, by = 2, length.out = 12),
    stringsAsFactors = FALSE
  )

  out <- syncfmri:::.smooth_and_detect_zones(
    sw = sw,
    event_seconds = c(8, 20),
    smoothing_window = 3L,
    threshold_quantile = 0.75,
    min_consecutive_windows = 2L,
    merge_gap_windows = 1L,
    event_max_distance_seconds = Inf,
    event_max_distance_seconds = Inf
  )

  expect_true(nrow(out$zones) >= 1)
  expect_true(all(out$zones$zone_peak_fisher_z > 0))
  expect_true(all(out$zones$zone_n_windows >= 2))
  expect_true(any(out$sw$in_zone %in% TRUE))
})

test_that("zone detection handles no-zone case without error", {
  sw <- data.frame(
    run_id = rep("run-1", 6),
    window_center_tp = seq_len(6),
    fisher_z = rep(-0.2, 6),
    subject = rep("sub-213", 6),
    session = rep("ses-01", 6),
    group = rep("OA", 6),
    source_file = rep("dummy.tsv", 6),
    time_seconds = seq(0, by = 2, length.out = 6),
    stringsAsFactors = FALSE
  )

  out <- syncfmri:::.smooth_and_detect_zones(
    sw = sw,
    event_seconds = numeric(0),
    smoothing_window = 3L,
    threshold_quantile = 0.75,
    min_consecutive_windows = 2L,
    merge_gap_windows = 1L,
    event_max_distance_seconds = Inf,
    event_max_distance_seconds = Inf
  )

  expect_equal(nrow(out$zones), 0)
  expect_true(all(c("smoothed_fisher_z", "zone_id", "in_zone") %in% names(out$sw)))
})
