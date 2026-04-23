test_that("peak detection finds local maxima with prominence threshold", {
  sw <- data.frame(
    run_id = rep("run-1", 9),
    window_center_tp = seq_len(9),
    fisher_z = c(0.1, 0.4, 0.2, 0.1, 0.8, 0.2, 0.1, 0.5, 0.1),
    subject = rep("sub-101", 9),
    session = rep("ses-01", 9),
    group = rep("YA", 9),
    source_file = rep("dummy.tsv", 9),
    time_seconds = seq(0, by = 2, length.out = 9),
    stringsAsFactors = FALSE
  )

  peaks <- syncfmri:::.detect_subject_peaks(
    sw = sw,
    event_seconds = c(8, 20),
    prominence_quantile = 0.5,
    min_distance_seconds = 4,
    min_height = NA_real_,
    positive_only = TRUE,
    event_max_distance_seconds = Inf,
    default_min_distance_seconds = 4
  )

  expect_true(nrow(peaks) >= 1)
  expect_true(all(peaks$peak_fisher_z > 0))
  expect_true(all(is.finite(peaks$prominence)))
  expect_true(all(is.finite(peaks$event_distance_seconds)))
})
