test_that("sliding windows do not cross run boundaries", {
  x <- rep(seq_len(20), 2)
  y <- x

  runs <- data.frame(
    run_id = c("run-1", "run-2"),
    start_tp = c(1L, 21L),
    end_tp = c(20L, 40L)
  )

  out <- compute_sliding_connectivity(
    x = x,
    y = y,
    run_boundaries = runs,
    window_size_tp = 10L,
    step_size_tp = 1L,
    min_points_per_window = 8L,
    regress_pair_mean_signal = FALSE
  )

  expect_true(all(out$window_start_tp <= 20L | out$window_start_tp >= 21L))
  expect_true(all(out$window_end_tp <= 20L | out$window_end_tp >= 21L))
  expect_true(!any(out$window_start_tp <= 20L & out$window_end_tp >= 21L))
})

test_that("windows below min valid points are NA", {
  x <- seq_len(30)
  y <- seq_len(30)

  x[5:10] <- NA_real_
  y[5:10] <- NA_real_

  runs <- data.frame(
    run_id = "run-1",
    start_tp = 1L,
    end_tp = 30L
  )

  out <- compute_sliding_connectivity(
    x = x,
    y = y,
    run_boundaries = runs,
    window_size_tp = 12L,
    step_size_tp = 1L,
    min_points_per_window = 10L,
    regress_pair_mean_signal = FALSE
  )

  expect_true(any(is.na(out$fisher_z)))
  expect_true(any(out$n_valid_tp < 10L))
})

test_that("nuisance regression does not force perfect anti-correlation", {
  set.seed(101)
  n <- 120
  nuisance <- as.numeric(scale(sin(seq(0, 6, length.out = n)) + rnorm(n, sd = 0.2)))
  x <- 0.6 * nuisance + rnorm(n, sd = 0.7)
  y <- 0.5 * nuisance + rnorm(n, sd = 0.7)

  runs <- data.frame(
    run_id = c("run-1", "run-2"),
    start_tp = c(1L, 61L),
    end_tp = c(60L, 120L)
  )

  out <- compute_sliding_connectivity(
    x = x,
    y = y,
    run_boundaries = runs,
    window_size_tp = 12L,
    step_size_tp = 1L,
    min_points_per_window = 8L,
    regress_pair_mean_signal = TRUE,
    nuisance_signal = nuisance
  )

  valid_r <- out$pearson_r[!is.na(out$pearson_r)]
  expect_gt(length(valid_r), 0)
  expect_false(all(valid_r < -0.995))
})
