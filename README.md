# syncfmri

An R package for sliding window fMRI functional connectivity analysis workflows.

## Installation

```r
# install.packages("pak")
pak::pak("jbartolotti/syncfmri")
```

## Development setup

```r
# install.packages(c("devtools", "roxygen2", "testthat", "usethis"))
devtools::document()
devtools::test()
devtools::check()
```

Open the project in RStudio using syncfmri.Rproj.

## Current functionality (v0)

- Discover ROI timecourse TSV files in a BIDS derivative.
- Compute run-aware sliding-window Pearson connectivity and Fisher-z values.
- Enforce run boundaries so windows do not cross run transitions.
- Drop windows below minimum valid paired timepoints.
- Preserve blank/censored rows from ROI TSV input as NaN samples.
- Optionally regress each ROI on a file-level nuisance signal (mean across numeric ROI columns) and use residuals.
- Write subject/group TSV + JSON sidecars + RDS outputs.
- Create a group heatmap (time x subject, color = Fisher-z) using viridis.
- Detect local maxima peaks using prominence and minimum distance rules.
- Write subject/group peak tables and peak-density tables.
- Create subject-level multi-panel figures (ROI activity, connectivity, peaks, events).
- Create group peak-density plots for all subjects and by subgroup.

## Example run

```r
library(syncfmri)

cfg <- syncfmri_default_config()
cfg$timecourse_derivative <- "gimmefMRI"
cfg$roi_x <- "LR_postJHipp_200"
cfg$roi_y <- "LR_mPFC_200"
cfg$run_boundaries <- data.frame(
	run_id = c("run-1", "run-2", "run-3"),
	start_tp = c(1L, 360L, 688L),
	end_tp = c(359L, 687L, 1012L)
)
cfg$tr_seconds <- 2
cfg$window_size_tp <- 12L
cfg$min_points_per_window <- 8L
cfg$regress_pair_mean_signal <- TRUE
cfg$peak_prominence_quantile <- 0.75
cfg$peak_min_distance_seconds <- NA_real_  # defaults to one window length in seconds
cfg$peak_positive_only <- TRUE
cfg$subgroup_column <- "group"
cfg$subgroup_levels <- c("YA", "OA")
cfg$peak_density_bin_seconds <- 4

result <- syncfmri_run_pipeline(
	bids_root = "P:/IRB_STUDY00149390_A015/MR_Data/Connectivity/natfMRI/A015_BIDS",
	config = cfg
)

# Regenerate only the group figure from previously saved windows:
syncfmri_run_pipeline(
	bids_root = "P:/IRB_STUDY00149390_A015/MR_Data/Connectivity/natfMRI/A015_BIDS",
	config = cfg,
	plot_only = TRUE
)
```

Study-specific wrapper scripts should be maintained outside this package repository.
