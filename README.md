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
- Optionally regress each ROI on pair-mean signal and use residuals.
- Write subject/group TSV + JSON sidecars + RDS outputs.
- Create a group heatmap (time x subject, color = Fisher-z) using viridis.

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

result <- syncfmri_run_pipeline(
	bids_root = "P:/IRB_STUDY00149390_A015/MR_Data/Connectivity/natfMRI/A015_BIDS",
	config = cfg
)
```

The minimal wrapper script is available at `inst/scripts/run_syncfmri_wrapper.R`.
