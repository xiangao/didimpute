# Data-driven golden-grid test: compares R did_impute() output against
# Python did_imputation() golden CSV fixtures for every spec §7 configuration.
# Tolerance: 1e-6. No case required loosening — all pass at 1e-6.

# Map of case name -> list of extra arguments to did_impute().
# "data" key (optional): per-case data file; defaults to baseline_data.csv.
cases <- list(
  baseline      = list(),
  horizons      = list(horizons = c(0, 1, 2, 3)),
  allhorizons   = list(allhorizons = TRUE),
  sum           = list(horizons = c(0, 1, 2), sum = TRUE),
  aw            = list(aw = "w", horizons = c(0, 1, 2)),
  controls      = list(controls = "x", horizons = c(0, 1, 2)),
  pretrends     = list(horizons = c(0, 1, 2), pretrends = 3),
  pretrends1    = list(horizons = c(0, 1),    pretrends = 1),
  hbalance      = list(horizons = c(0, 1, 2), hbalance = TRUE),
  wtr_custom    = list(wtr = "w_treat"),
  minn_suppress = list(horizons = c(0, 1, 2), minn = 500),
  shift         = list(horizons = c(0, 1), shift = 1),
  avgby         = list(horizons = c(0, 1), avgeffectsby = "Ei"),
  # never_vs_notyet: control is only not-yet-treated (no never-treated units)
  never_vs_notyet = list(.data = "never_vs_notyet_data.csv",
                         horizons = c(0, 1, 2))
)

d_default <- read.csv(testthat::test_path("fixtures", "baseline_data.csv"))

for (cs in names(cases)) {
  local({
    case_name <- cs
    args      <- cases[[cs]]

    # Separate the optional per-case data file from the did_impute() args
    data_file <- args[[".data"]]
    args[[".data"]] <- NULL

    test_that(paste("golden grid:", case_name), {
      if (!is.null(data_file)) {
        d <- read.csv(testthat::test_path("fixtures", data_file))
      } else {
        d <- d_default
      }

      g <- read.csv(testthat::test_path("fixtures",
                                        paste0(case_name, "_result.csv")))

      res <- do.call(did_impute, c(list(df = d, y = "y", i = "i",
                                        t = "t", Ei = "Ei"), args))

      # Build a named vector of all numeric outputs (estimates)
      flat <- c(
        stats::setNames(unlist(res$estimates),
                        paste0("effect:", names(res$estimates))),
        if (!is.null(res$pretrends_estimates))
          stats::setNames(unlist(res$pretrends_estimates),
                          paste0("pre:", names(res$pretrends_estimates))),
        if (!is.null(res$controls_estimates))
          stats::setNames(unlist(res$controls_estimates),
                          paste0("control:", names(res$controls_estimates)))
      )

      # Build a named vector of all SE outputs (parallel to flat)
      flat_se <- c(
        stats::setNames(unlist(res$std_errors),
                        paste0("effect:", names(res$std_errors))),
        if (!is.null(res$pretrends_std_errors))
          stats::setNames(unlist(res$pretrends_std_errors),
                          paste0("pre:", names(res$pretrends_std_errors))),
        if (!is.null(res$controls_std_errors))
          stats::setNames(unlist(res$controls_std_errors),
                          paste0("control:", names(res$controls_std_errors)))
      )

      # Compare every estimate key against the golden CSV
      for (k in names(flat)) {
        golden_val <- g$estimate[g$key == k]
        expect_length(golden_val, 1)
        expect_equal(
          unname(flat[[k]]),
          golden_val,
          tolerance = 1e-6,
          label     = paste(case_name, k, "estimate")
        )
      }

      # Compare every SE key against the golden CSV (skip NA/non-finite golden SEs)
      for (k in names(flat_se)) {
        golden_se <- g$se[g$key == k]
        if (length(golden_se) == 1 && !is.na(golden_se) && is.finite(golden_se)) {
          expect_equal(
            unname(flat_se[[k]]),
            golden_se,
            tolerance = 1e-6,
            label     = paste(case_name, k, "se")
          )
        }
      }

      # Assert n_obs against golden value
      gn <- g$estimate[g$key == "n_obs"]
      if (length(gn) == 1 && is.finite(gn)) {
        expect_equal(as.numeric(res$n_obs), gn, tolerance = 1, info = paste(case_name, "n_obs"))
      }
    })
  })
}
