test_that("event_plot builds a ggplot with the right x-positions", {
  skip_if_not_installed("ggplot2")
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2), pretrends = 2)
  p <- event_plot(res)
  expect_s3_class(p, "ggplot")
  xs <- sort(unique(p$data$x))
  expect_equal(xs, c(-2, -1, 0, 1, 2))      # pre at -k, effects at +k
})

test_that("event_plot with plot_type='rarea' still returns a ggplot", {
  skip_if_not_installed("ggplot2")
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2), pretrends = 2)
  p <- event_plot(res, plot_type = "rarea")
  expect_s3_class(p, "ggplot")
})

test_that("event_plot with together=TRUE returns a ggplot with correct x-positions", {
  skip_if_not_installed("ggplot2")
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2), pretrends = 2)
  p <- event_plot(res, together = TRUE)
  expect_s3_class(p, "ggplot")
  xs <- sort(unique(p$data$x))
  expect_equal(xs, c(-2, -1, 0, 1, 2))
})

test_that("event_plot with effects only (no pretrends) does not error", {
  skip_if_not_installed("ggplot2")
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2))
  p <- event_plot(res)
  expect_s3_class(p, "ggplot")
  xs <- sort(unique(p$data$x))
  expect_equal(xs, c(0, 1, 2))
})

test_that("event_plot with pretrends only (no horizons) does not error", {
  skip_if_not_installed("ggplot2")
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", pretrends = 2)
  p <- event_plot(res)
  expect_s3_class(p, "ggplot")
  xs <- sort(unique(p$data$x))
  expect_equal(xs, c(-2, -1))
})

test_that("event_plot via explicit dict args works and CIs use qnorm(1-alpha/2)*se", {
  skip_if_not_installed("ggplot2")
  pre_est  <- list(pre1 = 0.1, pre2 = -0.05)
  pre_se   <- list(pre1 = 0.02, pre2 = 0.03)
  eff_est  <- list(tau0 = 0.5, tau1 = 0.6)
  eff_se   <- list(tau0 = 0.05, tau1 = 0.07)
  p <- event_plot(pretrends = pre_est, pretrends_std = pre_se,
                  effects = eff_est, effects_std = eff_se)
  expect_s3_class(p, "ggplot")
  xs <- sort(unique(p$data$x))
  expect_equal(xs, c(-2, -1, 0, 1))
  # CI half-width check: err = qnorm(0.975) * se ≈ 1.96 * se
  cv <- qnorm(0.975)
  row0 <- p$data[p$data$x == 0, ]
  expect_equal(row0$err, cv * 0.05, tolerance = 1e-10)
})

test_that("together mode zero-fills the side missing SEs (Python parity)", {
  skip_if_not_installed("ggplot2")
  # effects have SEs, pretrends do not
  pre_est  <- list(pre1 = 0.1, pre2 = -0.2)
  eff_est  <- list(tau0 = 1.0, tau1 = 2.0)
  eff_se   <- list(tau0 = 0.3, tau1 = 0.4)
  p <- event_plot(pretrends = pre_est, pretrends_std = NULL,
                  effects = eff_est, effects_std = eff_se,
                  together = TRUE)
  expect_s3_class(p, "ggplot")
  # no NA in err: the pretrends side must be zero-filled, not NA
  expect_false(any(is.na(p$data$err)))
  # pretrends-side err (x < 0) should be 0, effects-side err (x >= 0) should be non-zero
  expect_true(all(p$data$err[p$data$x < 0] == 0))
  cv <- qnorm(0.975)
  expect_equal(p$data$err[p$data$x == 0], cv * 0.3, tolerance = 1e-10)
  expect_equal(p$data$err[p$data$x == 1], cv * 0.4, tolerance = 1e-10)
})
