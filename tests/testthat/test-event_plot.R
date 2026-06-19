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
