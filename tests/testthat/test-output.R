test_that("print returns object invisibly and shows 'didimpute'", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1))
  expect_output(print(res), "didimpute")
  ret <- withVisible(print(res))
  expect_false(ret$visible)
})

test_that("summary returns data.frame with term/estimate/std.error", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1))
  s <- summary(res)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("term", "estimate", "std.error") %in% names(s)))
  expect_true(all(c("tau0", "tau1") %in% names(res$estimates)))
  s <- summary(res)
  expect_true(all(c("tau0", "tau1") %in% s$term))
})

test_that("V is a KxK covariance matrix for pretrends+controls run", {
  set.seed(99)
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  d$ctrl <- rnorm(nrow(d))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1),
                    pretrends = 2, controls = "ctrl", minn = 0)
  # K = 2 effects + 2 pretrends + 1 control = 5
  K <- length(res$estimates) + length(res$pretrends_estimates) + length(res$controls_estimates)
  expect_true(is.matrix(res$V))
  expect_equal(dim(res$V), c(K, K))

  # sqrt(diag(V)) must equal reported SEs in documented order
  reported_ses <- c(unlist(res$std_errors),
                    unlist(res$pretrends_std_errors),
                    unlist(res$controls_std_errors))
  expect_equal(unname(sqrt(diag(res$V))), unname(reported_ses), tolerance = 1e-8)

  # V must be symmetric
  expect_true(isSymmetric(res$V))

  # dimnames must match estimand names
  est_names <- c(names(res$estimates), names(res$pretrends_estimates),
                 names(res$controls_estimates))
  expect_equal(rownames(res$V), est_names)
  expect_equal(colnames(res$V), est_names)

  # Off-diagonal entries should generally be non-zero (cross-covariances exist)
  off_diag <- res$V[upper.tri(res$V)]
  expect_true(any(abs(off_diag) > 1e-12),
              info = "V should have non-zero off-diagonal entries (cross-covariances)")
})

test_that("V is NULL when nose=TRUE", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", nose = TRUE)
  expect_null(res$V)
})

test_that("summary with pretrends includes pretrend rows", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2),
                    pretrends = 3)
  s <- summary(res)
  expect_true(any(startsWith(s$term, "pre")))
})
