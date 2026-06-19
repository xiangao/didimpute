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
})

test_that("summary with pretrends includes pretrend rows", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2),
                    pretrends = 3)
  s <- summary(res)
  expect_true(any(startsWith(s$term, "pre")))
})
