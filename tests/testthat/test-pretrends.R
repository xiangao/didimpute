test_that("pretrends estimates and SEs match python", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  g <- read.csv(test_path("fixtures", "pretrends_result.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2), pretrends = 3)
  for (k in g$key[startsWith(g$key, "pre:")]) {
    nm <- sub("^pre:", "", k)
    expect_equal(unname(res$pretrends_estimates[[nm]]), g$estimate[g$key == k], tolerance = 1e-6)
    expect_equal(unname(res$pretrends_std_errors[[nm]]), g$se[g$key == k], tolerance = 1e-6)
  }
})

test_that("nose=TRUE returns pretrends estimates but NULL SEs", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2),
                    pretrends = 3, nose = TRUE)
  expect_false(is.null(res$pretrends_estimates))
  expect_null(res$pretrends_std_errors)
  expect_null(res$std_errors)
})
