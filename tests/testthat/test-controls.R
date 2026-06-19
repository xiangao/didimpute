test_that("controls: estimates and SEs match python", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  g <- read.csv(test_path("fixtures", "controls_result.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei",
                    controls = "x", horizons = c(0, 1, 2))
  cc <- g[g$key == "control:x", ]
  expect_equal(unname(res$controls_estimates[["x"]]), cc$estimate, tolerance = 1e-6)
  expect_equal(unname(res$controls_std_errors[["x"]]), cc$se, tolerance = 1e-6)
  for (k in g$key[startsWith(g$key, "effect:")]) {
    nm <- sub("^effect:", "", k)
    expect_equal(unname(res$estimates[[nm]]), g$estimate[g$key == k], tolerance = 1e-6)
    expect_equal(unname(res$std_errors[[nm]]), g$se[g$key == k], tolerance = 1e-6)
  }
})
