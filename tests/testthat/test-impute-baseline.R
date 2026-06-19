test_that("baseline ATE matches the Python golden grid", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  g <- read.csv(test_path("fixtures", "baseline_result.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei")
  gold_ate <- g$estimate[g$key == "effect:tau_ate"]
  expect_equal(unname(res$estimates[["tau_ate"]]), gold_ate, tolerance = 1e-6)
})
