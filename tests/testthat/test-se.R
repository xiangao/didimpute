se_check <- function(case, args) {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  g <- read.csv(test_path("fixtures", paste0(case, "_result.csv")))
  res <- do.call(did_impute, c(list(df = d, y = "y", i = "i", t = "t", Ei = "Ei"), args))
  for (k in g$key[startsWith(g$key, "effect:")]) {
    nm <- sub("^effect:", "", k)
    expect_equal(unname(res$std_errors[[nm]]),
                 g$se[g$key == k], tolerance = 1e-6, info = paste(case, nm))
  }
}
test_that("baseline ATE SE matches python", se_check("baseline", list()))
test_that("horizon SEs match python",       se_check("horizons", list(horizons = c(0, 1, 2, 3))))
