library(didimpute)

golden_check <- function(case, args) {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  g <- read.csv(test_path("fixtures", paste0(case, "_result.csv")))
  res <- do.call(did_impute, c(list(df = d, y = "y", i = "i", t = "t", Ei = "Ei"), args))
  for (k in g$key[startsWith(g$key, "effect:")]) {
    nm <- sub("^effect:", "", k)
    expect_equal(unname(res$estimates[[nm]]),
                 g$estimate[g$key == k], tolerance = 1e-6, info = paste(case, nm))
  }
  gn <- g$estimate[g$key == "n_obs"]
  if (length(gn) == 1 && is.finite(gn)) {
    expect_equal(as.numeric(res$n_obs), gn, tolerance = 1, info = paste(case, "n_obs"))
  }
}

test_that("horizons match python",    golden_check("horizons",    list(horizons = c(0, 1, 2, 3))))
test_that("allhorizons match python", golden_check("allhorizons", list(allhorizons = TRUE)))
test_that("sum match python",         golden_check("sum",         list(horizons = c(0, 1, 2), sum = TRUE)))
test_that("aw match python",          golden_check("aw",          list(aw = "w", horizons = c(0, 1, 2))))

test_that("n_obs is finite and positive for horizons case", {
  d <- read.csv(test_path("fixtures", "baseline_data.csv"))
  res <- did_impute(d, y = "y", i = "i", t = "t", Ei = "Ei", horizons = c(0, 1, 2, 3))
  expect_true(is.finite(res$n_obs) && res$n_obs > 0)
})
