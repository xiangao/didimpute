test_that("unimplemented options raise the same errors as Python", {
  d <- read.csv(testthat::test_path("fixtures", "baseline_data.csv"))

  # Time-interacted controls not yet implemented
  expect_error(
    did_impute(d, "y", "i", "t", "Ei", timecontrols = "x"),
    "Time and unit interacted controls are not added yet\\."
  )

  # Leaveoneout not yet implemented
  expect_error(
    did_impute(d, "y", "i", "t", "Ei", leaveoneout = TRUE),
    "Leaveoneout standard errors and options are not added yet\\."
  )

  # hetby not yet implemented
  expect_error(
    did_impute(d, "y", "i", "t", "Ei", hetby = "x"),
    "Hetby and project options are not added yet\\."
  )

  # wtr + horizons cannot be combined
  expect_error(
    did_impute(d, "y", "i", "t", "Ei", wtr = "w", horizons = c(0, 1)),
    "User provided weights and horizons options can not be combined\\."
  )
})
