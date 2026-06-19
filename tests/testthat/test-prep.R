test_that("delta detection picks modal integer gap", {
  expect_equal(.detect_delta(c(1L, 2L, 3L, 5L)), 1)
  expect_error(.detect_delta(c(1.5, 2.5)), "integer")
})

test_that("untreated mask and Rel_time follow the BJS definition", {
  df <- data.frame(i = c(1, 1, 2, 2), t = c(1L, 2L, 1L, 2L),
                   Ei = c(2, 2, NA, NA), y = c(1, 2, 3, 4))
  p <- .prep_data(df, "y", "i", "t", "Ei", controls = character(),
                  fe = NULL, cluster = "", aw = NULL, sum = FALSE, shift = 0)
  d <- p$df[order(p$df$i, p$df$t), ]
  expect_equal(d$untreated, c(1, 0, 1, 1))       # (i=1, t=2) is treated
  expect_equal(p$fe, c("i", "t"))                 # default FE
  expect_equal(p$cluster, "i")                    # default cluster
  expect_equal(d$Rel_time[d$i == 1 & d$t == 2], 0)  # (t - Ei + shift) / delta = 0
})
