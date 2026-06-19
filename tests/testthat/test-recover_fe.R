test_that("single FE recovers weighted group means", {
  df <- data.frame(g = c("a","a","b"), combined = c(1, 3, 5), w = c(1, 1, 2))
  fe <- recover_fe(df, fe_cols = "g", combined_col = "combined", weight_col = "w")
  expect_equal(unname(fe$g["a"]), 2)   # mean(1,3)
  expect_equal(unname(fe$g["b"]), 5)
})

test_that("two FEs converge to an additive decomposition", {
  # combined = u[i] + v[t] exactly; recovery should reconstruct up to a constant shift
  df <- data.frame(
    i = c("1","1","2","2"), t = c("1","2","1","2"),
    combined = c(10+1, 10+2, 20+1, 20+2), w = 1
  )
  fe <- recover_fe(df, c("i","t"), "combined", "w")
  recon <- fe$i[df$i] + fe$t[df$t]
  expect_equal(unname(recon), df$combined, tolerance = 1e-6)
})
