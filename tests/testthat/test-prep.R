test_that("delta detection picks modal integer gap", {
  expect_equal(.detect_delta(c(1L, 2L, 3L, 5L)), 1)
  expect_error(.detect_delta(c(1.5, 2.5)), "integer")
})

test_that("aw weights are renormalized to mean 1 when sum = FALSE", {
  df <- data.frame(i=c(1,1,2,2), t=c(1L,2L,1L,2L), Ei=c(2,2,NA,NA),
                   y=c(1,2,3,4), wt=c(1,1,3,3))
  p <- .prep_data(df, "y","i","t","Ei", controls=character(),
                  fe=NULL, cluster="", aw="wt", sum=FALSE, shift=0)
  # wei = wt * N / sum(wt); N=4, sum(wt)=8 -> factor 0.5
  expect_equal(sort(p$df$wei), sort(c(0.5,0.5,1.5,1.5)))
  expect_equal(sum(p$df$wei), nrow(p$df))   # renormalized total = N
})

test_that("rows with NA in a required column are dropped; NA in Ei is kept", {
  df <- data.frame(i=c(1,1,2,2,3), t=c(1L,2L,1L,2L,1L),
                   Ei=c(2,2,NA,NA,NA), y=c(1,NA,3,4,5))   # row 2 has NA y; rows with NA Ei are never-treated
  p <- .prep_data(df, "y","i","t","Ei", controls=character(),
                  fe=NULL, cluster="", aw=NULL, sum=FALSE, shift=0)
  expect_equal(nrow(p$df), 4)                 # the NA-y row dropped, NA-Ei rows kept
  expect_true(all(!is.na(p$df$y)))
  expect_equal(sum(p$df$untreated == 1), 4)   # only (1,t=2) would be treated but it's the dropped row -> all remaining untreated
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
