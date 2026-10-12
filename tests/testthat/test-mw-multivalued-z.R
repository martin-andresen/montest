test_that("MW accepts a multivalued instrument and detects a defier contrast", {
  set.seed(2)
  n <- 8000
  d <- data.table::data.table(x = rnorm(n), z = sample(0:2, n, replace = TRUE), u = runif(n))
  ## no violation: D increasing in z for everyone
  d[, d := as.integer(u < 0.2 + 0.2 * z)]
  d[, y := rnorm(.N) + 0.5 * d]
  ok <- montest(y ~ x | d ~ z, data = d, condition = "MW", progress = FALSE)
  expect_gt(ok$minp[["p.holm"]], 0.01)

  ## defiers only in the 1 -> 2 contrast
  d[, d := as.integer(data.table::fifelse(z == 0, u < 0.3,
                      data.table::fifelse(z == 1, u < 0.5, u < data.table::fifelse(x > 0, 0.2, 0.7))))]
  d[, y := rnorm(.N) + 0.5 * d + (x > 0) * d]
  bad <- montest(y ~ x | d ~ z, data = d, condition = "MW", progress = FALSE)
  expect_lt(bad$minp[["p.holm"]], 1e-6)
  ## the violation is absent once the Z = 2 rows are dropped
  sub <- montest(y ~ x | d ~ z, data = d[z < 2], condition = "MW", progress = FALSE)
  expect_gt(sub$minp[["p.holm"]], 0.01)
})
