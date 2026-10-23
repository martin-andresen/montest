test_that("priority_te divides by the variance and, with w, by the weight", {
  d <- data.table::data.table(
    sample = rep(1:2, each = 6),
    pred = c(-3:2, 4:-1) / 2, pred_o = c(2:-3, -1:4) / 3,
    v = rep(c(1, 2, 3), 4), vo = rep(c(2, 3, 1), 4),
    w = rep(c(0.5, 1, 2), 4)
  )
  base <- data.table::copy(d)
  priority_te(base, "pred", "v", "pred_o", "vo")
  expect_equal(base$pred_p, d$pred / d$v)
  expect_equal(base$pred_o_p, d$pred_o / d$vo)

  ## a constant weight leaves the key unchanged
  d[, one := 1]
  cst <- data.table::copy(d)
  priority_te(cst, "pred", "v", "pred_o", "vo", w = "one")
  expect_equal(cst$pred_p, base$pred_p)

  ## a varying weight enters the key as tau / (w v), also for the other half's key
  wt <- data.table::copy(d)
  priority_te(wt, "pred", "v", "pred_o", "vo", w = "w")
  expect_equal(wt$pred_p, d$pred / (d$v * d$w))
  expect_equal(wt$pred_o_p, d$pred_o / (d$vo * d$w))
})

test_that("MW honours target, for binary and multivalued treatments", {
  set.seed(3)
  n <- 6000
  d <- data.table::data.table(x = rnorm(n), z = sample(0:2, n, replace = TRUE), u = runif(n))
  ## valid IV
  d[, d := as.integer(u < 0.2 + 0.2 * z)]
  d[, y := rnorm(.N) + 0.5 * d]
  ## defiers in the 1 -> 2 contrast when x > 0
  v <- data.table::copy(d)
  v[, d := as.integer(data.table::fifelse(z == 0, u < 0.3,
                      data.table::fifelse(z == 1, u < 0.5, u < data.table::fifelse(x > 0, 0.2, 0.7))))]
  v[, y := rnorm(.N) + 0.5 * d + (x > 0) * d]
  for (tg in c("all", "overlap")) {
    ok <- montest(y ~ x | d ~ z, data = d, condition = "MW", target = tg, progress = FALSE)
    expect_gt(ok$minp[["p.holm"]], 0.01)
    bad <- montest(y ~ x | d ~ z, data = v, condition = "MW", target = tg, progress = FALSE)
    expect_lt(bad$minp[["p.holm"]], 1e-4)
  }

  ## three-valued treatment: interior sets use the same IPW contrast under both targets
  m <- data.table::data.table(x = rnorm(n), z = rbinom(n, 1, 0.5), u = runif(n))
  m[, d := (u < 0.7 - 0.2 * z) + (u < 0.35 - 0.1 * z)]
  m[, y := rnorm(.N) + 0.3 * d]
  for (tg in c("all", "overlap")) {
    r <- montest(y ~ x | d ~ z, data = m, condition = "MW", target = tg, progress = FALSE)
    expect_true(is.finite(r$minp[["p.holm"]]))
  }
})

test_that("MW cells default to minarm = minarm_clusters = 20, other conditions keep 10 and 5", {
  set.seed(5)
  n <- 1500
  x <- rnorm(n); v <- rnorm(n)
  z <- rbinom(n, 1, pmin(pmax(plogis(3.5 * x), 0.02), 0.98))
  dd <- as.integer(v + z > -0.5) + as.integer(v + z > 0.7)
  y <- as.integer(0.5 * dd + 0.5 * x + 0.5 * v + rlogis(n) > 0.5)
  d <- data.frame(Y = y, D = dd, Z = z, x = x)
  GRF <- list(num.trees = 100L, num.threads = 1)
  run <- function(...) suppressWarnings(montest(Y ~ x | D ~ Z, data = d, condition = "MW", progress = FALSE, pool = "none",
                                                Zparameters = GRF, Qparameters = GRF, Cparameters = GRF, ...))$results
  def <- run()
  m20 <- run(minarm = 20L, minarm_clusters = 20L)
  m10 <- run(minarm = 10L, minarm_clusters = 5L)
  expect_equal(def$p.raw, m20$p.raw)
  expect_equal(def$coef, m20$coef)
  expect_false(isTRUE(all.equal(def$coef, m10$coef)))
})
