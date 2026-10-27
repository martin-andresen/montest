make_fsy_data <- function(n = 2000, seed = 1) {
  set.seed(seed)
  X <- rnorm(n)
  Z <- rbinom(n, 1, .5)
  D1 <- rbinom(n, 1, plogis(-.5 + 1.2 * Z))
  D2 <- rbinom(n, 1, plogis(-1 + .8 * D1 + .3 * Z))
  Y <- rnorm(n) + D1 + D2
  data.frame(Y, X, Z, D1, D2)
}

test_that("FSDY and FSSY run, with bins as the yval margin", {
  skip_on_cran()
  dat <- make_fsy_data()
  for (cn in c("FSDY", "FSSY")) {
    r <- seqtest(Y ~ X | D1 + D2 ~ Z, dat, condition = cn, progress = FALSE, select = "none",
                 pool = c("sample", "zmargin", "dval", "equation"))
    expect_s3_class(r, "seqtest")
    expect_true(!is.null(r$Wlookup))
    expect_gte(data.table::uniqueN(r$results$yval), 2L)
  }
})

test_that("FSDY/FSSY bin estimates equal the first stage of the summed cells; pooling gives FSD/FSS", {
  skip_on_cran()
  dat <- make_fsy_data(8000, seed = 2)
  dat$b <- cut(dat$Y, stats::quantile(dat$Y, 0:4 / 4), include.lowest = TRUE, labels = FALSE)
  fs <- function(q) mean(q[dat$Z == 1]) - mean(q[dat$Z == 0])
  pl <- c("sample", "zmargin", "dval", "equation")
  for (cn in c("FSDY", "FSSY")) {
    mom <- function(b) if (cn == "FSDY") (dat$D1 - dat$D2) * (dat$b == b) else (dat$D1 + dat$D2 - 1) * (dat$b == b)
    r <- seqtest(Y ~ 1 | D1 + D2 ~ Z, dat, condition = cn, progress = FALSE, select = "none", pool = pl)
    est <- r$results[order(yval)]$coef
    expect_equal(est, vapply(1:4, function(b) fs(mom(b)), numeric(1)), tolerance = 1e-4)
    rp <- seqtest(Y ~ 1 | D1 + D2 ~ Z, dat, condition = cn, progress = FALSE, select = "none")
    expect_equal(rp$results[train == FALSE]$coef, mean(est), tolerance = 1e-4)
  }
})

test_that("FSDY/FSSY need Y and binary treatments, and are in 'all'", {
  dat <- make_fsy_data(500)
  expect_error(seqtest(~ X | D1 + D2 ~ Z, dat, condition = "FSDY", progress = FALSE), "requires an outcome")
  dat$D2 <- dat$D2 + rbinom(nrow(dat), 1, .3)
  expect_error(seqtest(Y ~ X | D1 + D2 ~ Z, dat, condition = "FSSY", progress = FALSE), "binary")
})
