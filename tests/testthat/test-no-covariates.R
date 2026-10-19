make_nocov_data <- function(n = 2500, seed = 1) {
  set.seed(seed)
  d <- data.frame(Z = sample(0:3, n, TRUE), cl = sample(seq_len(n / 5), n, TRUE))
  d$D <- rbinom(n, 1, 0.2 + 0.15 * d$Z)
  d$Y <- rnorm(n) + d$D
  d
}

test_that("no covariates and no selection: global estimates, no sample split", {
  d <- make_nocov_data()
  out <- montest(Y ~ 1 | D ~ Z, data = d, condition = "simple", select = "none",
                 progress = FALSE)
  expect_false(out$options$local)
  expect_equal(nrow(out$results), 1L)
  expect_false(out$results$train)
  expect_true(is.finite(out$results$t))
  expect_true(all(is.finite(out$minp)))
  expect_error(montestplot(out), "without a subgroup search")
})

test_that("no covariates with selection over margins: sample is split", {
  d <- make_nocov_data()
  out <- montest(Y ~ 1 | D ~ Z, data = d, condition = "simple", progress = FALSE)
  expect_true(out$options$local)
  expect_true(any(out$results$train))
  expect_true(any(!out$results$train))
})

test_that("empty-covariate spellings and options work", {
  d <- make_nocov_data()
  ref <- montest(Y ~ 1 | D ~ Z, data = d, select = "none", progress = FALSE)
  zero <- montest(Y ~ 0 | D ~ Z, data = d, select = "none", progress = FALSE)
  expect_equal(zero$results$coef, ref$results$coef)
  for (dr in c(TRUE, FALSE)) {
    out <- montest(Y ~ 1 | D ~ Z, data = d, doubly.robust = dr, progress = FALSE)
    expect_true(nrow(out$results) >= 1L)
  }
  expect_no_error(montest(Y ~ 1 | D ~ Z, data = d, cluster = "cl", select = "none", progress = FALSE))
  expect_no_error(montest(Y ~ 1 | D ~ Z, data = d, parametric = TRUE, progress = FALSE))
  expect_warning(montest(Y ~ 1 | D ~ Z, data = d, shrink = 0.5, progress = FALSE), "shrink")
})

test_that("KR without covariates runs, split and unsplit", {
  set.seed(2)
  n <- 3000
  d <- data.frame(Z = rbinom(n, 1, 0.5))
  d$D <- rbinom(n, 1, 0.3 + 0.3 * d$Z)
  d$Y <- rnorm(n) + d$D
  expect_no_error(suppressMessages(montest(Y ~ 1 | D ~ Z, data = d, condition = "KR", progress = FALSE)))
  out <- suppressMessages(montest(Y ~ 1 | D ~ Z, data = d, condition = "KR", select = "none", progress = FALSE))
  expect_false(any(out$results$train))
})

test_that("MW and AHS always search, also without covariates (Y is a feature)", {
  set.seed(3)
  n <- 4000
  d <- data.frame(Z = rbinom(n, 1, 0.5))
  d$D <- rbinom(n, 1, 0.3 + 0.3 * d$Z)
  d$Y <- rnorm(n) + d$D
  for (cc in c("MW", "AHS", "simple_MW")) {
    cond <- if (cc == "simple_MW") c("simple", "MW") else cc
    out <- suppressMessages(montest(Y ~ 1 | D ~ Z, data = d, condition = cond, progress = FALSE))
    expect_true(out$options$local, info = cc)
    expect_true(any(out$results$train), info = cc)
  }
})

test_that("nothing to test with binary D and Z and no covariates still stops", {
  set.seed(4)
  n <- 500
  d <- data.frame(Z = rbinom(n, 1, 0.5))
  d$D <- rbinom(n, 1, 0.3 + 0.3 * d$Z)
  d$Y <- rnorm(n)
  expect_error(montest(Y ~ 1 | D ~ Z, data = d, condition = "simple", progress = FALSE), "Nothing to test")
})

test_that("the removed `local` argument errors informatively", {
  d <- make_nocov_data()
  expect_error(montest(Y ~ 1 | D ~ Z, data = d, local = FALSE, progress = FALSE), "removed")
})
