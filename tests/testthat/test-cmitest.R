test_that("cmitest statistic matches a hand computation for a single cell", {
  set.seed(1)
  n <- 1000
  m <- rnorm(n, -0.1)
  r <- cmitest(m, rep(1, n), reps = 50, seed = 1)
  t_hand <- sqrt(n) * mean(m) / sqrt((mean(m^2) - mean(m)^2) * (1 + 0.05))
  expect_equal(r$stat, pmax(-t_hand, 0)^2)
  expect_equal(r$nboxes, 1L)
})

test_that("cmitest is invariant to rescaling the moments", {
  set.seed(2)
  n <- 800
  x <- rnorm(n)
  m <- rnorm(n, -0.05)
  a <- cmitest(m, x, reps = 100, seed = 3)
  b <- cmitest(5 * m, x, reps = 100, seed = 3)
  expect_equal(a$stat, b$stat)
  expect_equal(a$pval, b$pval)
})

test_that("maxboxes and min_obs control the number of boxes", {
  set.seed(3)
  n <- 1500
  x <- rnorm(n)
  m <- rnorm(n)
  expect_lte(cmitest(m, x, maxboxes = 20, reps = 20)$nboxes, 20)
  expect_lt(cmitest(m, x, min_obs = 50, reps = 20)$nboxes, cmitest(m, x, min_obs = 15, reps = 20)$nboxes)
  expect_equal(cmitest(m, x, r1 = 3, reps = 20)$nboxes, sum(2 * 1:3))
})

test_that("cmitest handles discrete, mixed and several conditioning variables", {
  set.seed(4)
  n <- 1200
  d <- sample(0:3, n, TRUE)
  x <- rnorm(n)
  x2 <- rnorm(n)
  m <- cbind(rnorm(n), rnorm(n))
  expect_s3_class(cmitest(m, d, reps = 50), "cmitest")
  expect_s3_class(cmitest(m, cbind(x, d), reps = 50, stat = "CvM", S = "S3"), "cmitest")
  expect_s3_class(cmitest(m, cbind(x, x2, d), reps = 50), "cmitest")
  expect_warning(cmitest(m, cbind(x, x2, rnorm(n), rnorm(n)), reps = 20), "columns")
  expect_error(cmitest(m, x[-1]))
})

test_that("cmitest is reproducible with a seed and rejects a clear violation", {
  set.seed(5)
  n <- 2000
  x <- rnorm(n)
  m <- rnorm(n) - 0.8 * (x > 0.5)
  a <- cmitest(m, x, reps = 100, seed = 11)
  b <- cmitest(m, x, reps = 100, seed = 11)
  expect_equal(a$pval, b$pval)
  expect_true(a$reject)
  expect_true(all(a$worst$t < 0))
})

test_that("include_all adds the whole-sample box and detects a diffuse violation", {
  set.seed(6)
  n <- 2000
  x <- rnorm(n)
  m <- rnorm(n, -0.07)
  a <- cmitest(m, x, reps = 100, seed = 1)
  b <- cmitest(m, x, reps = 100, seed = 1, include_all = TRUE)
  expect_equal(b$nboxes, a$nboxes + 1)
  expect_equal(unname(b$boxes_per_scale["all"]), 1)
  expect_gte(b$stat, a$stat)
  one <- cmitest(m, rep(1, n), reps = 20, seed = 1, include_all = TRUE)
  expect_equal(one$nboxes, 2)
})
