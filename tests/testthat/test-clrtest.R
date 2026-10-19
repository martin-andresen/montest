test_that("biweight kernel is a density with support on [-1, 1]", {
  u <- seq(-1.5, 1.5, length.out = 30001)
  expect_equal(sum(montest:::.clr_kl(u)) * diff(u)[1], 1, tolerance = 1e-3)
  expect_equal(montest:::.clr_kl(c(-1.2, 1.01)), c(0, 0))
})

test_that("local linear estimate recovers a linear function exactly", {
  set.seed(1)
  x <- runif(500)
  y <- 2 + 3 * x
  z <- c(0.3, 0.5, 0.7)
  D <- outer(z, x, function(a, b) b - a)
  ghat <- montest:::.clr_llre_from(montest:::.clr_kl(D / 0.2), D, y)
  expect_equal(ghat, 2 + 3 * z, tolerance = 1e-8)
})

test_that("bandwidth follows the hand formula", {
  set.seed(2)
  n <- 400
  x <- rnorm(n)
  y <- sin(x) + rnorm(n, sd = 0.5)
  xx <- (x - mean(x)) / sd(x)
  X <- cbind(1, xx, xx^2, xx^3, xx^4)
  b <- solve(crossprod(X), crossprod(X, y))
  s2 <- mean((y - X %*% b)^2)
  d2 <- 2 * b[3] + 6 * b[4] * xx + 12 * b[5] * xx^2
  q <- quantile(xx, c(0.1, 0.9), names = FALSE)
  w <- as.numeric(xx >= q[1] & xx <= q[2])
  h_hand <- 2.036 * sd(x) * (s2 / sum(d2^2 * w))^(1 / 5)
  expect_equal(montest:::.clr_bandwidth(y, x, FALSE), h_hand)
  expect_equal(montest:::.clr_bandwidth(y, x, TRUE), h_hand * n^(1 / 5) * n^(-1 / 3.5))
})

test_that("clrtest does not reject a valid null and rejects a clear violation", {
  set.seed(3)
  n <- 1500
  x <- rnorm(n)
  y0 <- -0.3 + rnorm(n)                                  # mean < 0 everywhere
  r0 <- clrtest(y0, x, rnd = 2000)
  expect_false(r0$reject)
  y1 <- -0.3 + 1.2 * (x > 0.5) + rnorm(n)                # positive mean on x > 0.5
  r1 <- clrtest(y1, x, rnd = 2000)
  expect_true(r1$reject)
  expect_lt(r1$pval, 0.05)
})

test_that("clrtest handles several inequalities, AIS and reproducibility", {
  set.seed(4)
  n <- 1000
  x <- rnorm(n)
  y <- cbind(-0.4 + rnorm(n), -0.4 + 1.5 * (x < -0.5) + rnorm(n))
  a <- clrtest(y, x, rnd = 1000, seed = 5, levels = c(0.9, 0.95))
  b <- clrtest(y, x, rnd = 1000, seed = 5, levels = c(0.9, 0.95))
  expect_equal(a$pval, b$pval)
  expect_true(all(a$n_selected <= a$n_grid))
  expect_equal(length(a$bounds), 2)
  expect_equal(unname(a$det), as.integer(a$bounds < 0))
  noais <- clrtest(y, x, rnd = 1000, seed = 5, ais = FALSE)
  expect_equal(noais$n_selected, noais$n_grid)
})

test_that("x can be a vector or a matrix with one column per inequality", {
  set.seed(6)
  n <- 600
  x <- rnorm(n)
  y <- cbind(rnorm(n, -0.5), rnorm(n, -0.5))
  v <- clrtest(y, x, rnd = 500, seed = 1)
  m <- clrtest(y, cbind(x, x), rnd = 500, seed = 1)
  expect_equal(v$stat, m$stat)
  expect_error(clrtest(y, x[-1]))
})
