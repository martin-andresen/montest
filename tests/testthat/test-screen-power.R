mk <- function(t, rows = seq_along(t)) data.table::data.table(sel_t = t, train_row_id = rows, id = seq_along(t))

mk_u <- function(infl, G = 100, df = 90) {
  ## influence vectors with exact se: se^2 = sum(infl^2) * G/df
  lapply(infl, function(v) list(G = G, df = df, se = sqrt(sum(v^2) * G / df), u = v))
}

test_that("power screen: one strong cell keeps one", {
  d <- mk(c(-4, 0.2, 0.5))
  expect_equal(nrow(select_groups_screen(d, screen = "power")), 1L)
})

test_that("power screen: no negative t keeps one (keep-one fallback)", {
  d <- mk(c(0.3, 1, 2))
  expect_equal(nrow(select_groups_screen(d, screen = "power")), 1L)
})

test_that("power screen: two equal strong independent cells keeps two", {
  d <- mk(c(-3, -3))
  expect_equal(nrow(select_groups_screen(d, screen = "power")), 2L)
})

test_that("power screen: highly correlated strong cells keep one", {
  set.seed(1)
  base <- rnorm(60)
  nm <- as.character(1:60)
  u <- mk_u(list(setNames(base + rnorm(60, sd = 0.05), nm), setNames(base + rnorm(60, sd = 0.05), nm)))
  d <- mk(c(-3, -3))
  expect_equal(nrow(select_groups_screen(d, screen = "power", u = u)), 1L)
  ## independent influence vectors: keeps both
  u2 <- mk_u(list(setNames(rnorm(60), nm), setNames(rnorm(60), nm)))
  expect_equal(nrow(select_groups_screen(d, screen = "power", u = u2)), 2L)
})

test_that("cell_corr_from_u reproduces correlations and checks se", {
  nm <- as.character(1:40)
  set.seed(2)
  v1 <- setNames(rnorm(40), nm); v2 <- setNames(0.5 * v1 + sqrt(0.75) * rnorm(40), nm)
  u <- mk_u(list(v1, v2))
  R <- cell_corr_from_u(u)
  expect_equal(R[1, 2], stats::cor(v1, v2, method = "pearson") * 1, tolerance = 0.05)
  expect_equal(diag(R), c(1, 1), tolerance = 1e-8)
  u[[1]]$se <- u[[1]]$se * 2
  expect_warning(r <- cell_corr_from_u(u), "do not match")
  expect_null(r)
})

test_that("power_union is increasing in signal and restores the RNG", {
  set.seed(5); s0 <- .Random.seed
  p1 <- power_union(c(2, 2), diag(2) * 2, 0.025)
  expect_identical(.Random.seed, s0)
  p2 <- power_union(c(3, 3), diag(2) * 2, 0.025)
  expect_gt(p2, p1)
  expect_equal(power_union(2, matrix(2), 0.05), 1 - stats::pnorm(stats::qnorm(0.95), 2, sqrt(2)))
})

test_that("montest runs with screen = 'power' and guards unsupported settings", {
  skip_on_cran()
  data.table::setDTthreads(2)
  set.seed(7)
  n <- 3000
  x <- rnorm(n); z <- sample(0:4, n, TRUE); u <- runif(n); bad <- x > 0
  p <- 0.2 + 0.1 * z; p[bad] <- c(0.3, 0.4, 0.5, 0.35, 0.45)[z[bad] + 1]
  d <- as.integer(u < p); y <- rnorm(n) + 0.5 * d + bad * d
  dat <- data.frame(y, x, d, z)
  out <- montest(y ~ x | d ~ z, data = dat, condition = "MW", progress = FALSE,
                 select = c("zmargin", "equation"), screen = "power",
                 Rparameters = list(num.threads = 2))
  expect_true(all(is.finite(out$minp)))
  expect_true(any(out$results$train))
  expect_true(any(!out$results$train))
  expect_error(montest(y ~ x | d ~ z, data = dat, condition = "MW", progress = FALSE,
                       screen = "power", stack = FALSE), "not available")
})
