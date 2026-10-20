test_that("montest runs with default pool/select", {
  set.seed(1)
  n <- 2000
  data <- data.frame(Z = rbinom(n, 1, 0.5), X1 = rnorm(n), X2 = rnorm(n))
  data$D <- rbinom(n, 1, 0.3 + 0.4 * data$Z)
  data$Y <- rnorm(n) + data$D
  out <- montest(Y ~ X1 + X2 | D ~ Z, data = data,
                 condition = "simple", progress = FALSE)
  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("a margin cannot be in both pool and select, and fits carry no adaptive outputs", {
  set.seed(1)
  n <- 2000
  data <- data.frame(Z = sample(0:2, n, TRUE), X1 = rnorm(n), X2 = rnorm(n))
  data$D <- rbinom(n, 1, 0.3 + 0.15 * data$Z)
  data$Y <- rnorm(n) + data$D
  expect_error(montest(Y ~ X1 + X2 | D ~ Z, data = data, condition = "simple",
                       pool = "zmargin", select = "zmargin", progress = FALSE),
               "must not share margins")
  out <- montest(Y ~ X1 + X2 | D ~ Z, data = data, condition = "simple", progress = FALSE)
  expect_false(any(grepl("^adaptive", names(out))))
})

test_that("sample cannot be selected over", {
  set.seed(1)
  n <- 1000
  data <- data.frame(Z = rbinom(n, 1, 0.5), X1 = rnorm(n))
  data$D <- rbinom(n, 1, 0.3 + 0.4 * data$Z)
  data$Y <- rnorm(n) + data$D
  expect_error(montest(Y ~ X1 | D ~ Z, data = data, condition = "simple", select = "sample",
                       progress = FALSE), "cannot be in `select`")
  expect_error(montest(Y ~ X1 | D ~ Z, data = data, condition = "simple", select = c("condition", "sample"),
                       progress = FALSE), "cannot be in `select`")
})
