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
