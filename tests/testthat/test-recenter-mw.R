test_that("recenter_propensity leaves rows without FWL/AIPW components (MW) on their own score", {
  set.seed(1)
  n <- 60
  score <- rnorm(n)
  Vr <- rnorm(n); Ur <- rnorm(n); tau <- rnorm(n); v <- runif(n, 0.2, 0.3)
  cl <- seq_len(n)
  mw <- rep(c(TRUE, FALSE), length.out = n)
  Vm <- replace(Vr, mw, NA); Um <- replace(Ur, mw, NA); vm <- replace(v, mw, NA)
  mixed <- crv1_mean(score, rep(1, n), cl, rank_adj = 0, recenter_propensity = TRUE,
                     resid_treat = Vm, resid_outcome = Um, tau = tau, v = vm)
  expect_true(is.finite(mixed$coef))
  expect_equal(mixed$N, n)
  ## manual: recentered score for non-MW rows (Vbar over those rows), own score for MW rows
  Vbar <- mean(Vr[!mw])
  s_exp <- ifelse(mw, score, tau + (Vr - Vbar) / v * (Ur + tau * Vbar))
  expect_equal(mixed$coef, mean(s_exp))
  ## all-MW group: plain mean and no extra parameter
  allmw <- crv1_mean(score, rep(1, n), cl, rank_adj = 0, recenter_propensity = TRUE,
                     resid_treat = rep(NA_real_, n), resid_outcome = rep(NA_real_, n), tau = tau, v = rep(NA_real_, n))
  plain <- crv1_mean(score, rep(1, n), cl, rank_adj = 0)
  expect_equal(allmw$coef, plain$coef)
  expect_equal(allmw$df, plain$df)
  ## no MW rows: identical to the unchanged formula
  full <- crv1_mean(score, rep(1, n), cl, rank_adj = 0, recenter_propensity = TRUE,
                    resid_treat = Vr, resid_outcome = Ur, tau = tau, v = v)
  Vb <- mean(Vr)
  expect_equal(full$coef, mean(tau + (Vr - Vb) / v * (Ur + tau * Vb)))
})
