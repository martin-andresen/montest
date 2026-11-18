test_that("recenter_propensity leaves rows without FWL/AIPW components (MW) on their own score", {
  set.seed(1)
  n <- 60
  score <- rnorm(n)
  Vr <- rnorm(n); Ur <- rnorm(n); tau <- rnorm(n); v <- runif(n, 0.2, 0.3)
  cl <- seq_len(n)
  mw <- rep(c(TRUE, FALSE), length.out = n)
  Vm <- replace(Vr, mw, NA); Um <- replace(Ur, mw, NA); vm <- replace(v, mw, NA)
  mixed <- crv1_mean(score, rep(1, n), cl, rank_adj = 0, recenter_propensity = "hajek",
                     resid_treat = Vm, resid_outcome = Um, tau = tau, v = vm)
  expect_true(is.finite(mixed$coef))
  expect_equal(mixed$N, n)
  ## manual: Hajek-normalized score for non-MW rows (normalizer over those rows), own score for MW rows
  hs <- hajek_score(tau[!mw], Vr[!mw], Ur[!mw], v[!mw], NULL, binary = FALSE, infl = TRUE)
  s_exp <- score
  s_exp[!mw] <- hs$score
  expect_equal(mixed$coef, mean(s_exp))
  ## all-MW group: plain mean and no extra parameter
  allmw <- crv1_mean(score, rep(1, n), cl, rank_adj = 0, recenter_propensity = "hajek",
                     resid_treat = rep(NA_real_, n), resid_outcome = rep(NA_real_, n), tau = tau, v = rep(NA_real_, n))
  plain <- crv1_mean(score, rep(1, n), cl, rank_adj = 0)
  expect_equal(allmw$coef, plain$coef)
  expect_equal(allmw$df, plain$df)
  ## no MW rows: plain Hajek score
  full <- crv1_mean(score, rep(1, n), cl, rank_adj = 0, recenter_propensity = "hajek",
                    resid_treat = Vr, resid_outcome = Ur, tau = tau, v = v)
  expect_equal(full$coef, mean(hajek_score(tau, Vr, Ur, v, NULL, binary = FALSE, infl = TRUE)$score))
})
