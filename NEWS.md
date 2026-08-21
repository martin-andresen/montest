1# montest (development version)

* Initial CRAN submission.
* Fixed: with fixed effects (`has_FE`) and a binary instrument, `Z.hat` is
  estimated as an FE mean plus a fitted FE-residual, an unconstrained
  recombination that can land outside `[0,1]` -- previously this could hit a
  hard error from the closed-form `e(X)*(1-e(X))` propensity check (used
  whenever `doubly.robust = TRUE`, or `doubly.robust = FALSE` with
  `target = "all"`). Binary instruments estimated with fixed effects now
  always use the continuous/fitted-variance score representation instead
  (same one already used for multivalued instruments and `linearZ = TRUE`),
  where `Z.hat` enters the score only additively and a separately fitted
  `Var(Z|X,FE)` is used in place of the closed form -- targeting the same
  estimand, just without requiring a bounded propensity. That fitted
  variance nuisance is itself now floored (with a warning) rather than
  hard-erroring on small negative excursions near a zero-variance boundary
  (e.g. small/near-degenerate fixed-effect cells), since these are expected
  estimation noise from an unconstrained nonparametric fit of a
  nonnegative-by-construction target, not evidence of misspecification.
* Added: `drop_singletons` (default `TRUE`) and `drop_novar_Z` (default
  `TRUE`), two new `montest()` arguments that drop observations up front,
  before any nuisance estimation, when fixed effects are present.
  `drop_singletons` removes true fixed-effect singletons (an FE level with
  only one observation, following reghdfe's own default and implemented via
  fixest's `fixef.rm = "singletons"`) -- these affect degrees of freedom for
  every nuisance, not just `Z`'s. `drop_novar_Z` removes observations in FE
  cells with no within-cell variation in the instrument `Z` -- these are
  exactly zero-information for identifying the coefficient on `Z` in the
  classical FWL sense, and were the primary source of the "conditional
  variance v(X) for continuous-Z rows" clip warnings noted above. Deliberately
  asymmetric: an analogous lack of variation in `D` is *not* dropped, since
  such a cell still contributes a genuine zero-covariance data point to the
  classical estimator rather than being uninformative.
* Changed: `doubly.robust`'s default is now `NULL`, resolved from
  `parametric` instead of being unconditionally `TRUE`: `TRUE` whenever
  `parametric = FALSE` (the overall default), `FALSE` whenever
  `parametric = TRUE`. AIPW's orthogonality protection matters most when
  nuisances are flexibly/ML-estimated; under `parametric = TRUE` the
  nuisance functional form is trusted by construction, so AIPW buys
  comparatively little while still adding a second nuisance function's
  estimation noise on top of the singly-robust/FWL score. Passing
  `doubly.robust` explicitly still overrides this regardless of
  `parametric`.
* Added: `fml.varZ`, a new `montest()` argument giving a one-sided formula
  for the covariates used by the conditional-variance nuisance
  `v(X)=Var(Z|X,FE)` alone, separately from `fml.Z`'s covariates for `Z`'s
  own conditional mean. Defaults to `fml.Z` (which itself defaults to
  `fml`'s main X part), so leaving it unset reproduces the previous
  behavior exactly.
* Renamed: the `fe_rank` column in `forest_test()`'s and `CART_test()`'s
  output (`results`/`train`/`global`) is now `dof_rank`. It was never
  literally "the rank of the FE" -- it's the degrees-of-freedom correction
  from already-estimated parameters generally: under `parametric = TRUE`
  it includes the rank of both `X` and the FE, under the semiparametric
  default it includes only the FE rank, and it excludes the rank of an FE
  that clustering happens to be on. `fe_rank_adj`/`fe_rank_conservative`
  (the arguments controlling whether/how this is computed) are unchanged.
* Added: `recenter_test` (default `TRUE`), a new `montest()` argument.
  `normalize.Z` mean-shifts `Z.hat` only once, at whatever `sample`/
  `margins` level it operates on -- not separately within an adaptively-
  found leaf/subgroup, or within a "global" cell computed at a finer
  grouping than `normalize.Z` used. When `TRUE`, the FINAL (test-side and
  global) score average for AIPW (`doubly.robust=TRUE`) and the singly-
  robust score under `target="all"` is locally re-centered instead -- the
  propensity/treatment-residual term is re-demeaned using only that
  specific cell's own rows, the same finite-sample-bias correction
  `doubly.robust=FALSE & target="overlap"` already gets via its own (more
  thorough) centered/with-intercept construction. Implemented via a new
  `crv1_mean(recenter_propensity=TRUE)` mode that algebraically substitutes
  a locally-shifted `e' = e + mean(Z-e)` into the score formula, keeping
  the plug-in `tau` baseline (unlike `center=TRUE`, which drops it
  entirely) -- see `crv1_mean()`'s own docs for the derivation and a
  documented limitation: this does not re-estimate `tau` itself for the
  cell, so any local bias in the causal forest's own plug-in CATE average
  over that specific cell's rows (e.g. from adaptive selection correlating
  which rows land there with their predicted effect) is not corrected.
  Never applied to the training-side adaptive search, only to the final
  evaluation, and never to cells already using the `target="overlap" &
  doubly.robust=FALSE` centered path.
* Renamed and changed: `normalize.Z` is now `stabilize.scores`. For the
  continuous representation, behavior is unchanged (`Z.hat` is still
  mean-shifted per sample/margin group). For the binary/closed-form
  representation, the previous additive mean-shift of `Z.hat` is replaced by
  a Hajek/self-normalized-IPW stabilization: the score's actual IPW weight
  `Z/e(X)-(1-Z)/(1-e(X))` is rescaled per row so that treated rows'
  `mean(1/e(X))` and control rows' `mean(1/(1-e(X)))` each equal exactly 1
  within their own sample/margin group, folded into `make_scores_vec()`'s
  `v = e(1-e)` via a new `Z.stab` argument. `e(X)` itself is only clipped for
  these rows, no longer shifted -- the additive mean-shift was a less
  targeted fix than stabilizing the weight the score actually divides by,
  and is superseded rather than kept alongside it. Not applied under
  `doubly.robust=FALSE & target="overlap"`, same exclusion as the existing
  propensity clip. `recenter_test`'s local re-centering is not yet updated to
  match (still an additive correction) -- a known follow-up.
