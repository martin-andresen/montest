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
