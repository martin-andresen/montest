1# montest (development version)
* Changed: `montest()` no longer accepts `"sample"` in `select` (error; `select = "all"` now means every margin except `sample`). Selecting a sample direction ranks the two directions' candidates together, so whether a direction is carried forward depends on the other half's outcomes, which are the held-out data of the first direction. Leave `sample` in neither to test the directions separately, or in `pool`.
* Removed: the adaptive choice between pooling and selecting over a margin (the same margin in `pool` and `select`) and its outputs `$adaptive_dims`, `$adaptive_design`, `$adaptive_scores`, `$adaptive_pool`, `$adaptive_select`. It fitted every pool/select design and kept the winner on the training rows; simulations showed no clear gain over forced pooling or selection, and a choice made from both halves is not independent of the held-out tests. The code and a note are archived in `dev/adaptive/` and in the git tag `archive/adaptive-pool-select`. `forest_test_core()` is renamed `forest_test()`.
* Removed: the `studentize` argument of `montest()` (sorting on `pred / sqrt(pred_var)`); `priority` (default `TRUE`) replaces it. `studentize_te()` is gone.
* Changed: `montest()` drops `"sample"` from `pool` (with a message) whenever `select` contains a margin with at least two cells for the data and conditions (e.g. several conditions, or `select = "zmargin"` with a multivalued instrument). The cells selected by both sample directions were pooled across halves, but that choice uses both halves' outcomes, each of which is the held-out half of the other direction; the directions are now tested separately. Without effective selection, pooling the halves is unchanged.
* Added: `screen = "power"` in `montest()` (forest search, `stack = TRUE`): within each group of cells that differ only in the `select` margins, keep the top-j cells (by training t) that maximize the predicted probability that at least one is rejected by Holm on the held-out half. Held-out statistics are predicted from the training t-statistics (variance inflated by `1 + phi`, `options(montest.power_phi = 1)`) and the training-side correlation between cells, computed from cluster-level influence vectors; no prior. Opt-in; the default stays `"stepdown"`. Adds `mvtnorm` to Imports.
* Changed: defaults of `pool` and `select` in `montest()`. `pool` now defaults to `c("sample", "zmargin", "dval", "yval", "equation")` and `select` to `c("condition", "block")`. `block`, when present, is selected over by default. A margin the user specifies in one of the two is dropped from the default of the other; a margin in both `pool` and `select` is an error.
* Changed: `montest()` now supports no covariates. The `local` argument is removed (supplying it is an error) and is inferred: the sample split and subgroup search run whenever `fml`/`fml.C` has covariates or `condition` includes `"MW"`/`"AHS"` (which use Y as a feature, so they also run without covariates). With no covariates and only `"simple"`/`"KR"`, the sample is split only if `select` selects over margin cells (cells are selected on one half and reported on the other, with the causal forest on a constant feature); with no selection (e.g. `select = "none"`) there is no model selection and no split, and the global estimates within each margin cell are reported (the former `local = FALSE`). `shrink`, `priority` and `testtype = "CART"` are ignored without covariates. The old `Y ~ X` requirement ("The main X part of `fml` is empty") is gone; use `Y ~ 1 | D ~ Z`.
* Changed: `condition = "MW"` now accepts a multivalued instrument. As for the other conditions, `zmargin` stacks adjacent contrasts (Z = z_k vs z_{k+1}), each a binary-instrument MW problem, so the guard requiring a binary Z is removed. The MW paragraph of the help file is updated (multivalued D, `mw_sets`).
* Added: `minarm` (default 10) and `minarm_clusters` (default 5) in `montest()`: in the local (forest/CART) tests with a binary instrument, every candidate subgroup must contain at least `minarm` rows and `minarm_clusters` clusters with the instrument equal to 1 and with it equal to 0 (in both the training and estimation orderings). Subgroups without instrument variation identify nothing and gave degenerate standard errors and t-statistics of -50 or less. For CART the guard applies to the training leaves only. Margin cells that cannot satisfy `minsize`/`minarm` are skipped with a warning (previously an error); the error remains if no cell is feasible.
* Added: `memcheck` (default `TRUE`) in `montest()`: before the margin index is crossed with the data, the size of the crossed data is estimated and `montest()` stops with an error (suggesting `stack = FALSE` or fewer margins) if it would not fit in free memory. `memcheck = FALSE` skips the check.
* Changed: `seqtest()` with `"KRDY2"`/`"MWDY2"` now drops rows with missing values in D1, D2 or Y (with a message, for all requested conditions) instead of stopping.
* Added: `priority.var` (default `"model"`) in `montest()`. With `priority = TRUE`, the score variance `v(X) = E[(S - tau(X))^2 | X]` is now computed in closed form, without a forest, for rows with a binary instrument and an outcome `Q` taking two values in the margin cell (binarized-treatment cells of `"simple"`/`"AHS"`, and `"KR"`): `v = s1^2/e + s0^2/(1-e)` with arm variances `D^2 p_z (1 - p_z)` from `Z.hat`, `Q.hat` and the unshrunk predicted effect (plus a between-arm term when `doubly.robust = FALSE`). This removes the heavy-tailed `1/v(Z|X)^2` forest target and, with fixed effects, the full-sample fixed-effect means of the squared scores; `linearZ`, `linearD`, `"MW"` and continuous `Q` still use the forest. `priority.var = "scores"` restores the previous behaviour.
* Changed: the variance nuisances (`Var(Z|X)` behind `fml.varZ`, and the score variance behind `priority`) are now always estimated with regression forests, independently of `parametric`; they no longer add to the X-rank degrees-of-freedom penalty under `parametric = TRUE`. Fixed effects are still partialled out linearly before the forest.
* Added: `priority` (default `TRUE`, switched off automatically without covariates in the forest formula or with `testtype = "CART"`) in `montest()`: the subgroup search sorts on `pred / v(X)` (Neill, 2012, JRSS-B), with `v(X) = E[(S - tau(X))^2 | X]` estimated by a regression forest on the unshrunk predictions (own half out-of-bag, other half fitted on the search half only, so sample-part purity is preserved); `shrink` applies to `pred` before the division. New arguments `priority.floor` (floor of the variance as a fraction of its median, default 0.01), `fml.C` (covariates of the causal forest/search; default: X part of `fml`), `fml.varS` (covariates of the variance model; default `fml.C`) and `Sparameters`. `montestplot()` relabels the x axis.
* Fixed: `montest()` (and `seqtest()`) crashed or, worse, silently returned wrong numbers when a column of `data` had the same name as a variable used inside the package's own `data.table` calls (e.g. `y`, `weight`, `cluster`, `pred_buf`, `D1`, `sample`). Such columns are now renamed to internal aliases for the duration of the call (formulas and the `weight`/`cluster`/`block` arguments are rewritten, and names are mapped back in the output); the list of affected names is generated by `data-raw/reserved_names.R`.
* Added: with `stack = TRUE`, `progress = TRUE` now shows one titled progress bar per time-consuming step (nuisance fits, causal forests), titled like the `$time` entries; with `stack = FALSE` a single bar runs over the margin blocks.

* Added: `block` argument to `montest()` (for `condition = "KR"`): a column indexing independent problems stacked in one dataset. The supports, sets A, one-sided screen and `Q` are determined within each block, cells only use their own block's rows, and `"block"` is accepted in `pool`/`select`. Requires `cluster` (unit id). `seqtest()` uses it for the new condition `"KRDY2"`, the intersection of KR with D1 as treatment and (D2, Y) as outcome and KR with D2 as treatment and Y as outcome, so that cells can be pooled or selected across the two sets (`pool = "block"`/`select = "block"`). With `block`, `Z.hat` and the variance nuisance are fitted once on the first block and copied to the others when the blocks are row-aligned copies (so a unit's duplicate cannot be in-bag for its own out-of-bag prediction).

* Added: `stack` (default `TRUE`, the previous behaviour) and `progress`. With `stack = FALSE`, `montest()` no longer crosses the margin index (condition, dval, yval, equation) with the data but tests one block of the margin index at a time (still stacked across `zmargin`), shows a progress bar, appends all train and test estimates, and finally applies the `select`/`screen` rule jointly on the training side, dropping test rows that were not selected. `pool` is restricted to `"zmargin"`/`"sample"`; works with `testtype = "forest"` and `"CART"` and with `local = FALSE`. `seqtest()` passes it through `...`.
* Fixed: `testtype = 'CART'` stopped with 'margin column yval is NA' whenever several conditions were stacked (e.g. `condition = c('simple', 'KR')`), because `yval` is `NA` for the conditions without outcome sets. `NA` is now treated as an ordinary margin level.

* Added: `local` (default `TRUE`). With `local = FALSE`, the sample split and the local subgroup search (`forest_test()`/`CART_test()`) are skipped entirely and the estimates are computed globally within each margin cell (centered where the local test-side moment is), stored in `$results` and feeding `$minp`. `pool` still applies; `select`/`shrink`/`testtype` have no role; no causal forest is fit when `doubly.robust = FALSE`.
* Changed: `$global` now uses the same centered (with-intercept) moment as the local test side for `doubly.robust = FALSE` rows, so it can differ slightly from before.
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
* Added: the FINAL (test-side and global) score average for AIPW
  (`doubly.robust=TRUE`) and the singly-robust score under `target="all"` is
  now locally re-corrected within adaptively-found subgroups/cells, matching
  whichever correction `stabilize.scores` applied at the design/global level
  (see below) -- the same finite-sample-bias correction
  `doubly.robust=FALSE & target="overlap"` already gets via its own (more
  thorough) centered/with-intercept construction. Implemented via a new
  `crv1_mean(recenter_propensity=TRUE)` mode, with a `recenter_binary`
  switch selecting the local analogue: for the continuous representation, a
  locally-shifted `e' = e + mean(Z-e)` substituted into the score formula;
  for the binary/closed-form representation, a fresh Hajek/self-normalized-
  IPW rescaling of `v` recomputed from that cell's own rows. `montest()`
  sets `recenter_binary` automatically from whichever representation is
  actually in use -- see `stabilize.scores`'s docs and `crv1_mean()`'s own
  docs for the full derivation and a documented limitation: this does not
  re-estimate `tau` itself for the cell, so any local bias in the causal
  forest's own plug-in CATE average over that specific cell's rows (e.g.
  from adaptive selection correlating which rows land there with their
  predicted effect) is not corrected. Never applied to the training-side
  adaptive search, only to the final evaluation, and never to cells already
  using the `target="overlap" & doubly.robust=FALSE` centered path.
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
  propensity clip. There is no separate argument controlling the local
  test-side re-correction above -- it is driven entirely by
  `stabilize.scores`, matching whichever correction is in use.
