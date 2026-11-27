1# montest (development version)

* Changed: with fixed effects the causal forest and the subgroup search now split on the covariates as they are, not on their residuals from the fixed effects. New argument `C.features`: `"raw"` (default), `"residualized"` (the earlier behaviour) or `"both"`. The nuisances of Z, Q and the score variance still use residualized covariates. Results change for calls with fixed effects and search covariates.
* Fixed: a covariate that the fixed effects absorb (e.g. a state characteristic with state fixed effects) entered the nuisance forests, and the search, as its residual, which is numerical noise that the forests split on, so that subgroups described in that variable were meaningless. Such a covariate is now dropped from the nuisances and, with `C.features = "residualized"`, enters the search as it is; absorbed in some margin cells only, its residual there is set to zero. Under `parametric = TRUE`, a formula naming only absorbed covariates no longer stops with `RHS of assignment to existing column ... is zero length`.
* Fixed: with a few hundred covariate columns (typically with fixed effects) `montest()` stopped with `object '__xf_res_..._resid' not found` or `Missing columns: ...`: `data.table` ran out of spare column slots and the new columns were lost. Slots are now reserved according to the width of the covariate formulas.
* Changed: with several instruments `montest()` stops if `data` already has a column it creates itself when stacking (`unit__`, `inst__`, `Zstack__`, and with `Zcontrol = TRUE` `Zc1__`, ...). Such a column was overwritten, so that e.g. `cluster = "unit__"` would silently have been replaced by row numbers.

# montest 1.0.0

The version number marks releases; the GitHub commit ID identifies the exact code, see `?montest_version`.

* Added: `montest_version()` identifies the running code (`version`, git commit `sha`, `dirty` for uncommitted changes under `load_all()`, and `code_md5`, a checksum of the loaded functions). `montest()`, `seqtest()`, `montest_seeds()` and `seqtest_seeds()` store it as `$version`. The commit is `NA` after installing from a local folder.
* Added: `montest_seeds()` runs one `montest()` call over a set of seeds fixed in advance (default `1:51`) and combines the minimum p-values of each correction method as `p_agg = min(1, 2 * median(p))` over seeds (failed seeds count as p = 1), next to the share of seeds rejecting. With `file`, runs are saved after every seed and an interrupted call resumes; it stops if the formula, arguments, data, seeds, package code or relevant package versions have changed. Seeds run one at a time.
* Added: `seqtest_seeds()`, the same for `seqtest()`: `p_agg` combines its top-level `$minp` (all cells of all conditions as one family) and `p_agg_condition` the minimum p-values of each condition by itself (not corrected across conditions). A saved file is resumed only by the function that created it.
* Added: several instruments in `montest()`, `Y ~ X | D ~ Z1 + Z2` (condition `"simple"`, `"KR"` or `"MW"`, one per call). Each instrument is one block of a stacked problem (block label = instrument name) with the other instruments as covariates (`Zcontrol = TRUE`, default; `FALSE` tests each alone). By default each instrument is tested on its own and all tests are corrected as one family; `pool = "instrument"` pools and `select = "instrument"` selects across instruments (aliases of `"block"`). Instruments are binned separately (`Zcut` may be a named list) and may have different numbers of levels (margins `zmargin` are built per instrument; a two-valued instrument in a stack with a multivalued one shows up with `zmargin = 1`); D and Y are binned once. Z.hat, v(X) and the overlap checks are computed per instrument. New `$instruments` in the result. Previously only the first instrument was used, silently. Fixed effects are interacted with the instrument block (singletons, cells without instrument variation and the absorbed rank are determined per instrument). `"simple"` allows `linearD`/`linearZ`; its global test (`fml.C = ~1`) reports `$joint`, the joint first-stage Wald F with the covariance of the instruments' coefficients; with `Zcontrol = FALSE` the coefficients come from one regression on all instrument residuals (covariates and fixed effects partialled out only).
* Added: with several instruments `montest()` stops if an instrument is constant within every combination of the others (e.g. their interaction, or a complete set of dummies) while `Zcontrol = TRUE`; such a call used to run and report uninterpretable results for that instrument. Documented that several instruments test vector monotonicity (each instrument raises D at every value of the others and of X), not partial or Imbens-Angrist monotonicity, and that `Zcontrol = FALSE` is implied by the assumptions only for instruments that are independent given X.
* Added: one binning syntax, `Zcut`, `Dcut` and `Ycut` in `montest()` and `D1cut`, `D2cut`, `D2Tcut` and `Ycut` in `seqtest()`. Each takes `"none"`, `"equisized:K"`, `"equidistant:K"` (via `binarize_var`) or a vector of cut points (left-open, right-closed bins, coded 0, 1, ...); the default is `"equisized:4"`, as before. `"none"` is refused above 100 distinct values. `montest()` returns `$binning`, the spec used per variable and, for cut points, the intervals the codes stand for. Cutting D or Z does not relax the checks that use the true support (`"AHS"` still needs a binary D in the data).
* Deprecated: `Zsubsets`, `Dsubsets`, `Ysubsets`, `gridtypeZ`, `gridtypeD` and `gridtypeY`. They still work (`Dsubsets = K, gridtypeD = g` is `Dcut = "g:K"`) and warn once per session; passing one together with the new argument is an error.
* Added: `seqtest()` bins D1, D2 and Y by role, so they can be passed unbinned. D1 is always a treatment (`D1cut`); D2 is binned by `D2Tcut` as a treatment (second block of `"KRDY2"`/`"MWDY2"`, default as D1) and by `D2cut` as an outcome or conditioning variable W (default `D2Tcut` for the KR conditions, `"none"` for the MW conditions). Binning W changes power only; binning a treatment changes the exclusion restriction tested (Z may affect W only through the binned treatment). `montest()` is told not to bin again.
* Changed: `"KRDY2"` and `"MWDY2"` now bin D1 and D2 as treatments to `"equisized:4"` by default (they used them as passed), like `"KRD"`, `"KRDY"`, `"MWD"` and `"MWDY"`; results differ only for variables with more than four values (`D1cut = "none"`, `D2Tcut = "none"` restore the old behaviour). In `montest()` with `block`, D and Y may not need binning (`Dcut`/`Ycut = "none"`, or a count at least the number of distinct values).
* Fixed: `cluster` with a non-numeric-mode identifier (e.g. a string country code) errored
  ("Clusters must be able to be coerced to a numeric vector.", from `grf`'s own
  `clusters=` validation) only in the forest-fitting path used by `"simple"`/`"KR"`/`"AHS"`/
  `"MW"` nuisance and causal-forest fitting (`fit_models()`), which handed the raw cluster
  column straight to `grf::causal_forest`/`grf::regression_forest` instead of converting it
  first. Every other cluster-touching path in the package (`forest_test()`, `CART_test`,
  `fe_rank()`, `make_group_folds()`) already converts with `as.integer(factor(cluster,
  exclude = NULL))` before use; `fit_models()` now does the same, matching the existing
  convention instead of introducing a new one.
* Changed: all reductions on the path to forest inputs, cutoffs and test statistics (Hajek normalizers, propensity/outcome means, weighted quantiles, forest-search cumulative sums, cluster sums in the final test) now use order- and platform-independent summation (`ksum`, `kcumsum`, internal). R's `sum()`/`mean()`/`cumsum()` accumulate in `long double`, which differs between x86 (80 bit) and Apple Silicon (64 bit); forest/cutoff selection amplified those last-bit differences into different p-values. Results can differ from earlier versions in the last digits, and in rare cases in the selected subgroup. Bitwise agreement across platforms is still not guaranteed where `grf`'s C++ random number draws, BLAS, or compiler floating-point contraction differ.
* Added: `stabilize.scores = "hajek"` in `montest()`, a self-normalized (Hajek) replacement for the additive mean-shift of `Z.hat`, and it is now the default (`stabilize.scores = "hajek"`). The additive mean-shift of `Z.hat` (the previous default, `TRUE`) has been removed: `TRUE` is now an alias for `"hajek"` and `"shift"` is an error. Under `parametric = TRUE` the normalizers are pooled across the two sample halves. In simulations (n = 1500, binary Z, four propensity shapes) it was valid, and its power was equal or higher for KR and the simple first-stage test (up to 3 points), with no clear change for MW. For a binary instrument the arm weights `Z/e` and `(1-Z)/(1-e)` are each divided by their own weighted mean within every sample/margin cell; for a continuous instrument (or a binary one with fixed effects) the weight `(Z-e)/v` is divided by the cell mean of `(Z-e)^2/v`. Only the weight is normalized. On the test side the normalizers are re-estimated within the tested cell and the standard error uses the linearized influence function, with degrees of freedom reduced by the number of estimated normalizers; the `"MW"` contrast `Z/e - (1-Z)/(1-e)` is normalized the same way. Not applied to `doubly.robust = FALSE` with `target = "overlap"`. `FALSE` is also accepted as `"none"`.
* Fixed: a tested subset in which every score is exactly zero (zero mean and zero standard error, e.g. a multivalued-treatment `"MW"` subset whose pseudo-outcome vanishes) has `t = 0/0`; its p-value is now 1 instead of `NaN`, which was silently dropped from rejection counts.
* Changed: `minarm` and `minarm_clusters` are removed from `montest()` (passing them is an error) and replaced by `miness` (default 10) for all conditions. In the local (forest/CART) tests a candidate subgroup must now have effective sample size at least `miness` in each instrument arm, where the effective size of an arm is `(sum w)^2 / sum_c (sum_{i in c} w_i)^2` for the inverse probability weights `w` (`1/e` for `Z = 1`, `1/(1-e)` for `Z = 0`) and `c` the clusters. It equals the row count when the propensity is flat and is close to 1 when a few rows carry almost all of the weight, which row counts could not see. With a propensity near 0 or 1 for some rows, inverse probability weighted scores (especially `"MW"`, which has no outcome nuisance, and KR with `doubly.robust = FALSE`) over-rejected, in some designs up to 10-80% at the 5% level.
* Changed: `aipw.clip` now defaults to `"auto"`: `0.5 / sqrt(n / 2)` for `n` complete rows, between 1e-3 and 0.2 (about 0.022 at n = 1000), instead of 1e-3. An estimated propensity closer to 0 or 1 than that is mostly estimation noise at half-sample size, and its inverse dominates the weights. A numeric `aipw.clip` still overrides it (for the propensity and for the floor on the variance nuisance); `NULL` still means no clipping. The floor on the fitted variance nuisance `v(X)` (continuous instrument, fixed effects) stays at 1e-3 under `"auto"`. In simulations with propensities near 0 and 1, five propensity shapes, n = 500-4000 and binary and continuous outcomes, `miness = 10` with the auto clip gave rejection rates near or below the nominal level for `"MW"` and `"KR"` and cut the over-rejection of `"KR"` with `doubly.robust = FALSE` (28 of 40 cells above 10% before, 6 after); power at n = 1500 was unchanged or higher. Remaining weak spots: a rare instrument (propensity 0.05) at small n, where some tests have no eligible subgroup and error.
* Changed: `condition = "MW"` now follows `target`. Its scores are the plain inverse probability weighted contrasts `Q (Z/e - (1-Z)/(1-e))` instead of `e(1-e)` times that, so `target = "all"` (the default, as for KR) averages the local effects with equal weights, and `target = "overlap"` weights rows by `e(1-e)` (asymptotically the previous MW behaviour). Results for `"MW"` change under the default. The interior-value sets `{ghat > 0}`, `{ghat < 0}` are unaffected by the rescaling.
* Changed: with `priority = TRUE` the sort key is `pred / (w v(X))`, with `w` the row weight of the mean the search and test use (user weight and `target`), the Neill (2012) priority for the weighted t-statistic. Unchanged without weights and with `target = "all"`; rows with `doubly.robust = FALSE` keep the unweighted key.
* Added: condition `"FSS"` in `seqtest()`, the flipped first stage difference: the first stage sum `E[D1|Z=1]-E[D1|Z=0] + E[D2|Z=1]-E[D2|Z=0] >= 0`, tested like `"FSD"` for the constructed treatment `D1 + D2` (binary D1 and D2 only; not run by default). It is the first stage difference for D1 and 1 - D2 and has power when D2 moves against D1 among compliers, e.g. a programme that is rarely repeated.
* Added: conditions `"FSDY"` and `"FSSY"` in `seqtest()`, `"FSD"` and `"FSS"` localized in pre-specified bins of Y (binary D1 and D2, requires Y, not run by default). Each bin gives the summed discordant (`FSDY`) or concordant (`FSSY`) KRDY cells, i.e. the first stage of `(D1 - D2) 1{Y in b}` or `(D1 + D2 - 1) 1{Y in b}` is nonnegative; the bins are the `yval` margin, so `pool`/`select = "yval"` work across them, and pooling all bins gives FSD/FSS. Implemented by a new internal `kr_sums` argument of `montest()` (`condition = "KR"`, binary D, no `block`).
* Changed: `gridtypeY` in `montest()` now defaults to `"equisized"` (equal-count quantile bins), matching `gridtypeD` and `gridtypeZ`; it was `"equidistant"`. Pass `gridtypeY = "equidistant"` for the old behaviour.
* Changed: `montest()` no longer accepts `"sample"` in `select` (error; `select = "all"` now means every margin except `sample`). Selecting a sample direction ranks the two directions' candidates together, so whether a direction is carried forward depends on the other half's outcomes, which are the held-out data of the first direction. Leave `sample` in neither to test the directions separately, or in `pool`.
* Removed: the adaptive choice between pooling and selecting over a margin (the same margin in `pool` and `select`) and its outputs `$adaptive_dims`, `$adaptive_design`, `$adaptive_scores`, `$adaptive_pool`, `$adaptive_select`. It fitted every pool/select design and kept the winner on the training rows; simulations showed no clear gain over forced pooling or selection, and a choice made from both halves is not independent of the held-out tests. The code and a note are archived in `dev/adaptive/` and in the git tag `archive/adaptive-pool-select`. `forest_test_core()` is renamed `forest_test()`.
* Removed: the `studentize` argument of `montest()` (sorting on `pred / sqrt(pred_var)`); `priority` (default `TRUE`) replaces it. `studentize_te()` is gone.
* Changed: `montest()` drops `"sample"` from `pool` (with a message) whenever `select` contains a margin with at least two cells for the data and conditions (e.g. several conditions, or `select = "zmargin"` with a multivalued instrument). The cells selected by both sample directions were pooled across halves, but that choice uses both halves' outcomes, each of which is the held-out half of the other direction; the directions are now tested separately. Without effective selection, pooling the halves is unchanged.
* Added: `screen = "power"` in `montest()` (forest search, `stack = TRUE`): within each group of cells that differ only in the `select` margins, keep the top-j cells (by training t) that maximize the predicted probability that at least one is rejected by Holm on the held-out half. Held-out statistics are predicted from the training t-statistics (variance inflated by `1 + phi`, `options(montest.power_phi = 1)`) and the training-side correlation between cells, computed from cluster-level influence vectors; no prior. It is the default screen; with `testtype = "CART"` or `stack = FALSE`, where the influence vectors are not available, the default falls back to `"stepdown"` and an explicit `screen = "power"` is an error. Adds `mvtnorm` to Imports.
* Changed: defaults of `pool` and `select` in `montest()`. `pool` now defaults to `c("sample", "zmargin", "dval", "yval", "equation")` and `select` to nothing: `condition` and `block` are in neither, so each condition and block is tested on its own and corrected as one family (Holm), and the sample halves stay pooled. In simulations with `condition = c("simple", "KR")`, selecting across conditions lost 8-12 points of power to this default because selection forces the sample directions to be tested separately. A margin the user specifies in one of the two is dropped from the default of the other; a margin in both `pool` and `select` is an error.
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
