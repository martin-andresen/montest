# Brief: Mourifié–Wan (MW) tests for multivalued ordered treatments in montest / seqtest

**For:** Claude Code, plan mode, in the local `montest` repo (github.com/martin-andresen/montest).
**Owner:** Martin Eckhoff Andresen.
**Status:** Theory and procedure agreed in a design session. Nothing has been implemented yet. Plan first, pilot small, then run anything large.

## 0. Working rules (user preferences)

- Ask clarifying questions before large changes. The open design questions are in Section 7.
- Pilot before running anything large: small n, few trees, one scenario.
- Do not kill interactive R or Stata sessions. Do not run broad `pkill R` or `taskkill` commands.
- The user is on Windows. Parallel code uses PSOCK clusters, never `mclapply`. Reinstall with `devtools::install()` or `load_all()` before simulating, and make sure no stale `montest` object in the global environment masks the package.
- Toy DGPs and small simulations are the preferred way to validate.

## 1. Context: what exists now

- `montest()` (R/montest.R) builds pseudo-outcomes `Q` per condition and margin cell and stacks them. It fits nuisance and causal forests via `fit_models()` / `estimate_conditional_mean()` in R/functions.R, then runs the sorted-cutoff search in `forest_test()` (train/test halves, `pool`, `select`, multiplicity corrections).
- `seqtest()` (R/seqtest.R) is a thin wrapper for sequential treatments: Z, then D1, then D2, then Y. Its conditions include FS, FSD, KR, KRD, KRDY2, MW, MWD, MWDY and MWDY2.
- **Kwon–Roth (KR)** already supports multivalued D. It stacks one row per (dval, outcome set A) cell. The interior-value bug is fixed: for `dval > dmin`, Q = 1{D ≥ dval} − 1{Y ∈ A, D = dval}, and with `joint = FALSE` it uses singletons plus their complements at interior values.
- **MW today is binary-only.** It uses the IPW scores s1 = D(Z/p − (1−Z)/(1−p)) and s0 = (1−D)((1−Z)/(1−p) − Z/p), with W (the outcome, or (D2, Y)) entered as a forest feature. Global MW (`local = FALSE`) reduces to the first stage, so it has no power by construction. Global MW was dropped from the simulation tables for this reason.
- Global KR calls need `pool = "none"`. Otherwise the margins are pooled and the test is degenerate.

## 2. Goal

Implement MW-type conditional tests for **ordered multivalued D ∈ {0,…,J}**, and binary or multivalued Z handled by Z margins as now. The tests should use **data-driven outcome sets** learned on the training half. The point is to get KR-strength conditions at interior values without enumerating subsets of W-cells: the cost is 2(J+1) pseudo-outcomes instead of one per (dval, A) cell.

## 3. Theory, compact (details in KR_conditions_Martin.tex)

The pieces:

- W is the outcome vector (Y, or (D2, Y) in seqtest).
- For value m, s_m = 1{D = m}·(Z/p − (1−Z)/(1−p)), and g_m(w, x) = E[s_m | W = w, X = x], which is the Z-induced change in the density of (D = m, W = w).
- π_{≥m}(x) = E[1{D ≥ m}·(Z/p − (1−Z)/(1−p)) | X = x], and likewise π_{>m}.

Under IV validity and monotonicity, for every m and x:

- **Inflow:** E[max(g_m, 0) | X] ≤ π_{≥m}(X).
- **Outflow:** E[max(−g_m, 0) | X] ≤ π_{>m}(X).

At m = 0 the inflow bound becomes the sign restriction g_0 ≤ 0. At m = J the outflow bound becomes g_J ≥ 0. These are the binary MW conditions. At interior m the bounds are budgets, not sign restrictions, and they are equivalent to KR over all sets B of W-cells (the sup over B is attained at B+ = {g_m > 0} and B− = {g_m < 0}).

**Pseudo-outcomes with estimated sets** B̂+ = {ĝ_m > 0} and B̂− = {ĝ_m < 0}:

- Q_in = 1{D ≥ m} − 1{D = m}·1{W ∈ B̂+}. The test is that E[Q_in · IPW-contrast | X] ≥ 0.
- Q_out = 1{D > m} + 1{D = m}·1{W ∈ B̂−}. The test is that E[Q_out · IPW-contrast | X] ≥ 0.

Any fixed set gives a valid test, because these are KR moments for a particular B. Estimating the sets only affects power, **as long as B̂± is learned on the training half only**. Never compare the plug-in E[max(ĝ, 0)] with π directly: that is biased upward and over-rejects.

At the endpoints, Q_in at m = 0 and Q_out at m = J reproduce the existing binary MW test with a learned set. Keep the existing binary MW path as the default for J = 1 unless the regression check (Section 6) shows the two paths are identical.

## 4. Procedure per value m (the agreed design)

1. **Learn g_m.** On the training half, fit one regression forest of s_m on (W, X), using inner folds for out-of-bag or cross-fit predictions. Predict ĝ_m for every unit in both halves; test-half predictions come from the training-half model. One forest per m serves both sides; there is no separate inflow and outflow fit.
2. **Stack by side.** Duplicate rows with `side ∈ {in, out}` as a new margin, like `equation` for binary MW. Drop the inflow side at m = 0 if trivial and the outflow side at m = J (check which ones are trivial). Cluster on the unit `id`, because rows are stacked.
3. **Build Q.** Use Q_in or Q_out from Section 3 with B̂± from step 1.
4. **Fit the causal forest** of Q on Z, **with X only** as features. W is now absorbed in the set, so it must not be a feature. Use the same nuisance machinery as KR rows, including the doubly robust scores.
5. **Run the existing `forest_test()` cutoff search**, then pool or select across (m, side, zmargin, treatment) with the existing corrections (Holm, CCT and so on).

Optional threshold: B̂+ = {ĝ_m > c·se} to reduce noise in the set. Make it a tuning option with default c = 0.

## 5. seqtest integration

- MW, MWD, MWDY and MWDY2 should accept multivalued D1 and D2 through the same machinery: W = Y for D2 and W = (D2, Y) for D1.
- MWDY2 (use case 2, the first-treatment conditions) with multivalued D1 should become the "data-driven-set KR" for D1 with W = (D2, Y).
- Note from the theory: the use case 2 intersection of conditions is not sharp when D1 is multivalued (there is an LP counterexample). This is not a blocker for implementation; it is a documentation note.

## 6. Validation plan

1. **Unit tests** (tests/testthat):
   - Q_in and Q_out on a hand-built toy table.
   - Endpoint reduction: at m = 0 and m = J the new path gives the same Q as the binary MW with the same set.
   - With J = 1, results are numerically close to the current MW.
2. **Size** (Panel B DGPs in sim_KR.R): MV0_tracks, MV0_boundary, MV1 and MV2. Rejection rates should be at or below 5% (FSD over-rejects in MV1 and MV2 and KR does not; the new MW should behave like KR).
3. **Power:** MV3, MV3_X and MV3_D2, compared against local KR. The expectation is that power is similar or better, at a fraction of the time.
4. **Timing:** add a "seconds per test" row, comparing the new MW against KR with `joint = TRUE` and `joint = FALSE`.
5. Pilot each step with a small n and few reps before the full run.

## 7. Open design questions (ask the user before deciding)

1. Where to fit ĝ_m: inside `montest()` before stacking, cached per m, or as a new helper in functions.R?
2. Should `side` be its own margin, or reuse `equation`?
3. Option names, for example `mw_sets = c("datadriven", "binary")`, `mw_threshold = 0`, `mw_inner_folds`, `mw_priority`.
4. Should step 1 use doubly robust scores instead of IPW for s_m when there are X or fixed effects? (MW with fixed effects is currently disallowed.)
5. Ranking in the cutoff search: keep the current ordering, or use the Neill (2012) optimal priority −τ̂(x)/σ̂²(x)? Scores may be reused only with honest splitting.
6. How does `local = FALSE` behave? A global version with a learned set is no longer degenerate. Decide whether to offer it (with `pool = "none"`).
7. Multivalued Z: keep the Z-margin stacking as it is?

## 8. Code locations

Line numbers come from an older snapshot, so verify them before editing.

- **R/montest.R:**
  - condition parsing and validation;
  - margin-index construction and the KR `A_specs` loop (~1870–2035);
  - joins with the stacked data;
  - Q construction (~2195–2220, next to the KR interior fix);
  - the calls to `fit_models`;
  - pool and select (~2834 and ~2873).
- **R/functions.R:** `estimate_conditional_mean`, `fit_models`, `forest_test`, `make_group_folds`, `binarize_var`.
- **R/seqtest.R:** test lists and guards (`test %in% ...`), MW Q construction, and expansion into treatment and MWD rows. The current MW block assumes binary D.
- **sim_KR.R:** Panels A, B and C, `run_mt`, `run_st` and `write_meat()`.

## 9. References

- Mourifié and Wan (2017, REStat)
- Kitagawa (2015, Econometrica)
- Kwon and Roth (2024)
- Sun (2023, JoE): multivalued D, endpoint conditions only
- Kédagni and Mourifié (2020, Biometrika)
- Farbmacher, Guber and Klaassen (2022, JBES)
- Neill (2012, JRSS-B): linear-time subset scanning
- Project note: KR_conditions_Martin.tex (the MW section and Proposition 3)
