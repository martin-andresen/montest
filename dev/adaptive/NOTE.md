# Adaptive choice between pooling and selecting over a margin — archived (October 2026)

Removed from `montest()` after being disabled by default. A principled version may be worth revisiting.

## What it was
A margin listed in both `pool` and `select` meant "choose adaptively". `forest_test()` fitted every pool/select combination of those margins (2^k designs, `forest_test_core()` each), scored each design on the train rows, and kept the winner. The score was the smallest Holm-adjusted one-sided p-value of the screen survivors, `log(m) + min log p` (later variant: `m` = all candidate train cells; see the git history around `design_score()`), with a tie-break toward pooling. Outputs: `$adaptive_dims`, `$adaptive_design`, `$adaptive_scores`, `$adaptive_pool`, `$adaptive_select`.

- Code: `forest_test_adaptive.R` (verbatim legacy wrapper). Last package commit containing it: git tag `archive/adaptive-pool-select`.
- Prototype, never in the package: `weights_prototype.R` — train-side power-optimal weighted Bonferroni over the candidate cells of all designs (the copy it came from is `../montest_wts`, outside the repo; it can be deleted once this folder is confirmed).
- Simulation scripts: `sims/adaptive_choice_sim.R`, `sims/weights_vs_forced.R`. Results: `results/*.rds`; the paper and its extracted text: `refs/`. (`*.rds`, `*.pdf`, `*.txt` are in `.gitignore`, so `results/` and `refs/` live only on this machine unless added with `git add -f`.)
- Still in the package and independent of this feature: `screen = "power"` (greedy selection of the set of cells that maximizes the predicted held-out rejection probability, using influence-vector correlations, `train_u`, `cell_corr_from_u()`, `power_union()`).

## Why it was removed
1. **No clear gain.** Rejection rates at 5% (Holm min p), condition `MW`, 5-valued Z, binary D, margins = 3 zmargins x 2 equations, n = 3000, 200 replications, violations among x > 0 only:
   - concentrated (violation in one z-step): forced pool 0.31-0.33, forced select (stepdown screen) 0.275-0.40, old adaptive 0.315-0.37, weighted Bonferroni 0.25;
   - diffuse: forced pool 0.13-0.23, forced select 0.08-0.13, old adaptive 0.08-0.13, weighted Bonferroni 0.065;
   - null: at most 1% rejections for every arm (too conservative to detect size distortion).
   The adaptive default never beat the better fixed design. Monte Carlo s.e. is about 0.03, so many differences are within noise.
2. **Leakage across sample halves.** A design chosen from the train rows of both halves is not independent of the held-out tests: each half is the held-out half of the other direction. A single common design that is leak-free would need a design independent of both halves' outcomes (a third fold). Per-direction designs/weights with a fixed budget (alpha/2 each) and sample tested separately are leak-free.
3. **Sample halves can be pooled only without effective selection.** Cells selected by both directions were pooled across halves; that agreement uses both halves' outcomes. `montest()` therefore now drops `sample` from `pool` (with a message) whenever `select` contains a margin with at least two cells.
4. **Pooling is itself a row-level screen.** Forests are fit per margin cell and per half, and a pooled cell sorts stacked rows by cell-specific predicted effect with one cutoff, so cells with nothing are largely excluded already. The gap between pooling and selecting is therefore small except when the violation is diffuse (pooling wins) or concentrated with a strong signal (selecting wins).

## Theory notes
- Rubin, Dudoit & van der Laan (2006), `refs/rubin_dudoit_vanderlaan_2006.pdf`: optimal cutoffs c_m = log(lambda)/d_m + d_m/2 maximizing expected true positives under expected-false-positive control; estimating the shifts from the same data as the test collapses to the common cutoff; sample splitting (estimate cutoffs on a hold-out) can lose to the common cutoff unless shifts are pooled/smoothed.
- Dobriban, Fortney, Kim & Owen (2015, Biometrika, arXiv 1504.02935): Bayes p-value weights maximizing the expected number of rejections of weighted Bonferroni (flat prior on the means, variance inflated by 1 + phi, phi = 1 default); weighted Bonferroni controls FWER under any dependence; the problem is non-convex but solvable. It does not model dependence between tests and optimizes expected rejections, not P(at least one rejection).
- Our own objective: P(at least one held-out rejection) under a correlated predictive distribution N(-t, (1 + phi) R), R from cluster-level influence vectors; for any-rejection the Holm step-down adds nothing over Bonferroni at the first rejection.

## Findings from the weighted-Bonferroni prototype
Equal weights over the candidates did worse than optimized weights in the concentrated scenario (0.20 vs 0.25), more shrinkage (phi = 4) moved the weights toward equal weights (0.205), and neither beat forced pooling. Confound left open: the weights arm tests sample halves separately while the forced arms pool them; the missing arms are forced pool/select with sample separate.

## Ideas not tried
- Continuous weights inside `screen = "power"`; pooled cell as its own candidate in the same family (needs covariances between a pooled cell and its subcells via cluster-summed scores at each cell's own cutoff).
- Sample-separate forced-pool and forced-select benchmarks in the MW setup.
- Pooling/smoothing the shift estimates across cells (Rubin et al.'s suggestion); modelling the upward bias of the best-cutoff t-statistic.

## Restoring
Check out the tag (`git show archive/adaptive-pool-select:R/functions.R`), or paste `forest_test_adaptive.R` back as the wrapper around `forest_test_core()` and relax the overlap error in `montest()`. For the weights prototype also add `defer`/`ret_u` arguments to the wrapper's `run_core()` and the `p.wbonf` branch in `montest()` step 7.
