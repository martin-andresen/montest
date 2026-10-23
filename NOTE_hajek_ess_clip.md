# MW over-rejection: Hajek weights, effective arm size, propensity clip (October 2026)

Status: investigation only. Nothing here is in the package except the MW default `minarm = minarm_clusters = 20`. The experimental code lives outside the repo (scratchpad copy `montest_proto`: centered MW test via `options(montest.mw_center)`, effective-arm-size rule via `options(montest.arm_ess)`).

## Setting
Boundary null: binary Z with propensity `clip(plogis(3.5 x), .02, .98)`, 3-valued D, valid IV, condition `"MW"`, `target = "all"`, n = 1500, 300 replications (Monte Carlo s.e. about 1.3 points), 5% nominal.

## Findings
1. **Cause.** The forest propensity hits the clip bound for a few rows, giving control rows weights near 1/aipw.clip = 1000. In half-samples of 750 rows the effective size of the control-arm weights, (sum w)^2 / sum w^2, is 3.7 and 10.8. The IPW test is then driven by one or two rows. Row counts (`minarm`) do not see this.
2. **Hajek / recentering.** Arm-normalised IPW (Hajek) and the package's with-intercept centered slope are numerically close (toy leaf: corr .98-.998, same s.d.). Both remove the dependence on the level of Q and are better than plain IPW with the true propensity, but both fail worse than IPW when the estimated propensity hits the clip (centered: 31% boundary rejection at old `minarm`, vs 9.7% for IPW).
3. **Centered test is valid only for X-only leaves.** For pointwise MW rows (binary D, endpoints) the forest uses Y as a feature, so leaves depend on Y. Arm normalisation then compares conditional probabilities, not the joint density change, and power collapses (endpoint power 0.89 -> 0.27). Interior rows (sets learned first, leaves in X) improve (3-category violation 0.56 -> 0.82).
4. **`stabilize.scores` is not the cause** (off: IPW still 10.7%).
5. **Effective-size arm rule** (`(sum w)^2 / sum w^2 >= k` per arm, in place of `minarm`/`minarm_clusters`): at the default clip it declares the test infeasible in 28-95% of replications (k = 10-40), which is the correct diagnosis. With a floor of 0.02 it is valid: boundary null 0.7% (k = 10), 0% (k = 20); endpoint power 0.73-0.75 (vs 0.89 with the repo default); interior power unchanged (MW .56, KR .62).
6. **Clip level matters for IPW.** Floor .05: boundary rejection 13% (k = 10); IPW clip .02 with `minarm` 10/5: 14.7%. The floor .02 equals the design's true minimum propensity, so it is an oracle choice. The centered test is insensitive to the floor (0% at .05).

## Repo default (kept)
MW `minarm = minarm_clusters = 20`, clip 1e-3: boundary null 4.6%, endpoint power .89, interior (3-category) power .56.

## Open
- A floor that scales with n (e.g. about n^-1/2 per half), with the effective-size rule replacing the two counts, was not tested.
- Making endpoint/binary-D MW rows set-based (learned sets, X-only leaves) would let the centered test, or a doubly robust score, apply to every MW row.
- One design, one n; not checked: KR with `doubly.robust = FALSE` in the boundary design.
