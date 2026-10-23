# MW over-rejection: Hajek weights, effective arm size, propensity clip (October 2026)

Status: the effective-size arm rule (`miness = 10`, replacing `minarm`/`minarm_clusters`) and the auto propensity clip (`aipw.clip = "auto"` = 0.5 / sqrt(n/2), clamped to [1e-3, .2]) were implemented in the package afterwards (see NEWS.md). The rest of this note is the investigation behind them; the centered-MW prototype and the `montest.arm_ess` / `montest.clip_rate` options exist only in the scratchpad copy.

## Setting
Boundary null: binary Z with propensity `clip(plogis(3.5 x), .02, .98)`, 3-valued D, valid IV, condition `"MW"`, `target = "all"`, n = 1500, 300 replications (Monte Carlo s.e. about 1.3 points), 5% nominal.

## Findings
1. **Cause.** The forest propensity hits the clip bound for a few rows, giving control rows weights near 1/aipw.clip = 1000. In half-samples of 750 rows the effective size of the control-arm weights, (sum w)^2 / sum w^2, is 3.7 and 10.8. The IPW test is then driven by one or two rows. Row counts (`minarm`) do not see this.
2. **Hajek / recentering.** Arm-normalised IPW (Hajek) and the package's with-intercept centered slope are numerically close (toy leaf: corr .98-.998, same s.d.). Both remove the dependence on the level of Q and are better than plain IPW with the true propensity, but both fail worse than IPW when the estimated propensity hits the clip (centered: 31% boundary rejection at old `minarm`, vs 9.7% for IPW).
3. **Centered test is valid only for X-only leaves.** For pointwise MW rows (binary D, endpoints) the forest uses Y as a feature, so leaves depend on Y. Arm normalisation then compares conditional probabilities, not the joint density change, and power collapses (endpoint power about .87 -> .27, same seeds). Interior rows (sets learned first, leaves in X) improve (3-category violation 0.56 -> 0.82).
4. **`stabilize.scores` is not the cause** (off: IPW still 10.7%).
5. **Effective-size arm rule** (`(sum w)^2 / sum w^2 >= k` per arm, in place of `minarm`/`minarm_clusters`): at the default clip it declares the test infeasible in 28-95% of replications (k = 10-40), which is the correct diagnosis. With a floor of 0.02 it is valid: boundary null 0.7% (k = 10), 0% (k = 20); interior power unchanged (MW .56, KR .62). (An earlier version of this note reported a drop in endpoint power, 0.89 -> 0.73; that compared different seed blocks and was wrong, see the recipe test below.)
6. **Clip level matters for IPW.** Floor .05: boundary rejection 13% (k = 10); IPW clip .02 with `minarm` 10/5: 14.7%. The floor .02 equals the design's true minimum propensity, so it is an oracle choice. The centered test is insensitive to the floor (0% at .05).

## Repo default (kept)
MW `minarm = minarm_clusters = 20`, clip 1e-3: boundary null 4.6% (500 replications), endpoint power about .86 (.855 in 400 replications; .89 in an earlier 100-replication block), interior (3-category) power .56.

## Recipe test: effective-size rule + floor c0 / sqrt(n/2) (scratchpad copy, option `montest.arm_ess = k`, `montest.clip_rate = c0`)
Boundary null, 300 replications, rejection at 5%:

| n | MW k=10, c0=1 | MW k=20, c0=1 | KR k=10, c0=1 | KR k=20, c0=1 |
|---|---|---|---|---|
| 800 | 3.3% | 1.3% | 0.3% | |
| 1500 | 5.0% | 1.0% | 2.7% | 0.7% |
| 4000 | 4.3% | 0.3% | | |

MW, n = 1500, k = 10, c0 = 0.5: 2.7%. Floor values used: .05 (n=800), .037 (n=1500), .022 (n=4000).
Power (n = 1500, identical seeds, 400 replications, endpoint violation): repo default .855; scratchpad with nothing on .855; effective-size rule only .853; floor only .855; both .853. Interior designs (3-category violation): MW .56, KR .62, unchanged; nulls 0%. So the recipe is about as valid and as powerful as the MW `minarm = 20` default in this design, with one threshold instead of two counts plus an MW special case. k = 10 is close to nominal (5.0% at n = 1500); k = 20 is conservative.
Caveat: one DGP; the fixed clip .05 failed at n = 1500 (13%) but c_n = .05 at n = 800 did not, consistent with clipping bias that grows with n unless the floor shrinks.

## Broader designs (scratchpad suite; October 2026)
Designs: propensity shapes const (.5), moderate plogis(1.5x), rare (.05), non-monotone plogis(1.5(x^2-1)), boundary (as above); n = 500, 800, 1500, 4000; binary and continuous Y; valid IV (size) and an exclusion violation at D = 2 (power, n = 1500). 60-150 replications per cell (Monte Carlo s.e. 2-4 points). Settings: repo defaults; recipe = effective-size arm rule k = 10 plus floor c0 / sqrt(n/2) with c0 = 1; half = same with c0 = 0.5. 40 comparable cells per condition.

| Condition | Setting | Mean rejection | Cells > 10% | Worst |
|---|---|---|---|---|
| MW | repo / c0=1 / c0=.5 | 0.8 / 1.1 / 0.6% | 0 / 0 / 0 | 8 / 8 / 5% |
| KR | repo / c0=1 / c0=.5 | 2.6 / 2.6 / 1.5% | 1 / 3 / 1 | 12 / 32 / 11% |
| KR, doubly.robust = FALSE | repo / c0=1 / c0=.5 | 26 / 5.3 / 4.6% | 28 / 8 / 6 | 78 / 25 / 23% |

- The repo default for KR with `doubly.robust = FALSE` over-rejects badly (up to 78%) whenever the propensity depends on x (moderate, non-monotone, boundary), a pre-existing problem; the effective-size rule plus floor removes most of it. Remaining problems: rare instrument (15-25%) and boundary with continuous Y at n = 4000 (23%).
- c0 = 1 breaks KR with a rare instrument at small n (32% at n = 500, continuous Y; floor .063 clips a true propensity of .05); c0 = 0.5 does not (8.5%).
- Rare instrument at n = 500: 21-26% of KR and MW runs error because no subgroup reaches effective size 10 (repo: KR 12-13%, MW 99%).
- Power (n = 1500): MW and KR equal or higher than repo defaults (boundary, binary Y: MW .25 vs .21, KR .63 vs .58). KR with `doubly.robust = FALSE` is lower than the repo in x-dependent designs, but the repo is invalid under the null there.
- Priority (60 replications, selected cells): not the cause of the MW/KR problems. For KR with `doubly.robust = FALSE` under repo defaults, `priority = FALSE` lowers over-rejection (43 -> 32, 35 -> 18, 65 -> 45%) without removing it; with the recipe the difference is 0-7%. The key for those rows is unweighted while the centered search weights by 1/v, a plausible mismatch (not tested). Clipping precedes both the scores and the priority variance (`stabilize_zhat()` clips Z.hat in place), so the order of operations is consistent.
- Not tested: cluster-aware effective size, c0 below 0.5, k other than 10, other outcome or treatment types.

## Open
- Making endpoint/binary-D MW rows set-based (learned sets, X-only leaves) would let the centered test, or a doubly robust score, apply to every MW row.
- One design, one n; not checked: KR with `doubly.robust = FALSE` in the boundary design.
