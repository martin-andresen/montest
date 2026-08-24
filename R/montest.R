#' Monotonicity and LATE assumptions tests using sample splitting and machine learning
#'
#' \code{montest()} searches for violations of monotonicity and LATE assumptions in
#' data-adaptive subsets of the sample and across various margins using parametric and semiparametric methods. It combines sample splitting,
#' cross-fitting, and generalized random forest (or CART-based subset search) to identify
#' regions of the covariate space or margins of the instrument or treatment where test
#' statistics are most negative, and then evaluates those regions in the held-out sample.
#' The function supports several testable conditions, including a simple test of whether
#' the first stage is negative \code{"simple"}, the Kwan and Roth (2026)/Sun (2023) conditions \code{"KR"} (collapsing to Balke and Pearl (1997) in the case of binary instrument and treatment),
#' the Mourifie and Wan (2017) conditions  \code{"MW"}, and the first stage conditional on Y-test from
#' Andresen, Huber and Sloczynski (2026) \code{"AHS"}. Multivalued instruments, treatments and outcomes
#' are expanded across margins and discretized into bins before estimation or treated continuously, depending on options. See Details.
#'
#' @param data A \code{data.frame} or \code{data.table} containing the analysis sample.
#'   Observations with missing values in any variables used by the call are dropped.
#' @param fml A fixest-style formula on the form Y~X|fe|D~Z that specifies the IV call to be tested. Importantly, the instrument Z should be coded so that higher values of Z imply weakly higher values of D for everyone - positive monotonicity. By default, the X part of the formula is used for both nuisance estimation and heterogeneous effects. X may contain the same variables as FE, indicating that some FE variables should also be used for heterogeneity.
#' @param fml.Z Optional: A one-sided formula for the nuisance of the instrument. Defaults to the same as the the general formula in fml
#' @param fml.varZ Optional: A one-sided formula for the conditional-variance nuisance
#'   \eqn{v(X)=Var(Z|X,FE)} alone (see \code{doubly.robust}/\code{target} for when this is
#'   estimated), separately from \code{fml.Z}'s covariates for \code{Z}'s own conditional
#'   mean -- the two need not use the same covariates (e.g. a more parsimonious variance
#'   model for stability). Defaults to \code{fml.Z} (which itself defaults to \code{fml}'s
#'   main X part), so leaving it unset reproduces the previous behavior of reusing the same
#'   covariates used for \code{Z}'s conditional mean.
#' @param fml.Q Optional: A one-sided formula used for the nuisance of the pseudo-outcome. Defaults to the same as the the general formula in fml
#' The formula may be one-sided and omit Y if testing only the simple first stage condition. Note that the exact functional form does not matter in the default case when \code{parametric=FALSE} because the command uses semiparametric methods.
#' @param parametric A boolean indicating whether nuisances should be estimated using the parametric functional form specified or using semiparametric methods (the default). In the latter case,
#' all fixed effects are residualized as specified in the FE part of the formula, while the functional form in the main part of the formula is ignored and determined by the corresponding regression forests for the nuisance parameters.
#' \code{parametric} governs nuisance estimation only.
#' Nuisance fragility that can arise under \code{parametric=TRUE} (e.g. a linear-probability-model propensity falling outside \eqn{[0,1]}) is caught at score-construction time by the validity checks
#' described under \code{aipw.clip}.
#' @param condition Character vector selecting which tests to run. Allowed values are any combination of
#'   \code{"simple"}, \code{"KR"} (Kwan-Roth conditions), \code{"MW"} (Mourifie and Wan), \code{"AHS"} (Andresen-Huber-Sloczynski), or \code{"all"}.
#'   If \code{Y} is omitted from the formula, only \code{"simple"} is allowed.
#' @param target a scalar equal to either "all" (default) or "overlap" to determine whether the function should target the average treatment effect (all) in the relevant subsample, or the overlap-adjusted / weighted ATE ("overlap").
#'   Applies to \code{"simple"}, \code{"KR"}, and \code{"AHS"} (not \code{"MW"}, whose own moment is untouched by \code{target} -- see \code{doubly.robust}), across both instrument representations (binary or
#'   continuous). These are genuinely different estimands whenever the underlying effect is heterogeneous in X: each row's score has conditional expectation equal to the *local* effect at that row's X, so
#'   \code{"all"} targets the plain (unweighted) average of those local effects, while \code{"overlap"} targets their variance-weighted average. When \code{doubly.robust=TRUE}, \code{target} changes only the
#'   weight used when aggregating scores into the final test statistic, not the per-row score formula itself. When \code{doubly.robust=FALSE}, \code{target} additionally selects which conditional-variance
#'   estimator normalizes each row's score (a fitted per-observation model for \code{"all"}, versus the raw empirical \eqn{(Z-\hat Z)^2} for \code{"overlap"}), so the per-row scores themselves differ between the two
#'   settings, not just their aggregation weight; \code{"overlap"} in that case coincides exactly with the classical Frisch-Waugh-Lovell / OLS partialling-out coefficient, at any level of aggregation (the full sample, or any subgroup subsequently averaged over).
#' @param doubly.robust Logical, default \code{NULL} (resolved from \code{parametric}: \code{TRUE}
#'   whenever \code{parametric = FALSE}, the overall default; \code{FALSE} whenever
#'   \code{parametric = TRUE}). If \code{TRUE}, scores are constructed using the doubly-robust
#'   (AIPW-style) moments, augmenting the residual term with the causal forest's predicted CATE.
#'   If \code{FALSE}, scores use the singly-robust/partialling-out (FWL) moment. "Singly robust" here means robust specifically to misspecification of the outcome nuisance
#'   \code{Q.hat}: the FWL score is consistent whenever the instrument's conditional mean \code{Z.hat} is correctly specified, for *any* \code{Q.hat} (which then only affects efficiency, not consistency) --
#'   but, unlike AIPW, offers no protection in the other direction. AIPW (\code{doubly.robust=TRUE}) is robust to
#'   misspecification of *either* nuisance, as long as the other is correct. AIPW's orthogonality protection matters most
#'   when nuisances are flexibly/ML-estimated (\code{parametric = FALSE}) -- under \code{parametric = TRUE} the nuisance
#'   functional form is trusted by construction, so AIPW buys comparatively little while still adding a second nuisance
#'   function's estimation noise, hence the conditional default above. Passing \code{doubly.robust} explicitly always
#'   overrides this regardless of \code{parametric}; see \code{target} for how the two interact when \code{doubly.robust=FALSE}.
#'   Applies to \code{"simple"}, \code{"KR"}, and \code{"AHS"}; has no effect on
#'   \code{"MW"}, which does not use this scoring machinery.
#' @param linearD Logical, default FALSE. If TRUE, the treatment D is scored on its raw/linear scale
#'   instead of being margin-stacked and binarized -- i.e. estimated with a linear, continuous,
#'   multivalued D rather than cut-and-stack. Only meaningful for \code{"simple"} and \code{"AHS"}, whose moment
#'   conditions support a slope-based representation; \code{"KR"} and \code{"MW"} are always indicator-based and
#'   always use the binarized representation, so \code{linearD = TRUE} is disallowed in combination with them. Downgraded to
#'   FALSE with a warning if D is binary in the data (no meaningful linear scoring in that case).
#' @param linearZ Logical, default FALSE. If TRUE, the instrument Z is scored on its raw/linear scale
#'   instead of being margin-stacked and binarized -- i.e. estimated with a linear, continuous,
#'   multivalued Z rather than cut-and-stack. Only meaningful for \code{"simple"} and \code{"AHS"}; \code{"KR"}/
#'   \code{"MW"} always use the binarized representation, so \code{linearZ = TRUE} is disallowed in combination
#'   with them. Downgraded to FALSE with a warning if Z is binary in the data.
#' @param inner.folds Optional integer giving the number of within-sample folds used for
#'   cross-fitting nuisance functions and, optionally, forest predictions. Set to
#'   \code{NULL} to disable the inner split. Defaults to NULL - nuisances and predictions from causal forests are fit out-of-bag. See option crossfit, which decides which parts this applies to.
#' @param crossfit Character vector of what parts of the procedure to cross-fit. Accepts "Z","Q","Y","C". If e.g. "Z" appears in crossfit, nuissances for Z are cross fit, either across outer sample part (if inner.folds==NULL), or within outer sample part across inner folds. If "Z" does not appear, OOB predictions are used. "C" is for the causal forest fit.
#' @param stabilize.scores Logical, default TRUE; applies a finite-sample stabilization
#'   to the estimated instrument nuisance, using whichever correction is appropriate for
#'   the representation actually in use:
#'   \itemize{
#'     \item Continuous representation (a genuinely multivalued instrument, \code{linearZ=TRUE},
#'       or any binary instrument estimated \emph{with} fixed effects -- see \code{aipw.clip}):
#'       \code{Z.hat} is mean-shifted within each sample/margin group so its group-average
#'       residual is exactly zero, correcting finite-sample/cross-fitting bias in the
#'       conditional-mean estimate itself. There is no IPW/arm structure for a continuous
#'       instrument to self-normalize, so this additive correction is the only one available.
#'     \item Binary representation, estimated \emph{without} fixed effects: \code{Z.hat} is
#'       clipped to \code{[aipw.clip, 1-aipw.clip]} (as before), but is otherwise left alone --
#'       instead, the closed-form conditional variance \eqn{v(X)=e(X)(1-e(X))} used in the score
#'       is rescaled per row by a Hajek/self-normalized-IPW constant: \eqn{c_1=}
#'       mean(1/e(X)) over treated rows and \eqn{c_0=}mean(1/(1-e(X))) over control rows,
#'       each computed within the same sample/margin group, so that
#'       mean(1/(e(X)*c_1))=1 exactly among treated rows and mean(1/((1-e(X))*c_0))=1 exactly
#'       among control rows. This targets the IPW weight \eqn{Z/e(X)-(1-Z)/(1-e(X))} that the
#'       score actually uses, rather than \code{Z.hat} itself -- a multiplicative correction on
#'       the scale the score is built on, replacing what used to be an additive mean-shift of
#'       \code{Z.hat} for this representation. Not applied when \code{doubly.robust=FALSE} and
#'       \code{target="overlap"} (see \code{target}), since that path already uses each row's own
#'       empirical variance instead of the closed form, and needs \code{Z.hat} untouched to
#'       reproduce the classical FWL/OLS coefficient exactly.
#'   }
#'   Both corrections are computed once, at whatever \code{sample}/\code{margins} level
#'   \code{stabilize.scores} operates on -- not separately within an adaptively-found
#'   leaf/subgroup nested inside that level, or within a "global" cell computed at a finer
#'   grouping. When \code{TRUE}, the FINAL (test-side and global) score average for AIPW
#'   (\code{doubly.robust=TRUE}) and the singly-robust score under \code{target="all"} redoes
#'   whichever correction applies -- the additive mean-shift, or the Hajek rescaling --
#'   using only that specific cell's own rows, the same finite-sample-bias correction
#'   \code{doubly.robust=FALSE & target="overlap"} already gets via its own (more thorough)
#'   centered/with-intercept construction (see \code{target}). This only re-corrects the
#'   propensity/IPW-weight term -- it does NOT re-estimate the causal forest's plug-in CATE
#'   prediction (\code{tau}) for that cell; if \code{tau}'s own average over exactly the rows
#'   selected into a subgroup carries local bias (e.g. from adaptive selection correlating
#'   which rows land there with their predicted effect), this does not correct for that.
#'   Never applied to the training-side adaptive search itself, only to the final evaluation,
#'   and never to cells that already use the \code{target="overlap" & doubly.robust=FALSE}
#'   centered path.
#' @param aipw.clip Positive scalar in \code{(0,1)}, default 1e-3, used to trim estimated propensity
#'   scores when augmented inverse-probability weighted scores are constructed and when stabilizing
#'   binary-instrument scores (see \code{stabilize.scores}).
#'   Also used, for the continuous representation (a genuinely multivalued instrument, \code{linearZ=TRUE}, or any
#'   binary instrument estimated with fixed effects), to floor the fitted conditional-variance nuisance \eqn{v(X)}
#'   at \code{aipw.clip} times a pooled variance scale. In both cases a structurally invalid value (a propensity
#'   outside \eqn{[0,1]}, or a variance below 0) triggers a hard error rather than being silently clipped, since
#'   that indicates nuisance misspecification; values merely close to the boundary are clipped with a warning.
#'   Exception: when \code{doubly.robust=FALSE} and \code{target=="overlap"}, the propensity is used only
#'   additively (never as a divisor) in the score for binary instruments, and is deliberately left unclipped so
#'   the resulting estimate reproduces the classical Frisch-Waugh-Lovell / OLS coefficient exactly -- see \code{target}.
#' @param drop_singletons Logical, default TRUE. If \code{TRUE} and fixed effects are
#'   present, observations belonging to a fixed-effect "singleton" (an FE level with only
#'   one observation -- possibly only after removing other singletons in a different FE
#'   dimension, under multiple FE) are dropped before any nuisance estimation. These
#'   observations are perfectly explained by their own FE dummy and contribute nothing to
#'   any regressor's FE-adjusted coefficient; per Correia (2016) and reghdfe's own default
#'   (\code{keepsingletons = FALSE}), leaving them in without a matching degrees-of-freedom
#'   correction biases cluster-robust standard errors downward. Coefficient point estimates
#'   are unaffected either way -- this only corrects inference. Implemented via fixest's own
#'   \code{fixef.rm = "singletons"} (iterative, multi-way-FE-aware) rather than a hand-rolled
#'   version.
#' @param drop_novar_Z Logical, default TRUE. If \code{TRUE} and fixed effects are present,
#'   observations in fixed-effect cells with no within-cell variation in the instrument
#'   \code{Z} are dropped before any nuisance estimation. Such a cell's residualized
#'   \code{Z} is exactly 0, so -- in the classical Frisch-Waugh-Lovell sense -- it
#'   contributes nothing to identifying the coefficient on \code{Z} or its variance; in the
#'   semiparametric pipeline these cells are additionally the primary source of the
#'   "conditional variance v(X) for continuous-Z rows" clip warnings described under
#'   \code{aipw.clip} (a variance model fit by pooling across FE cells inevitably leaks a
#'   small nonzero fitted deviation into a cell whose true variance is exactly zero). Unlike
#'   \code{drop_singletons}, this is specific to \code{Z} -- an analogous lack of variation
#'   in \code{D} (or \code{Y}) is not dropped, since such a cell still contributes a genuine
#'   zero-covariance data point to the classical estimator (correctly diluting a
#'   heterogeneous effect toward its population average) rather than being uninformative.
#' @param weight Optional character scalar naming a nonnegative weight variable.
#' @param cluster Optional character scalar naming a cluster identifier. Cluster-robust
#'   inference is used in forest-based testing. CART testing cannot be combined with cluster.
#' @param seed Integer random seed, default 10101, set to NULL to disable setting seed
#' @param minsize Integer minimum effective sample size or minimum cluster count required
#'   for subset search and testing. Default 50.
#' @param shrink Shrink predicted treatment effects using empirical bayes before sorting. Default 0: No shrinkage. 1: Full shrinkage
#' @param gridtypeY,gridtypeD,gridtypeZ Character strings controlling how continuous
#'   variables are discretized before stacking. Must be one of \code{"equisized"} or
#'   \code{"equidistant"}.
#' @param Ysubsets,Dsubsets,Zsubsets Integers giving the number of bins used when
#'   discretizing \code{Y}, \code{D}, and \code{Z}, respectively. Integer >1 required.
#' @param Y.res Logical; if \code{TRUE}, outcomes are residualized from X before tests that use
#'   outcome on the right hand side \code{MW,AHS}.
#' @param testtype Character string selecting the subset-search routine. Must be
#'   \code{"forest"} or \code{"CART"}. Default: Forest.
#' @param fe_rank_adj Logical; if \code{TRUE} (default), degrees of freedom for the
#'   final cluster-robust test statistics are reduced by the rank of the fixed effects
#'   (and, under \code{parametric = TRUE}, the rank of the parametric nuisance
#'   covariates) actually used within each tested cell. This always uses an exact
#'   computation and affects every reported statistic, regardless of \code{testtype}.
#' @param fe_rank_conservative Logical; only relevant when \code{testtype = "forest"}.
#'   Controls the cheap, approximate running fixed-effects-rank formula used purely to
#'   screen candidate cutoffs during the forest-based search (recomputing the exact
#'   rank at every candidate would be too slow). \code{TRUE} (default) sums each fixed-
#'   effect term's running category count independently, which over-counts (and is
#'   therefore conservative -- smaller effective degrees of freedom, wider intervals)
#'   whenever multiple fixed-effect terms share or interact categories. \code{FALSE}
#'   applies a cruder overlap correction that can under-count in the same situation.
#'   This setting never affects the final reported statistic for whichever cutoff is
#'   selected -- that value always comes from the exact computation controlled by
#'   \code{fe_rank_adj} -- only which cutoff the search prefers. Has no effect at all
#'   under \code{testtype = "CART"}, whose own cutoff search already uses the exact
#'   computation.
#' @param gridpoints Optional integer controlling the number of candidate cutoffs searched
#'   by the forest-based test. If \code{NULL}, all eligible cutoffs are considered.
#' @param min_n Integer minimum number of treated and untreated instrument observations
#'   required within each sample half and margin cell considered for residual variation.
#' @param pool Character vector controlling which dimensions are pooled when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dval"},
#'   \code{"yval"}, \code{"condition"}, \code{"equation"},
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. Margins (except "sample") can appear in both \code{pool} and \code{select}, implying adaptive pooling. Relevant margins that appear in neither are all tested, and tests are corrected for multiple hypothesis testing.
#' @param select Character vector controlling which dimensions are selected over when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dval"},
#'   \code{"yval"}, \code{"condition"}, \code{"equation"},
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. Margins (except "sample") can appear in both \code{pool} and \code{select}, implying adaptive pooling. Relevant margins that appear in neither are all tested, and tests are corrected for multiple testing.
#' @param screen Screening rule for deciding what determines a "promising" leaf or cell to carry forward to testing. May be "minimum","negative","nonpositive","stepdown","fg_relevant","none". Defaults to stepdown, described below.
#' @param cp,maxrankcp,alpha,prune Tuning parameters for the CART-based search
#'   routine. See Details.
#' @param Zparameters,Yparameters,Qparameters,Cparameters,Rparameters Named lists of
#'   additional arguments passed to the underlying estimation routines for different
#'   nuisance or target models. See regression_forest, causal_forest, feols and rpart for for details.
#' @param joint specifies that all Kwan-Roth conditions should be included in the test, not only those for which the subset A contains only one outcome value. Defaults to TRUE.
#'
#' @details
#' The rough steps of the montest algorithm is as follows
#'
#' \enumerate{
#'   \item Data is restricted to complete cases of the variables used, and the sample is split into two halves (optionally respecting clusters).
#'   \item Continuous or multivalued instruments, treatments and outcomes discretized into bins.
#'   \item The data is stacked across instrument margins, treatment margins, outcomes,
#'   equations, and test conditions. The outcome variable Q is defined, depending on condition and margins.
#'   \item Nuisance functions for as \code{Z.hat}, and the outcome \code{Q.hat} are estimated
#'   \item Separate causal forests of the outcome \code{Q}, on
#'   the instrument \code{Z} using features \code{X} (and optionally \code{Y} for MW and AHS conditions),
#'   treatment effects are predicted in and out of sample and scores constructed
#'   \item Each sample part (optionally within margins, depending on the options in \code{pool}) is sorted
#'   according to treatment effects, and the mean of scores is estimated numerically for all possible cutoffs
#'   in predicted treatment effects. Select the cutoff with the smallest t-statistic on the mean of scores
#'   Alternatively, subset selection can be done using a CART algorithm.
#'   \item Depending on the choices in \code{select}, promising subsets (using selection rules in \code{screen}) are evaluated  in the opposite sample half. If performing multiple tests,
#'   depending on the option \code{pool} and \code{select}, p-values are adjusted for multiple testing.
#' }
#'
#' \code{montest} supports weights, clustering, and fixed effects and allow for multivalued treatments and instruments
#' by binarizing instruments and treatments into quantile or equisized bins or treating them continuously. The command also tests for
#' one-sided monotonicity (within margins of the treatment and instrument) and if found, warns the user
#' and skips testing any trivially satisfied conditions.
#'
#' Consider a multivalued instrument with K values, a multivalued treatment with J values, a multivalued outcome with L values. \code{montest} can test the following families of conditions:
#' \enumerate{
#'  \item \code{condition="simple"} tests the sharp condition of a nonnegative first stage everywhere, which tests monotonicity,
#'   but not exclusion in addition to instrument exogeneity. There are a total of (J-1)(K-1) such conditions.
#'  \item \code{condition="KR"} tests the sharp Kwan-Roth (2026) / Sun (2023) conditions for instrument validity,
#'  which require exclusion and monotonicity in addition to instrument exogeneity. There are in total (K-1)(L + (J-1)(2^L-2))
#'   such conditions. If option joint=FALSE, montest tests only the (K-1)JK cellwise conditions with only one outcome value in the set A.
#'   With a binary treatment and instrument, these conditions collapse to the conditions from Balke and Pearl (1997). Outcome variable Y must be specified.
#'   \item \code{condition="AHS"} tests the non-sharp condition from Andresen-Huber-Sloczynski of a nonnegative
#'   first stage conditional on Y, which require monotonicity and exclusion in addition to instrument exogeneity.
#'   There are a total of (J-1)(K-1) such conditions.
#'   \item \code{condition = "MW"} tests the sharp conditions from Mourifie and Wan (2017), which tests monotonicity
#'   and exclusion conditional on instrument validity. This is only allowed for a binary treatment and a genuinely
#'   binary instrument (Z with at most 2 support points in the data -- an arbitrary margin/threshold cut of a
#'   multivalued Z does not correspond to their theory), and \code{fml} may not include fixed effects (MW's Q
#'   statistic uses Z.hat directly as a propensity weight with no AIPW/FWL orthogonalization protecting it).
#'   \code{linearZ}/\code{linearD} are disallowed in combination with \code{"MW"} for the same reason: its moment
#'   conditions are always indicator-based, never slope-based. There are a total of 2K such conditions.
#' }
#'
#' \code{parametric} and \code{doubly.robust} are orthogonal arguments (see their own entries above),
#' but not every combination is equally well-justified statistically:
#' \enumerate{
#'   \item \code{parametric=FALSE}, \code{doubly.robust=TRUE} (the default for both): the best-justified
#'   pairing for flexible nuisances. Cross-fit forests converge slower than \eqn{\sqrt{n}} and can be
#'   locally misspecified; a Neyman-orthogonal (doubly-robust) score is what keeps that first-order
#'   nuisance-estimation error out of the final estimate, letting it enter only at second order instead
#'   of biasing the result. Recommended unless there is a specific reason to deviate.
#'   \item \code{parametric=TRUE}, \code{doubly.robust=FALSE}: sensible when the nuisance functional form
#'   is genuinely trusted to be (close to) linear. A correctly-specified linear model does not have the
#'   slow-convergence problem doubly-robust scoring exists to protect against, so the AIPW correction adds
#'   little while a doubly-robust score built on in-sample, non-cross-fit nuisances loses the finite-sample
#'   protection cross-fitting would otherwise provide.
#'   \item \code{parametric=FALSE}, \code{doubly.robust=FALSE}: still consistent (the FWL score only
#'   requires \code{Z.hat} to be correctly specified, for any \code{Q.hat} -- see \code{doubly.robust}),
#'   but gives up both the orthogonality protection and the efficiency gain from using the causal forest's
#'   predicted CATE as a control variate, for no compensating benefit. Rarely worth choosing deliberately.
#'   \item \code{parametric=TRUE}, \code{doubly.robust=TRUE}: the most fragile combination. Nuisances are
#'   fit in-sample without cross-fitting, so the doubly-robust correction loses much of its usual
#'   theoretical justification (no sample-splitting protects against the same-sample fit's own overfitting),
#'   while keeping the fragility risk of an unconstrained parametric propensity (e.g. a linear-probability-
#'   model prediction outside \eqn{[0,1]}). This combination is not disallowed -- the concrete failure mode
#'   is caught by the validity checks described under \code{aipw.clip} (a hard stop on structurally invalid
#'   nuisance values, or a clip-and-warn near the boundary) -- but it is the pairing to reach for last, and
#'   only when the linear nuisance model is trusted to be well-specified with propensities safely bounded
#'   away from 0 and 1.
#' }
#'
#'
#' @return
#' A named list:
#'
#' \describe{
#'   \item{\code{results}}{A Matrix of train and test results by sample and margin cell.}
#'   \item{\code{margins}}{A Matrix of margins of condition, zmargin, dval, yval, equation that characterizes the cells used for testing.}
#'   \item{\code{global}}{Global mean test statistics over the same stratification.}
#'   \item{\code{grid}}{Stored cutoff grid evaluated by the forest test.}
#'   \item{\code{Xmeans}}{Weighted covariate means in the selected testing subset.}
#'   \item{\code{Xmeans_all}}{Weighted covariate means in the full sample.}
#'   \item{\code{XSD}}{Weighted covariate standard deviations in the full sample.}
#'   \item{\code{shares}}{Shares of the testing subset and full sample across pooled
#'   margins, when applicable.}
#'   \item{\code{minp}}{Minimum adjusted p-values and the Cauchy combination p-value when
#'   multiple hypotheses are tested.}
#'    \item{\code{time}}{A matrix of elapsed user, system, and total time by stage.}
#'   \item{\code{obs}}{A vector containing the number of observations \code{N} and, when
#'   clustering is used, the number of clusters \code{G}.}
#' }
#'

#'
#' @examples
#' \dontrun{
#'
#' #Generate data - Simulation DGP from Farbmacher et. al.
#' data=fct_datasim(setup="A",dgp=2,n=3000)
#'
#' # Simple monotonicity-style test
#' out <- montest(
#'   ~Xvar1+Xvar2+Xvar3|D~Z,
#'   data = data,
#'   condition = "simple")
#'
#' # Test multiple conditions, pooling evidence
#' out2 <- montest(
#'   Y~Xvar1+Xvar2+Xvar3|D~Z,
#'   data = data,
#'   X = c("Xvar1", "Xvar2", "Xvar3"),
#'   condition = c("simple","KR","MW"))
#' }
#'
#' @seealso montestplot LATEtest
#' @export

montest=function(fml,data,fml.Z=NULL,fml.Q=NULL,fml.varZ=NULL,condition=NULL,inner.folds=NULL,crossfit=NULL,
                 stabilize.scores=TRUE,aipw.clip=1e-3,drop_singletons=TRUE,drop_novar_Z=TRUE,weight=NULL,cluster=NULL,seed=10101,minsize=50L,
                 gridtypeY="equidistant",gridtypeD="equisized",gridtypeZ="equisized",stratify=TRUE,joint=TRUE,
                 Ysubsets = 4L, Dsubsets = 4L,Zsubsets=4L,Y.res=TRUE,testtype="forest",fe_rank_conservative=TRUE,fe_rank_adj=TRUE,
                 gridpoints=NULL,min_n=1L,pool=NULL,select=NULL,shrink=0,linearD=FALSE,linearZ=FALSE,target=NULL,
                 doubly.robust=NULL,
                 cp=0,maxrankcp=10L,Rparameters=list(),alpha=0.05,prune=TRUE,screen="stepdown",parametric=FALSE,
                 Zparameters=list(),Yparameters=list(),Qparameters=list(),Cparameters=list()
){

  ## Captured first, before any argument gets reassigned/resolved below, so
  ## it reflects exactly what the caller typed (standard R idiom, as in
  ## lm()/glm(): supports print()/update() on the returned object). See
  ## `options` in the return value for the resolved (post-default) values
  ## instead.
  mc <- match.call()

  time=rbind(start=proc.time())
  if (!is.null(seed)) set.seed(seed)


  ################### 1 CHECK INPUT #####################
  data <- data.table::as.data.table(data.table::copy(data))
  target <- if (is.null(target)) "all" else match.arg(target, c("all", "overlap"))

  ## `doubly.robust`'s default depends on `parametric`: AIPW's orthogonality
  ## protection (robust to misspecification of *either* nuisance) matters
  ## most when nuisances are flexibly/ML-estimated (parametric = FALSE, the
  ## overall default) -- so that case defaults to doubly.robust = TRUE. Under
  ## parametric = TRUE, the nuisance functional form is taken as given/
  ## trusted by construction, so AIPW is buying comparatively little (no
  ## protection against misspecification you've already assumed away) while
  ## still adding a second nuisance function's estimation noise on top of
  ## the singly-robust/FWL score -- so that case defaults to
  ## doubly.robust = FALSE instead, matching the classical parametric
  ## FWL/2SLS-style coefficient more directly. An explicit doubly.robust=
  ## TRUE/FALSE always overrides this regardless of `parametric`.
  doubly.robust <- if (is.null(doubly.robust)) !isTRUE(parametric) else doubly.robust

  ##Check formula and validate
  v <- validate_iv(fml, data)

  Y <- if (is.null(v$Y)) NULL else as.character(v$Y)[1L]
  D <- as.character(v$D)[1L]
  Z <- as.character(v$Z)[1L]
  X_forest <- v$X
  X=v$X

  FE <- v$FE

  X_expr_forest <- v$X_expr
  FE_expr <- v$FE_expr
  has_FE <- !is.null(FE_expr)

  ## Forest/search features are required.
  ## The main X part of fml is used for heterogeneous-effect / subset search.
  has_X_expr_forest <- !is.null(X_expr_forest) &&
    !identical(X_expr_forest, quote(1)) &&
    length(all.vars(X_expr_forest)) > 0L

  if (!has_X_expr_forest) {
    stop(
      "The main X part of `fml` is empty, so montest has no covariates for ",
      "forest/CART subset search.\n\n",
      "Please specify at least one variable before the first `|` in `fml`, e.g.\n",
      "  Y ~ X1 + X2 | fe1 + fe2 | D ~ Z\n\n",
      "If you only want fixed effects for identification, include any FE variable ",
      "also in the X part when it is meaningful as a moderator, e.g.\n",
      "  Y ~ year | group + year | D ~ Z\n",
      "or use `factor(year)` if you want year dummies as forest features.\n\n",
      "The FE part after `|` is absorbed for identification; it does not by itself ",
      "create forest/search features.",
      call. = FALSE
    )
  }

  X_expr_Z <- parse_one_sided_rhs(fml.Z, "fml.Z")
  X_expr_Q <- parse_one_sided_rhs(fml.Q, "fml.Q")

  if (is.null(X_expr_Z)) X_expr_Z <- X_expr_forest
  if (is.null(X_expr_Q)) X_expr_Q <- X_expr_forest

  has_X_expr_forest <- !is.null(X_expr_forest) && !identical(X_expr_forest, quote(1))
  has_X_expr_Z <- !is.null(X_expr_Z) && !identical(X_expr_Z, quote(1))
  has_X_expr_Q <- !is.null(X_expr_Q) && !identical(X_expr_Q, quote(1))


  if (is.null(X_expr_Z)) X_expr_Z <- X_expr_forest
  if (is.null(X_expr_Q)) X_expr_Q <- X_expr_forest

  ## `fml.varZ` controls the covariates used for the conditional-variance
  ## nuisance v(X) = Var(Z|X,FE) alone (see the `need_v_hat` block below),
  ## separately from `fml.Z`'s covariates for Z's own conditional mean --
  ## the two need not be the same set (e.g. a more parsimonious variance
  ## model for stability). Defaults to `fml.Z` (which itself already
  ## defaults to `fml`'s main X part), so leaving it unset reproduces
  ## exactly the previous behavior of reusing Z.hat's own covariates.
  X_expr_varZ <- parse_one_sided_rhs(fml.varZ, "fml.varZ")
  if (is.null(X_expr_varZ)) X_expr_varZ <- X_expr_Z

  ## Under parametric = TRUE, Z.hat/Q.hat's linear nuisance models are fit
  ## on the same rows used for testing (no cross-fitting), so the degrees-
  ## of-freedom correction must also account for X's rank, alongside the
  ## existing FE-rank correction (fe_rank_adj). Under parametric = FALSE,
  ## X is genuinely cross-fit, so no such penalty applies.
  x_rank_vars <- if (isTRUE(parametric)) {
    unique(c(all.vars(X_expr_Z), all.vars(X_expr_Q), all.vars(X_expr_varZ)))
  } else {
    character(0)
  }

  if (!is.null(inner.folds)) {
    if (!is.numeric(inner.folds) ||
        length(inner.folds) != 1L ||
        !is.finite(inner.folds) ||
        inner.folds != as.integer(inner.folds) ||
        inner.folds < 2L) {
      stop("inner.folds must be NULL or a single integer >= 2.", call. = FALSE)
    }
    inner.folds <- as.integer(inner.folds)
  }

  stopifnot(
    "linearD must be a single non-missing logical" =
      is.logical(linearD) && length(linearD) == 1L && !is.na(linearD),
    "linearZ must be a single non-missing logical" =
      is.logical(linearZ) && length(linearZ) == 1L && !is.na(linearZ),
    "doubly.robust must be a single non-missing logical" =
      is.logical(doubly.robust) && length(doubly.robust) == 1L && !is.na(doubly.robust)
  )

  ## Computed early (needs only doubly.robust/target, both already resolved)
  ## because it also gates the propensity clip/stabilization in the
  ## stabilize.scores block below, in addition to the v(X) block further
  ## down: whenever
  ## doubly.robust = FALSE & target == "overlap", binary-Z rows use the
  ## row-level empirical (Z-Z.hat)^2 instead of the closed-form e*(1-e),
  ## and that path needs an UNCLIPPED Z.hat to actually reproduce the
  ## classical FWL/OLS coefficient it's meant to match -- see
  ## make_scores_vec().
  need_ols_v <- !isTRUE(doubly.robust) && identical(target, "overlap")

  screen=match.arg(screen,c("stepdown","negative","nonpositive","minimum","none","fgk_relevant"))
  gridtypeY=match.arg(gridtypeY,c("equidistant","equisized"))
  gridtypeD=match.arg(gridtypeD,c("equidistant","equisized"))
  gridtypeZ=match.arg(gridtypeZ,c("equidistant","equisized"))
  if (is.null(crossfit)) {
    crossfit <- character()
  } else {
    crossfit <- match.arg(crossfit, c("Z", "Q", "C", "Y"), several.ok = TRUE)
  }
  stopifnot(shrink >= 0, shrink <= 1)
  testtype=match.arg(testtype,c("forest","CART"))
  if (testtype=="CART") shrink=0
  if ((is.null(cluster)==FALSE)&("CART" %in% testtype)) stop("Clustering not supported with testtype = CART.")
  if (testtype=="CART" && !is.null(gridpoints)) {
    warning(
      "`gridpoints` only applies to testtype = \"forest\"; it has no effect ",
      "with testtype = \"CART\" and will be ignored.",
      call. = FALSE
    )
  }

  if (!is.null(aipw.clip)) {
    stopifnot(
      is.numeric(aipw.clip),
      length(aipw.clip) == 1L,
      is.finite(aipw.clip),
      aipw.clip >= 0,
      aipw.clip < 1
    )
  }

  if (!is.numeric(minsize) || length(minsize) != 1L ||
      !is.finite(minsize) || minsize <= 0 || minsize != as.integer(minsize)) {
    stop("minsize must be a positive integer.")
  }

  minsize <- as.integer(minsize)
  if (is.null(condition)==TRUE) {
    if (is.null(Y)==TRUE) {
      condition="simple"
    } else condition="KR"
  }

  condition=match.arg(condition,c("simple","KR","MW","AHS","all"),several.ok=TRUE)
  if ("all" %in% condition) {
    condition=c("simple","KR","MW","AHS")
  }
  if (is.null(Y)==TRUE) {
    if (any(condition %in% c("KR","AHS","MW"))) {
      stop("Other conditions than (variants of) simple may not be used when Y is not specified. Specify the Y argument or use condition=simple")
    }
  }

  ## AHS's Q.hat fit (its i_yres rows, see the Q.hat-estimation block below)
  ## uses one extra regressor beyond X_expr_Q -- either Y itself or its FE/X-
  ## residualized version Y.res (y_name_rhs, built later from Y.res), which
  ## x_rank_vars otherwise has no way to know about. Folding Y itself into
  ## x_rank_vars adds exactly the +1 that extra regressor is worth (Y is
  ## virtually never collinear with the X's used for nuisance fitting, and
  ## its own rank contribution doesn't depend on whether the fitted model
  ## ends up using raw Y or Y.res), without needing to track the exact
  ## column name AHS's fit will end up using.
  if (isTRUE(parametric) && "AHS" %in% condition) {
    x_rank_vars <- unique(c(x_rank_vars, Y))
  }

  ## linearD/linearZ have no meaning for KR/MW -- their moment conditions are
  ## indicator-based, not slope-based, so they always use the margin/
  ## binarized representation regardless. Rather than silently ignoring an
  ## explicit linearD/linearZ = TRUE request for those conditions, fail
  ## loudly: an unmet explicit request should error, not be a no-op.
  if (isTRUE(linearD) && any(c("KR", "MW") %in% condition)) {
    stop(
      "linearD = TRUE has no effect for condition \"KR\"/\"MW\" (they always ",
      "use the binarized representation of D) and is disallowed in ",
      "combination with them. Drop \"KR\"/\"MW\" from `condition`, or set ",
      "linearD = FALSE.",
      call. = FALSE
    )
  }
  if (isTRUE(linearZ) && any(c("KR", "MW") %in% condition)) {
    stop(
      "linearZ = TRUE has no effect for condition \"KR\"/\"MW\" (they always ",
      "use the binarized representation of Z) and is disallowed in ",
      "combination with them. Drop \"KR\"/\"MW\" from `condition`, or set ",
      "linearZ = FALSE.",
      call. = FALSE
    )
  }

  check_integer_gt1 <- function(x, name) {
    if (!is.numeric(x) ||
        length(x) != 1L ||
        !is.finite(x) ||
        x != as.integer(x) ||
        x <= 1L) {
      stop(name, " must be a single integer larger than 1. Use e.g. 4L.", call. = FALSE)
    }

    as.integer(x)
  }

  Ysubsets <- check_integer_gt1(Ysubsets, "Ysubsets")
  Dsubsets <- check_integer_gt1(Dsubsets, "Dsubsets")
  Zsubsets <- check_integer_gt1(Zsubsets, "Zsubsets")

  ##VALIDATE POOL/SELECT
  if ((sum(pool=="none")==1)&(sum(pool=="all")==1)) stop("Do not specify both none and all in pool().")
  else if (sum(pool=="all")==1) pool=c("zmargin","dval","yval","condition","equation","sample")
  else if (sum(pool=="none")==1) pool=c()
  else if (is.null(pool)==FALSE) pool <- match.arg(
    pool,
    c("zmargin", "dval", "yval", "condition", "equation", "sample"),
    several.ok = TRUE
  )
  else pool=c("zmargin","dval","yval","sample")

  if ((sum(select=="none")==1)&(sum(select=="all")==1)) stop("Do not specify both none and all in select().")
  else if (sum(select=="all")==1) select=c("zmargin","dval","yval","condition","equation","sample")
  else if (sum(select=="none")==1) select=c()
  else if (is.null(select)==FALSE) select=match.arg(select,c("zmargin","dval","yval","condition","equation","sample"),several.ok=TRUE)
  else select="condition"

  if ("sample" %in% intersect(pool,select)) {
    stop("Sample cannot appear in both pool and select. Adaptive selection or pooling across sample halves might invalidate sample splitting.")
  }

  hat_vars <- grep("\\.hat$", names(data), value = TRUE)

  if (length(hat_vars)) {
    stop(
      "Variable names ending in '.hat' are reserved for internal use. ",
      "Please rename: ",
      paste(hat_vars, collapse = ", "),
      call. = FALSE
    )
  }
  if ("Q" %in% colnames(data)) stop("Data contains variable named Q, which is reserved for internal use. Please rename.")

  if ((length(D)!=1)|(!(D %in% colnames(data)))) {
    stop("Argument D must be the name of a single column in data")
  }

  if (is.null(weight)==FALSE) {
    wvar=weight
    if ((length(weight)!=1)|(!(weight %in% colnames(data)))) {
      stop("Argument weight must be the name of a single column in data")
    }
  } else wvar=NA_character_

  if (is.null(cluster)==FALSE) {
    clvar=cluster
    if ((length(cluster)!=1)|(!(cluster %in% colnames(data)))) {
      stop("Argument cluster must be the name of a single column in data")
    }
  } else clvar=NA_character_

  if ((length(Z)!=1)|(!(Z %in% colnames(data)))) {
    stop("Argument Z must be the name of a single column in data")
  }

  if (length(X)!=sum(X %in% colnames(data))) {
    stop("All variables in argument X must exist in data")
  }

  if (is.null(Y)==FALSE) {
    if (sum(!(Y %in% colnames(data)))>0) {
      stop("All entries in argument Y must be the name of a column in data")
    }
  }

  if (!any(condition %in% c("AHS", "MW", "KR")) && !is.null(Y)) {
    Y <- NULL
  }

  ## Downgrade linearD/linearZ if D or Z is binary in the data.
  ## Number of support points in original D and Z
  Zname=Z
  Dname=D

  J <- data.table::uniqueN(data[[D]])
  K <- data.table::uniqueN(data[[Z]])

  ## J/K get reassigned below to post-binning bin counts (length(Dsup)/
  ## length(Zsup)) once binarize_var() runs -- keep the TRUE support-point
  ## counts around under different names for checks that must reflect the
  ## original data, not how coarsely it happens to have been binned (e.g.
  ## the MW binary-instrument/treatment guards further down).
  J_true <- J
  K_true <- K

  if (isTRUE(linearZ) && K <= 2L) {
    warning(
      "linearZ = TRUE requested, but Z has only ", K, " support point(s) in ",
      "the data; using linearZ = FALSE (a binary/near-constant instrument ",
      "has no meaningful linear scoring).",
      call. = FALSE
    )
    linearZ <- FALSE
  }
  if (isTRUE(linearD) && J <= 2L) {
    warning(
      "linearD = TRUE requested, but D has only ", J, " support point(s) in ",
      "the data; using linearD = FALSE (a binary/near-constant treatment ",
      "has no meaningful linear scoring).",
      call. = FALSE
    )
    linearD <- FALSE
  }

  ## Conditions for which linearD/linearZ have an effect. KR/MW always use
  ## the margin/binarized representation of D and Z regardless of these
  ## flags -- their moment conditions are indicator-based, not slope-based.
  linear_conditions <- c("simple", "AHS")

  has_linear_conditions <- any(condition %in% linear_conditions)

  ## Conditions that always need margin/binary versions
  nonlinear_conditions <- setdiff(condition, linear_conditions)

  ## Need binarized Z if:
  ## - any non-linear-aware condition (KR/MW) is requested, or
  ## - simple/AHS are requested with linearZ = FALSE
  need_binarized_Z <-
    length(nonlinear_conditions) > 0L ||
    (has_linear_conditions && !linearZ)

  ## Need original/linear Z if simple/AHS are requested with linearZ = TRUE
  need_linear_Z <- has_linear_conditions && linearZ


  ## Need margin conditions if:
  ## - simple/AHS are requested with linearZ = FALSE (i.e. margin-Z rows exist)
  ## - KR/MW are requested (always margin-Z)
  ##
  ## These rows require the one-sided monotonicity pre-check (run on
  ## margin-Z rows only, i.e. data[z_is_linear == FALSE]), which itself
  ## needs the binarized treatment margin variable D.bin.
  has_margin_conditions <-
    (has_linear_conditions && !linearZ) ||
    any(c("KR", "MW") %in% condition)

  ## Need binarized D for Q construction if:
  ## - any non-linear-aware condition (KR/MW) is requested, or
  ## - simple/AHS are requested with linearD = FALSE
  need_binarized_D_for_Q <-
    length(nonlinear_conditions) > 0L ||
    (has_linear_conditions && !linearD)

  ## Need binarized D for the one-sided monotonicity screen.
  need_binarized_D_for_os <- has_margin_conditions

  ## Actually create D.bin if either Q construction or the one-sided screen needs it.
  need_binarized_D <- need_binarized_D_for_Q || need_binarized_D_for_os

  ## Need original/linear D if simple/AHS are requested with linearD = TRUE
  need_linear_D <- has_linear_conditions && linearD

  ##Validate select/pool choices for CART
  if (identical(testtype, "CART")) {
    if (is.null(pool)) pool <- character()
    if (is.null(select)) select <- character()

    pool <- unique(as.character(pool))
    select <- unique(as.character(select))

    if (length(pool) > 0L) {
      stop(
        "Pooling margins/sample incompatible with testtype = \"CART\". "
      )
    }

    allowed_select <- c("zmargin","dval","yval","condition","equation","sample")
    bad_select <- setdiff(select, allowed_select)
    if (length(bad_select) > 0L) {
      stop(
        "Invalid `select` for testtype = \"CART\". ",
        "For CART, `select` may only contain \"zmargin\",\"dval\",\"yval\",\"condition\",\"equation\",\"sample\". ",
        "Invalid entries: ",
        paste(bad_select, collapse = ", "),
        call. = FALSE
      )
    }

    overlap <- intersect(pool, select)
    if (length(overlap) > 0L) {
      stop(
        "Invalid `pool`/`select` combination for testtype = \"CART\". ",
        "CART does not allow overlap between `pool` and `select`. ",
        "Overlapping entries: ",
        paste(overlap, collapse = ", "),
        call. = FALSE
      )
    }

  }

  time <- rbind(time, "Check input" = proc.time())

  ###################### 2 Prepare data #########################3


  vars_forest <- all.vars(X_expr_forest)
  vars_Z      <- all.vars(X_expr_Z)
  vars_varZ   <- all.vars(X_expr_varZ)
  vars_Q      <- all.vars(X_expr_Q)
  vars_FE     <- if (has_FE) all.vars(FE_expr) else character()

  allvars <- unique(c(
    vars_forest,
    vars_Z,
    vars_varZ,
    vars_Q,
    vars_FE,
    Y, D, Z,
    weight,
    cluster
  ))
  allvars <- allvars[!is.na(allvars)]

  data <- data[, ..allvars]
  dropped <- sum(!complete.cases(data))

  if (dropped > 0L) {
    message("Note: dropped ", dropped,
            " observations with missing data on one or more input variables.")
  }

  data <- data[complete.cases(data)]

  if (!is.null(weight)) {
    bad_w <- data[[weight]] <= 0
    n_bad_w <- sum(bad_w)
    if (n_bad_w > 0L) {
      message(
        "Note: dropped ", n_bad_w,
        " observation(s) with zero or negative `weight`. A zero weight leaves ",
        "that row's nuisance fit unpopulated downstream; keeping it around would ",
        "either silently drop it again later or (once a sandwich weight is added ",
        "on top of a zero base weight) corrupt aggregate statistics instead of ",
        "being simply excluded."
      )
      data <- data[!bad_w]
    }
  }

  ## Two distinct FE-driven exclusions, run in this order (see drop_singletons/
  ## drop_novar_Z docs): (1) true singletons (FE levels with a single
  ## observation) are regressor-agnostic and affect degrees of freedom for
  ## every nuisance, not just Z's; (2) no-within-FE-variation-in-Z rows are
  ## specific to identifying gamma and to the numerical stability of the
  ## Z-variance nuisance. A true singleton is trivially also a no-variation-
  ## in-Z row, so running (1) first just avoids (2) rediscovering the same
  ## rows -- it is not required for correctness either way.
  n_dropped_singletons <- 0L
  n_dropped_novar_Z <- 0L

  if (has_FE && isTRUE(drop_singletons)) {
    idx_singleton <- drop_fe_singletons(data, FE_expr, placeholder_y = D)
    n_dropped_singletons <- length(idx_singleton)

    if (n_dropped_singletons > 0L) {
      message(
        "Note: dropped ", n_dropped_singletons,
        " observation(s) belonging to a fixed-effect singleton (an FE level ",
        "with only one observation, possibly only after removing other ",
        "singletons in a different FE dimension). These are perfectly ",
        "explained by their own FE dummy and contribute nothing to any ",
        "regressor's FE-adjusted coefficient or its standard error -- see ",
        "`drop_singletons`. Set `drop_singletons = FALSE` to keep them."
      )
      data <- data[-idx_singleton]
    }
  }

  if (has_FE && isTRUE(drop_novar_Z)) {
    idx_novar_Z <- drop_novar_Z_rows(data, Z, FE_expr, weight = weight)
    n_dropped_novar_Z <- length(idx_novar_Z)

    if (n_dropped_novar_Z > 0L) {
      message(
        "Note: dropped ", n_dropped_novar_Z,
        " observation(s) in fixed-effect cells with no within-cell variation ",
        "in `", Z, "`. Their residualized `", Z, "` is exactly 0, so they ",
        "contribute nothing to identifying the coefficient on `", Z, "` (or ",
        "its variance) in the classical FWL sense, and are the primary ",
        "source of the \"conditional variance v(X) for continuous-Z rows\" ",
        "clip warnings this option exists to avoid -- see `drop_novar_Z`. ",
        "Set `drop_novar_Z = FALSE` to keep them."
      )
      data <- data[-idx_novar_Z]
    }
  }

  n=nrow(data)
  if (is.null(cluster)==FALSE) {
    G <- data.table::uniqueN(data[[cluster]])
    if (G<=2*minsize) stop("Number of clusters is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
  } else {
    if (n<=2*minsize) stop("Number of observations is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
    G=NA
    }
  obs=c(N=n,G=G)
  data[, id_ := .I]


  ##OUTER SPLIT
  stopifnot(is.logical(stratify), length(stratify) == 1L, !is.na(stratify))
  stratify_explicit <- !missing(stratify)

  strat <- NULL

  if (isTRUE(stratify) && is.null(cluster)) {
    ## Stratifying a 2-fold split by Z's raw value only makes sense for a
    ## small number of support points -- make_group_folds() groups by exact
    ## value (by = Z), and its fallback for a too-small stratum
    ## (fold[ok] <- 1L, since a stratum of size 1 can't be split into 2
    ## non-empty folds) fires almost everywhere for a genuinely continuous/
    ## high-cardinality Z, collapsing the whole sample split into one
    ## degenerate all-in-one-fold assignment -- which then crashes
    ## fit_models() downstream once the opposite sample turns out empty.
    ## K_true <= 2 mirrors the same binary/continuous threshold used for the
    ## has_FE/linearZ representation choice elsewhere in this file.
    if (K_true <= 2L) {
      strat <- Z
    } else if (isTRUE(stratify_explicit)) {
      warning(
        "stratify = TRUE was requested, but Z has ", K_true, " support ",
        "point(s) in the data (not binary). Stratifying the sample split ",
        "by a continuous/high-cardinality Z can collapse it into a ",
        "degenerate, all-in-one-fold split, so an unstratified split is ",
        "used instead.",
        call. = FALSE
      )
    }
  }
  make_group_folds(data,K = 2,cluster_name = cluster, fold_col = "sample",verbose = FALSE,diag_prefix=NULL,strat_col=strat)

  ##OPTIONAL INNER SPLIT
  if (is.null(inner.folds)==FALSE) {
    make_group_folds(data,K = inner.folds,cluster_name = cluster,fold_col = "cf_fold",verbose = FALSE,by_col="sample",diag_prefix=NULL,strat_col=strat)
    foldname="cf_fold"
  } else {
    foldname=NULL
  }

  if (is.null(cluster)==TRUE) cluster="id_"

  time=rbind(time,"Prepare data"=proc.time())

  ############### 3 Discretize Z, D and Y into subsets ###############

  if (need_binarized_D) {
    ##bin treatment
    data <- binarize_var(
      data    = data,
      var     = D,
      ngroups = Dsubsets,
      gridtype = gridtypeD,
      wvar    = wvar,
      newvar=paste0(D,".bin")
    )
    Dbincol=paste0(D,".bin")
    Dsup=sort(unique(data[,get(Dbincol)]));J=length(Dsup)
  } else {
    J=Inf;Dbincol=NULL
  }

  if (need_binarized_Z){
    ##bin instrument
    data <- binarize_var(
      data    = data,
      var     = Z,
      ngroups = Zsubsets,
      gridtype = gridtypeZ,
      wvar    = wvar,
      newvar=paste0(Z,".bin")
    )
    Zbincol=paste0(Z,".bin")
    Zsup=sort(unique(data[,get(Zbincol)]));K=length(Zsup)
  } else {
    K=Inf;Zbincol=NULL
  }

  if ("KR" %in% condition) { ##bin outcome(s)
    data <- binarize_var(
      data    = data,
      var     = Y,
      ngroups = Ysubsets,
      gridtype = gridtypeY,
      wvar    = wvar,
      newvar = paste0(Y, ".bin")
    )
    Ybincol=paste0(Y,".bin")
    Ysup=sort(unique(data[,get(Ybincol)]));L=length(Ysup)
  }

  n=nrow(data)

  ## Use the TRUE support-point counts (J_true/K_true, captured before
  ## binarize_var() overwrote J/K with post-binning bin counts) -- checking
  ## the post-binning counts here would let e.g. Zsubsets=2 mask a genuinely
  ## multivalued Z and defeat this guard entirely.
  if (J_true>2&("MW" %in% condition)) stop("Multivalued treatment not supported with condition MW.")
  if (K_true>2&("MW" %in% condition)) {
    stop(
      "condition = \"MW\" requires a genuinely binary instrument Z (K <= 2 ",
      "support points in the data); Z has ", K_true, " here. Mourifie and Wan's ",
      "(2017) conditions are derived for a true binary instrument -- testing ",
      "them against an arbitrary margin/threshold cut of a multivalued Z ",
      "would not correspond to their theory. Use \"simple\", \"KR\", or ",
      "\"AHS\" for a multivalued instrument.",
      call. = FALSE
    )
  }
  if (has_FE&("MW" %in% condition)) {
    stop(
      "condition = \"MW\" does not support fixed effects in `fml`. MW's Q ",
      "statistic uses Z.hat directly as a propensity weight with no AIPW/FWL ",
      "orthogonalization; with FE present, Z.hat is unreliable in exactly ",
      "the way that machinery exists to protect against for the other ",
      "conditions. Drop the FE term from `fml`, or use \"simple\", \"KR\", ",
      "or \"AHS\" instead.",
      call. = FALSE
    )
  }
  if (J==2&K==2&is.null(X)==TRUE&!any(condition %in% c("AHS","MW","KR"))) {
    stop("Nothing to test with a binary treatment, a binary instrument, the simple first stage condition and no variables in X.")
  }

  time=rbind(time,"Binarize treatment, instrument and outcome"=proc.time())

  ######################## 6a STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
  margins=c()

  ## Z's representation is a single global choice for the whole call, not a
  ## per-row one: need_binarized_Z == !linearZ and need_linear_Z == linearZ
  ## always hold here (linearZ = TRUE is disallowed together with "KR"/"MW"
  ## earlier in this function, which is exactly what would otherwise let
  ## the two representations coexist in one stacked dataset). So there are
  ## only three cases -- linear, margin-stacked, or plain binary -- never a
  ## mix, and z_is_linear ends up the same value for every row.
  if (linearZ) {

    ## No margin stacking: Z used on its raw/linear scale throughout.
    data[, (Z) := get(Zname)]
    data[, zmargin := NA_real_]
    data[, z_is_linear := TRUE]

    if (!is.null(Zbincol) && Zbincol %in% names(data)) {
      data[, (Zbincol) := NULL]
    }

  } else if (K > 2L) {

    ## Margin-specific expansion: one row per (original observation,
    ## candidate zmargin threshold).
    ##
    ## - lowest Z gets second-lowest margin, so Zbin >= zmargin gives 0
    ## - highest Z gets highest margin, so Zbin >= zmargin gives 1
    ## - intermediate Z gets current and next-higher margins, giving one 1 and one 0
    stopifnot(!is.null(Zbincol), Zbincol %in% names(data))

    zvals <- sort(unique(data[[Zbincol]]))
    zvals <- zvals[!is.na(zvals)]

    if (length(zvals) < 2L) {
      stop("Need at least two non-missing Z values to construct zmargin.")
    }

    minZ <- zvals[1L]
    maxZ <- zvals[length(zvals)]

    zmap <- data[, {
      z0 <- get(Zbincol)
      k0 <- match(z0, zvals)

      if (is.na(z0) || is.na(k0)) {
        zm <- NA_real_
      } else if (z0 == minZ) {
        zm <- zvals[2L]
      } else if (z0 == maxZ) {
        zm <- zvals[length(zvals)]
      } else {
        zm <- c(z0, zvals[k0 + 1L])
      }

      .(rowid = .I, zmargin = zm)
    }, by = id_]

    data <- data[zmap$rowid]
    data[, zmargin := zmap$zmargin]
    data[, z_is_linear := FALSE]

    margins <- unique(c(margins, "zmargin"))

    data[, (Z) := as.integer(get(Zbincol) >= zmargin)]

    if (!is.null(Zbincol) && Zbincol %in% names(data)) {
      data[, (Zbincol) := NULL]
    }

  } else {

    ## Binary-Z case (K <= 2): a single cutpoint, no margin stacking needed.
    if (!is.null(Zbincol) && Zbincol %in% names(data)) {
      data[, (Z) := get(Zbincol)]
      data[, (Zbincol) := NULL]
    }

    data[, zmargin := NA_real_]
    data[, z_is_linear := FALSE]
  }

  data[, z_is_linear_raw := z_is_linear]
  data[, z_use_linear_score := z_is_linear]

  ## Continuous-Z rows (FWL, or AIPW with a fitted variance nuisance) are
  ## chosen whenever:
  ##   - has_FE: under FE, Z.hat is never a directly-fitted probability --
  ##     estimate_conditional_mean() always builds it as an FE mean
  ##     (Zbar_g) plus an ML fit of the FE-residual on Xtilde, recombined
  ##     additively. That recombination is unconstrained: nothing about it
  ##     guarantees Z.hat in [0,1], for a multivalued Z OR a genuinely
  ##     binary one. For binary Z it is *true* that Var(Z|X,FE) =
  ##     e(X,FE)*(1-e(X,FE)) exactly, as an identity about the estimand --
  ##     but that identity is no protection when the ESTIMATE of e(X,FE)
  ##     itself can land outside [0,1] (small/unbalanced FE cells are
  ##     exactly where this happens, and exactly where you have the least
  ##     evidence to fall back on). Using the closed form there means the
  ##     score construction needs a bound the nuisance estimator cannot
  ##     promise, and fails hard (validate_and_clip()'s hard_lower/
  ##     hard_upper) rather than merely losing precision.
  ##
  ##     The fitted-v(X) alternative targets the exact same Var(Z|X,FE) --
  ##     (Z-Z.hat)^2 has conditional expectation e(X,FE)*(1-e(X,FE)) under
  ##     correct specification, so this is not a different estimand, just a
  ##     different estimator of it -- but it only ever needs v(X) floored
  ##     near 0 (small negative excursions near a genuine near-zero-variance
  ##     boundary are treated as estimation noise, not misspecification --
  ##     see make_scores_vec()'s idx_linear branch), and Z.hat itself is
  ##     then used solely additively in (Z-Z.hat), never as a multiplicative/
  ##     bounded quantity. A one-sided, wide-slack requirement instead of a
  ##     two-sided one that typically sits close to a boundary in exactly
  ##     the low-overlap cells FE identification relies on.
  ##
  ##     Without FE, this tradeoff reverses: e(X) is fit directly (no FE
  ##     recombination), so it already has whatever bound its own model
  ##     provides (a closed-form logit/forest classifier stays in [0,1] by
  ##     construction), and e(X)*(1-e(X)) is then genuinely the cheaper,
  ##     lower-noise choice with no offsetting robustness benefit. So the
  ##     binary/closed-form path remains the default whenever !has_FE.
  ##   - linearZ = TRUE: the user explicitly asked for Z on its raw/linear
  ##     scale (see `linearZ` docs); already implies K_true > 2, since
  ##     linearZ is downgraded to FALSE above whenever Z is binary.
  ## `parametric` is deliberately NOT a trigger here: nuisance fragility
  ## (e.g. an LPM-based Z.hat falling outside [0,1]) is instead caught by
  ## `validate_and_clip()` at score-construction time, fully decoupling the
  ## parametric/semiparametric choice from the binary/continuous
  ## representation choice. `doubly.robust` (AIPW vs. FWL) is likewise an
  ## independent choice -- see the `fit_models()` calls below.
  if (has_FE || isTRUE(linearZ)) {
    data[, z_use_linear_score := TRUE]
  }

  ## ----------------------
  ## CREATE X tilde/forest matrix
  ## -----------------------

  X_forest <- NULL
  X_forest_info <- NULL

  if (has_X_expr_forest) {

    X_forest_info <- make_X_residualized_from_FE(
      DT          = data,
      x_expr      = X_expr_forest,
      fe_expr     = FE_expr,
      prefix      = "__xf",
      by          = margins,
      weight      = weight,
      fixest_opts = Rparameters
    )

    X_forest <- X_forest_info$x_names

    if (is.null(X_forest) || length(X_forest) == 0L) {
      stop(
        "Internal error: the main X formula is non-empty, but no forest feature ",
        "columns were created. This likely means `make_X_residualized_from_FE()` ",
        "does not handle the no-FE case correctly.",
        call. = FALSE
      )
    }
  }

  ## Under parametric = FALSE, estimate_conditional_mean() re-derives its own
  ## FE-residualized X (one feols() fit per covariate per margin group) via
  ## make_X_residualized_from_FE() whenever x_names is NULL. When Z.hat's/
  ## Q.hat's covariate formula is literally the same as the main formula's
  ## (fml.Z/fml.Q left at their default), that work exactly duplicates what
  ## X_forest_info just computed above -- reuse its columns instead of
  ## refitting. Gated on `fixest_opts` also matching: it doesn't affect the
  ## residualized *values* (feols_partial_out() only ever extracts residuals/
  ## fitted, never SEs), but there's no need to assume that here when it can
  ## just be checked.
  X_names_reuse_for_Z <- if (
    !isTRUE(parametric) &&
      identical(X_expr_Z, X_expr_forest) &&
      identical(Rparameters, Zparameters)
  ) X_forest_info$x_names else NULL

  X_names_reuse_for_Q <- if (
    !isTRUE(parametric) &&
      identical(X_expr_Q, X_expr_forest) &&
      identical(Rparameters, Qparameters)
  ) X_forest_info$x_names else NULL

  ## ---------------------------------------------------------------------
  ## Estimate Z.hat for each Z margin in stacked data
  ## ---------------------------------------------------------------------

  zhat <- paste0(Z, ".hat")

  Z_hat_info <- estimate_conditional_mean(
    DT = data,
    y_name = Z,
    x_expr = X_expr_Z,
    fe_expr = FE_expr,
    out_hat = zhat,
    by = margins,
    sample_var = "sample",
    weight = weight,
    cluster = cluster,
    parametric = parametric,
    foldname = foldname,
    crossfit = crossfit,
    crossfit_label = "Z",
    forest_opts = Zparameters,
    fixest_opts = Zparameters,
    x_names = X_names_reuse_for_Z,
    x_prefix = "__xz",
    keep_x = TRUE,
    return_residual = FALSE,
    partial_out_y_fe = TRUE
  )




  ## ---------------------------------------------------------------------
  ## Design-side finite-sample correction of Z.hat -- see stabilize_zhat()'s
  ## own docs (functions.R) for the full rationale (mean-shift for the
  ## continuous representation, Hajek/self-normalized IPW rescaling for the
  ## binary one) and the need_ols_v exclusion/warning. Factored out of
  ## montest() so other callers of fit_models()/make_scores_vec() (e.g.
  ## seqtest) can drive the same correction without duplicating it.
  ## ---------------------------------------------------------------------

  stabilize_zhat(
    data = data,
    Z = Z,
    margins = margins,
    parametric = parametric,
    weight = weight,
    aipw.clip = aipw.clip,
    has_FE = has_FE,
    need_ols_v = need_ols_v,
    stabilize.scores = stabilize.scores
  )

  ## ---------------------------------------------------------------------
  ## Conditional-variance nuisance v(X) = Var(Z|X). Continuous-Z rows always
  ## need this from an external source (there is no closed form). Binary-Z
  ## rows have a closed form (v = Z.hat*(1-Z.hat)) available for free, but
  ## whether that closed form is actually USED depends on the same
  ## doubly.robust/target logic as the continuous case -- see below.
  ##
  ## A genuine per-observation fit (`need_v_hat`, continuous rows only) is
  ## needed when doubly.robust = TRUE (the AIPW orthogonality argument
  ## requires a genuinely local v(X)) or target == "all" (a constant v
  ## would silently collapse "all" and "overlap" into the same estimand,
  ## since dividing every row by the same number doesn't change relative
  ## weighting). Binary rows never need a fitted model for this case: the
  ## closed form e*(1-e) already varies by X for free.
  ##
  ## Otherwise (doubly.robust = FALSE & target == "overlap", `need_ols_v`)
  ## every row -- binary and continuous alike -- gets its OWN raw
  ## (Z_i - Z.hat_i)^2 as v_i, not a value averaged/pooled across a group.
  ## This is exactly what classical FWL/OLS uses (its coefficient is a sum
  ## divided by a sum, never a per-row division by a pre-averaged variance),
  ## so `crv1_mean()`'s theta = sum(w*score)/sum(w) reproduces the FWL/OLS
  ## coefficient at whatever level it is later summed over: the full
  ## sample, one margins/sample block, or a leaf found afterwards by
  ## forest_test()/CART_test()'s own subgroup search. A pooled/broadcast
  ## constant would only match FWL/OLS at the specific group it was pooled
  ## over -- any finer subgroup nested inside that group (which is exactly
  ## what the search is built to find) would then be normalized by a
  ## variance that isn't its own, contaminating the effect-heterogeneity
  ## test with unrelated heterogeneity in Var(Z|X). This also means the
  ## row-level value never requires Z.hat to represent a valid probability,
  ## so it is robust to Z.hat straying outside [0,1] (e.g. an LPM fit on a
  ## rare/near-degenerate binary Z) in exactly the situation where the
  ## classical FWL/OLS coefficient this is meant to reproduce would also be
  ## robust to it -- see `make_scores_vec()`.
  ## ---------------------------------------------------------------------

  zvarhat <- paste0(Z, ".var.hat")
  z_use_lin_any <- any(data$z_use_linear_score, na.rm = TRUE)
  need_v_hat <- z_use_lin_any &&
    (isTRUE(doubly.robust) || identical(target, "all"))
  ## need_ols_v was already computed early (before the stabilize.scores
  ## block), since it also gates the propensity clip/stabilization there.

  if (need_v_hat || need_ols_v) {
    ## `data[[Z]]`/`data[[zhat]]` (plain `[[`, evaluated here, not inside a
    ## `data[i, j]` call) sidestep data.table's column-name masking of `j` --
    ## masking would otherwise silently resolve the bare symbol `Z` to the
    ## *column* named "Z" instead of the R variable holding that name,
    ## whenever the instrument column happens to literally be called "Z".
    z_resid_sq <- paste0(Z, "__resid_sq__")
    z_vals <- data[[Z]]
    zhat_vals <- data[[zhat]]
    data[, (z_resid_sq) := (z_vals - zhat_vals)^2]
  }

  if (need_v_hat) {
    ## `X_expr_varZ` (from `fml.varZ`) may name a different covariate set
    ## than `X_expr_Z`. Reuse Z_hat_info's already-residualized columns only
    ## when the two formulas actually agree (the default, `fml.varZ` unset)
    ## -- otherwise let estimate_conditional_mean() residualize `X_expr_varZ`
    ## itself fresh (x_names = NULL triggers this internally).
    X_names_reuse_for_varZ <- if (
      !isTRUE(parametric) && identical(X_expr_varZ, X_expr_Z)
    ) Z_hat_info$x_names else NULL

    ## Restricted to continuous-representation rows: binary rows never need
    ## a fitted model here (see comment block above).
    estimate_conditional_mean(
      DT = data,
      y_name = z_resid_sq,
      x_expr = X_expr_varZ,
      fe_expr = FE_expr,
      out_hat = zvarhat,
      by = margins,
      sample_var = "sample",
      weight = weight,
      cluster = cluster,
      parametric = parametric,
      foldname = foldname,
      crossfit = crossfit,
      crossfit_label = "Z",
      forest_opts = Zparameters,
      fixest_opts = Zparameters,
      x_names = X_names_reuse_for_varZ,
      x_prefix = "__xzv",
      keep_x = FALSE,
      return_residual = FALSE,
      partial_out_y_fe = TRUE,
      i = which(data$z_use_linear_score)
    )
  }

  if (need_ols_v) {
    ## Every row -- binary and continuous alike -- gets its own raw
    ## (Z_i-Z.hat_i)^2, no group averaging. See the comment block above for
    ## why this (rather than a pooled/broadcast constant) is what makes the
    ## resulting coefficient reproduce classical FWL/OLS at any level of
    ## aggregation, not just the level it would have been pooled at.
    data[, (zvarhat) := get(z_resid_sq)]
  }

  if (need_v_hat || need_ols_v) {
    data[, (z_resid_sq) := NULL]
  }

  ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
  if (any(condition %in% c("MW", "AHS")) && isTRUE(Y.res)) {

    y_name_rhs <- paste0(Y, ".res")
    yhat <- paste0(Y, ".hat")

    Y_res_info <- estimate_conditional_mean(
      DT = data,
      y_name = Y,
      x_expr = X_expr_Q,
      fe_expr = FE_expr,
      out_hat = yhat,
      out_resid = y_name_rhs,
      by = margins,
      sample_var = "sample",
      weight = weight,
      cluster = cluster,
      parametric = parametric,
      foldname = foldname,
      crossfit = crossfit,
      crossfit_label = "Y",
      forest_opts = Yparameters,
      fixest_opts = Yparameters,
      x_names = NULL,
      x_prefix = "__xy",
      keep_x = FALSE,
      return_residual = TRUE,
      partial_out_y_fe = TRUE
    )

    ## Optional cleanup: keep Y.res, drop nuisance fitted value.
    data[, (yhat) := NULL]

  } else {
    y_name_rhs <- Y
  }

  ##Check for onesided noncompliance
  ## Skipped entirely when has_FE = TRUE: this is a pooled/unconditional
  ## check (no FE, no X, no fitted quantities at all -- see
  ## test_one_sided_noncompliance()), so it can't see FE-localized
  ## noncompliance patterns and can misjudge a margin as trivially
  ## satisfied purely from cross-FE-group composition effects (a
  ## Simpson's-paradox-style aggregation bias), even though the surviving
  ## margins get scored by FE-aware machinery regardless. Below, both
  ## make_linear_blocks() and the KR block fall back to testing every
  ## margin directly (built from Zsup/Dsup) instead of filtering by os_res.
  if (has_margin_conditions && !has_FE) {
    if (is.null(Dbincol) || !(Dbincol %in% names(data))) {
      stop(
        "Internal error: one-sided margin conditions require `", paste0(D, ".bin"),
        "`, but binarized D was not created.",
        call. = FALSE
      )
    }

    os_res <- test_one_sided_noncompliance(
      data = data[z_is_linear == FALSE],
      D = Dbincol,
      Z = Z,
      zmargin_var = margins
    )

    osmargins <- margins

    ## Only KR is legitimately exempt here: it builds its margins from a
    ## different, finer table (os_res$exact) that can still have non-trivial
    ## rows even when os_res$threshold is entirely one-sided. MW builds its
    ## margins directly from this same os_res$threshold table (see the MW
    ## block below), so if this stop is skipped for MW and every threshold
    ## cell is one-sided, MW's own idx_blocks ends up empty and the run
    ## fails later with a confusing, unrelated error instead of this one.
    if (
      nrow(os_res$threshold[one_sided == FALSE]) == 0 &&
      !("KR" %in% condition)
    ) {
      stop(
        "One-sided noncompliance for all margins of Z and D - ",
        "all specified conditions are trivially satisfied."
      )
    }
  }

  ### ---------------------------------- #####
  ## Check for residual variation in Z
  ### ---------------------------------- #####

  byvars <- unique(c("sample", margins))
  byvars <- byvars[byvars %in% names(data)]

  z_col_work <- Z
  zhat_col <- zhat

  stopifnot(z_col_work %in% names(data))
  stopifnot(zhat_col %in% names(data))

  ## Common summaries
  data[, n_group := .N, by = byvars]

  data[
    ,
    nZ := data.table::uniqueN(.SD[[z_col_work]], na.rm = TRUE),
    by = byvars,
    .SDcols = z_col_work
  ]

  data[
    ,
    sdZ := stats::sd(.SD[[z_col_work]], na.rm = TRUE),
    by = byvars,
    .SDcols = z_col_work
  ]

  data[
    ,
    sd_res := stats::sd(.SD[[z_col_work]] - .SD[[zhat_col]], na.rm = TRUE),
    by = byvars,
    .SDcols = c(z_col_work, zhat_col)
  ]

  ## Initialize helper columns
  if (!("min_cell_Z" %in% names(data))) {
    data[, min_cell_Z := NA_integer_]
  }

  data[, bad_Z := FALSE]

  ## Continuous / linear raw-Z rows.
  ##
  ## Use raw z_is_linear here, not z_use_linear_score.
  ## z_use_linear_score is for score construction. This diagnostic is about
  ## whether the working Z has usable variation in the relevant sample ?? margin cell.
  if ("z_is_linear" %in% names(data)) {
    data[
      z_is_linear == TRUE,
      bad_Z := (
        n_group < min_n |
          is.na(sdZ) | sdZ == 0 |
          is.na(sd_res) | sd_res == 0
      ),
      by = byvars
    ]
  }

  ## Discrete / margin-Z rows.
  if (has_margin_conditions) {
    data[
      z_is_linear == FALSE,
      min_cell_Z := {
        ztab <- table(.SD[[z_col_work]], useNA = "no")
        if (length(ztab) == 0L) NA_integer_ else min(ztab)
      },
      by = byvars,
      .SDcols = z_col_work
    ]

    data[
      z_is_linear == FALSE,
      bad_Z := (
        n_group < min_n |
          nZ < 2 |
          is.na(min_cell_Z) | min_cell_Z < min_n |
          is.na(sd_res) | sd_res == 0
      ),
      by = byvars
    ]
  }

  ## If z_is_linear is absent for some reason, fall back to a generic residual-variation check.
  if (!("z_is_linear" %in% names(data))) {
    data[
      ,
      bad_Z := (
        n_group < min_n |
          is.na(sdZ) | sdZ == 0 |
          is.na(sd_res) | sd_res == 0
      ),
      by = byvars
    ]
  }

  if (length(margins) == 0L) {
    ## No margins: if any sample part fails, drop everything.
    drop_all <- data[, any(bad_Z, na.rm = TRUE)]

    if (isTRUE(drop_all)) {
      stop("At least one sample part has no residual variation in Z.")
    }

  } else {
    ## Identify margin cells where at least one sample part fails.
    bad_margins <- unique(
      data[bad_Z == TRUE, ..margins]
    )

    if (nrow(bad_margins) > 0L) {
      if (length(margins) == 1L) {
        msg <- paste(bad_margins[[margins]], collapse = ", ")
        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part has no usable variation in Z: ",
          msg
        )
      } else {
        bad_labels <- apply(bad_margins, 1L, function(r) {
          paste(paste(names(r), r, sep = "="), collapse = ", ")
        })

        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part has no usable variation in Z:\n",
          paste("  -", bad_labels, collapse = "\n")
        )
      }

      ## Anti-join on margins only: removes both sample parts.
      bad_margins[, drop_bad_Z__ := TRUE]

      data <- bad_margins[data, on = margins]
      data <- data[is.na(drop_bad_Z__)][, drop_bad_Z__ := NULL]
    }
  }

  helper_cols_Z <- intersect(
    c("n_group", "nZ", "sdZ", "sd_res", "min_cell_Z", "bad_Z"),
    names(data)
  )

  if (length(helper_cols_Z) > 0L) {
    data[, (helper_cols_Z) := NULL]
  }

  time=rbind(time,"Stack data for Z margins and estimate nuisance for Z"=proc.time())

  #######STACK ACROSS MARGINS ##########


  # --------------------------------------------------
  # Build one unified condition index
  #   Generic columns:
  #     condition, linear, equation, zmargin, dval, yval, Avals
  # --------------------------------------------------

  idx_blocks <- list()

  ## helper for simple / AHS: builds exactly one representation of D/Z,
  ## selected by linearD/linearZ (no more mixing multiple representations
  ## of the same condition-type within one call).
  make_linear_blocks <- function(cond_name) {
    if (!linearD && !linearZ) {
      ## ordinary margins: binarized Z and binarized D
      if (has_FE) {
        ## os_res is not computed under FE (see the one-sided-noncompliance
        ## block above) -- test every margin directly instead of filtering.
        tmp <- if (K > 2L && J > 2L) {
          data.table::CJ(zmargin = Zsup[-1L], dval = Dsup[-1L])
        } else if (K > 2L) {
          data.table::data.table(zmargin = Zsup[-1L])
        } else if (J > 2L) {
          data.table::data.table(dval = Dsup[-1L])
        } else {
          data.table::data.table()
        }

      } else {
        if (!exists("os_res")) {
          stop("Internal error: ordinary margin rows require `os_res`, but it was not computed.",
               call. = FALSE)
        }

        tmp <- os_res$threshold[one_sided == FALSE]

        keep <- intersect(c("zmargin", "dmargin"), names(tmp))
        tmp <- tmp[, ..keep]

        if ("dmargin" %in% names(tmp)) {
          data.table::setnames(tmp, "dmargin", "dval")
        }

        ## If D is binary, dval is not a meaningful margin for simple/AHS-none.
        ## Missing dval will later be treated as threshold 1 only inside Q construction.
        if (J == 2L && "dval" %in% names(tmp)) {
          tmp[, dval := NULL]
        }

        ## If Z is binary, zmargin is not a meaningful margin.
        if (K == 2L && "zmargin" %in% names(tmp)) {
          tmp[, zmargin := NULL]
        }
      }

      tmp[, condition := cond_name]
      tmp[, linear := "none"]

    } else if (linearD && !linearZ) {
      ## linear in D, margin in Z
      tmp <- if (K > 2L) {
        data.table::data.table(zmargin = Zsup[-1L])
      } else {
        data.table::data.table()
      }

      tmp[, condition := cond_name]
      tmp[, linear := "D"]

    } else if (!linearD && linearZ) {
      ## linear in Z, margin in D
      tmp <- if (J > 2L) {
        data.table::data.table(dval = Dsup[-1L])
      } else {
        data.table::data.table()
      }

      tmp[, condition := cond_name]
      tmp[, linear := "Z"]

    } else {
      ## linear in both Z and D
      tmp <- data.table::data.table()

      tmp[, condition := cond_name]
      tmp[, linear := "DZ"]
    }

    tmp
  }

  # simple
  if ("simple" %in% condition) {
    idx_blocks[["simple"]] <- make_linear_blocks("simple")
  }

  # AHS
  if ("AHS" %in% condition) {
    idx_blocks[["AHS"]] <- make_linear_blocks("AHS")
  }

  # MW
  if ("MW" %in% condition) {
    ## MW+FE is disallowed by the has_FE/"MW" %in% condition guard earlier
    ## in this function, so os_res is guaranteed to have been computed by
    ## the time this code runs -- check explicitly rather than relying on
    ## that guard staying in sync with this code across future edits (unlike
    ## simple/AHS/KR, MW has no has_FE fallback of its own, since it can
    ## never legitimately reach this point with has_FE = TRUE).
    if (has_FE || !exists("os_res")) {
      stop(
        "Internal error: MW margin rows require `os_res`, but it was not ",
        "computed (has_FE should have been disallowed for condition = \"MW\" ",
        "earlier).",
        call. = FALSE
      )
    }

    tmp <- os_res$threshold[one_sided == FALSE]

    keep <- intersect(c("zmargin", "dmargin", "direction"), names(tmp))
    tmp <- tmp[, ..keep]

    if ("dmargin" %in% names(tmp)) {
      data.table::setnames(tmp, "dmargin", "dval")
    }

    ## MW uses D.bin directly below. If D is binary, dval is not a meaningful
    ## MW margin and should not enter margin_index.
    if (J == 2L && "dval" %in% names(tmp)) {
      tmp[, dval := NULL]
    }

    ## If Z is binary, zmargin is not a meaningful MW margin.
    if (K == 2L && "zmargin" %in% names(tmp)) {
      tmp[, zmargin := NULL]
    }

    tmp[, condition := "MW"]

    eq_idx <- data.table::data.table(equation = 0:1, dummy__ = 1L)
    tmp[, dummy__ := 1L]

    tmp <- eq_idx[tmp, on = "dummy__", allow.cartesian = TRUE][
      , dummy__ := NULL
    ]

    if (nrow(os_res$threshold[one_sided == TRUE]) > 0L &&
        "direction" %in% names(tmp)) {
      tmp <- tmp[
        !(
          (direction == "Never-takers only"  & equation == 1L) |
            (direction == "Always-takers only" & equation == 0L)
        )
      ]
    }

    if ("direction" %in% names(tmp)) {
      tmp[, direction := NULL]
    }

    ## os_res$threshold[one_sided==FALSE] can be empty when every margin
    ## exhibits one-sided noncompliance, even though the top-level "all
    ## conditions trivially satisfied" stop above is skipped whenever KR is
    ## also requested (KR has its own, finer-grained escape hatch via
    ## os_res$exact that MW doesn't share). Catch that case explicitly here
    ## instead of silently producing a 0-row MW block that turns into
    ## NA-filled `condition` rows downstream via the margin_index join.
    if (nrow(tmp) == 0L) {
      stop(
        "One-sided noncompliance for all margins of Z and D - ",
        "condition = \"MW\" is trivially satisfied and there is nothing to test.",
        call. = FALSE
      )
    }

    idx_blocks[["MW"]] <- tmp
  }

  # KR
  if ("KR" %in% condition) {
    if (has_FE) {
      ## os_res is not computed under FE (see the one-sided-noncompliance
      ## block above) -- test every margin directly instead of filtering.
      tmp <- if (K > 2L) {
        data.table::CJ(zmargin = Zsup[-1L], dval = Dsup)
      } else {
        data.table::data.table(dval = Dsup)
      }

    } else {
      tmp <- os_res$exact[trivial_exact == FALSE]

      keep <- intersect(c("zmargin", "dval", "dmargin"), names(tmp))
      tmp <- tmp[, ..keep]

      if ("dmargin" %in% names(tmp)) {
        data.table::setnames(tmp, "dmargin", "dval")
      }

      ## If Z is binary, zmargin is not a meaningful KR margin.
      if (K == 2L && "zmargin" %in% names(tmp)) {
        tmp[, zmargin := NULL]
      }
    }

    tmp[, condition := "KR"]

    A_specs <- list()

    for (dv in Dsup) {
      if (joint == FALSE) {

        if (dv == min(Dsup)) {
          ## Match LATEtest Q0 singleton conditions:
          ## Q = -1(Y in {y}, D = dmin)
          A_list <- all_subsets(Ysup, min_size = 1L, max_size = 1L)

        } else if (dv == max(Dsup)) {
          ## Match LATEtest Q1 singleton conditions under current Q construction:
          ## Q = D - 1(Y in A, D = dmax)
          ##   = 1(D = dmax, Y notin A)
          ## so use complements of singletons.
          A_list <- lapply(Ysup, function(y) setdiff(Ysup, y))

        } else {
          ## For interior treatment values there is no binary-LATEtest analogue.
          ## If joint=FALSE is meant to be "singleton-style", use singleton sets.
          A_list <- all_subsets(Ysup, min_size = 1L, max_size = 1L)
        }

      } else if (dv == min(Dsup)) {

        A_list <- all_subsets(Ysup, min_size = 1L, max_size = 1L)

      } else if (dv == max(Dsup)) {

        A_list <- if (L >= 2L) {
          all_subsets(Ysup, min_size = 1L, max_size = L - 1L)
        } else {
          list()
        }

      } else {

        A_list <- all_subsets(Ysup, min_size = 1L, max_size = L)
      }

      if (length(A_list)) {
        A_specs[[as.character(dv)]] <- data.table::rbindlist(
          lapply(A_list, function(a) {
            lbl <- paste0("{", paste(a, collapse = ","), "}")
            data.table::data.table(
              dval = dv,
              yval = lbl,
              Avals = list(a)
            )
          }),
          use.names = TRUE,
          fill = TRUE
        )
      }
    }

    A_idx <- data.table::rbindlist(A_specs, use.names = TRUE, fill = TRUE)

    tmp <- A_idx[tmp, on = "dval", allow.cartesian = TRUE]

    idx_blocks[["KR"]] <- tmp
  }

  margin_index <- data.table::rbindlist(idx_blocks, use.names = TRUE, fill = TRUE)

  # Drop columns with no substantive variation in the margin index itself.
  # This prevents variables like linear = "none" only, or dval all NA, from
  # becoming margins later.
  drop_nonvarying_index_cols <- function(mi,
                                         always_keep = "condition",
                                         exclude = c("Avals")) {
    drop_cols <- setdiff(names(mi), c(always_keep, exclude))

    for (cc in drop_cols) {
      x <- mi[[cc]]
      n_nonmiss <- data.table::uniqueN(x, na.rm = TRUE)

      ## Drop if zero or one actual non-missing value.
      ## NA versus one real value is not treated as meaningful margin variation.
      if (n_nonmiss <= 1L) {
        mi[, (cc) := NULL]
      }
    }

    mi
  }

  margin_index <- drop_nonvarying_index_cols(margin_index)



  # --------------------------------------------------
  # Cross with zmargin-stacked data
  # --------------------------------------------------

  mi <- data.table::copy(margin_index)

  if ("zmargin" %in% names(data) && "zmargin" %in% names(mi)) {
    ## Make sure both join columns have the same type.
    data[, zmargin := as.numeric(zmargin)]
    mi[, zmargin := as.numeric(zmargin)]

    data <- mi[
      data,
      on = "zmargin",
      allow.cartesian = TRUE
    ]

  } else {
    data[, join_dummy__ := 1L]
    mi[, join_dummy__ := 1L]

    data <- mi[
      data,
      on = "join_dummy__",
      allow.cartesian = TRUE
    ]

    data[, join_dummy__ := NULL]
  }

  # --------------------------------------------------
  # Define margin variables
  # --------------------------------------------------

  nonredundant_margin_cols <- function(margin_index,
                                       data,
                                       exclude = c("Avals", "rowid", "join_dummy__")) {
    cand <- intersect(names(margin_index), names(data))
    cand <- setdiff(cand, exclude)

    cand[vapply(cand, function(cc) {
      x <- data[[cc]]

      ## Keep only if there are at least two actual non-missing values.
      ## NA versus one real value is not treated as meaningful variation.
      data.table::uniqueN(x, na.rm = TRUE) >= 2L
    }, logical(1))]
  }

  margins <- nonredundant_margin_cols(margin_index, data)

  ## ------------------------------------------------------------
  ## Drop margin cells with too few observations/clusters
  ## in either sample half
  ## ------------------------------------------------------------

  byvars <- c("sample", margins)

  if (!is.null(cluster)) {
    cell_size <- data[, .(
      n_obs = .N,
      n_clusters = data.table::uniqueN(.SD[[cluster]], na.rm = TRUE)
    ), by = byvars, .SDcols = cluster]
  } else {
    cell_size <- data[, .(
      n_obs = .N,
      n_clusters = .N
    ), by = byvars]
  }

  cell_size[, bad_size := n_obs < minsize | n_clusters < minsize]

  if (length(margins) == 0L) {
    if (cell_size[, any(bad_size)]) {
      warning(
        "Dropping all rows because at least one sample half has ",
        "fewer than minsize observations or clusters. ",
        "minsize = ", minsize,
        call. = FALSE
      )

      data <- data[0]

      if (exists("margin_index")) {
        margin_index <- margin_index[0]
      }
    }

  } else {
    bad_margins <- unique(cell_size[bad_size == TRUE, ..margins])

    if (nrow(bad_margins) > 0L) {
      bad_labels <- apply(bad_margins, 1L, function(r) {
        paste(paste(names(r), r, sep = "="), collapse = ", ")
      })

      warning(
        "Dropping ", nrow(bad_margins),
        " margin cell(s) because at least one sample half has ",
        "fewer than minsize observations or clusters. ",
        "minsize = ", minsize, ".\n",
        paste("  -", bad_labels, collapse = "\n"),
        call. = FALSE
      )

      bad_margins[, drop_bad_size__ := TRUE]

      data <- bad_margins[data, on = margins]
      data <- data[is.na(drop_bad_size__)][, drop_bad_size__ := NULL]

      if (exists("margin_index")) {
        margin_index <- bad_margins[margin_index, on = margins]
        margin_index <- margin_index[is.na(drop_bad_size__)][, drop_bad_size__ := NULL]
      }
    }
  }

  if (nrow(data) == 0L) {
    stop(
      "No remaining margin cells after minsize screening. ",
      "At least one sample half had fewer than minsize observations or clusters.",
      call. = FALSE
    )
  }


  # --------------------------------------------------
  # Create Q by condition
  # --------------------------------------------------

  Dcol <- Dbincol
  Ycol <- if (!is.null(Y)) paste0(Y, ".bin") else NULL
  Zcol <- Z

  data[, Q := NA_real_]

  linear_conditions <- c("simple", "AHS")

  is_linear_condition <- data[["condition"]] %in% linear_conditions
  is_MW <- data[["condition"]] == "MW"
  is_KR <- data[["condition"]] == "KR"

  if (!is.null(Dcol) && Dcol %in% names(data)) {
    Dvals__ <- sort(unique(stats::na.omit(data[[Dcol]])))
    D_is_binary__ <- length(Dvals__) <= 2L && all(Dvals__ %in% c(0, 1))
  } else {
    D_is_binary__ <- FALSE
  }


  # ---------------- simple / AHS: linear-D Q ----------------
  # Rows from a linear-eligible condition get raw D iff linearD = TRUE
  # (a single global choice -- no more per-row reconciliation needed).

  linear_D_rows <- is_linear_condition & linearD

  if (any(linear_D_rows)) {
    stopifnot(!is.null(Dname), Dname %in% names(data))

    data[
      linear_D_rows,
      Q := as.numeric(get(Dname))
    ]
  }


  # ---------------- simple / AHS: binarized-D Q ----------------
  # These are rows where Q should be 1(D >= dval).

  binarized_rows <- is_linear_condition & !linearD

  if (any(binarized_rows)) {
    stopifnot(!is.null(Dcol), Dcol %in% names(data))
    data[, dval_Q__ := NA_real_]

    if ("dval" %in% names(data)) {
      data[
        binarized_rows,
        dval_Q__ := as.numeric(dval)
      ]
    }

    if (D_is_binary__) {
      data[
        binarized_rows & is.na(dval_Q__),
        dval_Q__ := 1
      ]
    }

    if (any(is.na(data[binarized_rows, dval_Q__]))) {
      bad <- data[
        binarized_rows & is.na(dval_Q__),
        .N,
        by = intersect(
          c("condition", "linear", "dval", "zmargin", "yval", "equation"),
          names(data)
        )
      ]

      print(bad)

      stop(
        "dval is missing in rows that need binarized Q. ",
        "Because D is not binary, the margin index must provide dval for ",
        "ordinary simple/AHS rows with linearD = FALSE.",
        call. = FALSE
      )
    }

    data[
      binarized_rows,
      Q := as.numeric(get(Dcol) >= dval_Q__)
    ]

    data[, dval_Q__ := NULL]
  }

  # ---------------- MW ----------------

  if (any(is_MW)) {
    stopifnot(!is.null(Dcol), Dcol %in% names(data))
    stopifnot(!is.null(Zcol), Zcol %in% names(data))
    stopifnot(zhat %in% names(data))
    stopifnot("equation" %in% names(data))

    ## MW's Q statistic uses Z.hat directly as a propensity weight, with no
    ## other route through make_scores_vec()'s AIPW-branch clip/validation
    ## (simple/KR/AHS get that protection; MW previously had none at all).
    ## Apply the same hard-stop/clip-and-warn helper here.
    data[
      condition == "MW",
      (zhat) := validate_and_clip(
        get(zhat),
        hard_lower = 0, hard_upper = 1,
        floor = aipw.clip, ceiling = 1 - aipw.clip,
        label = "propensity e(X) for MW rows"
      )
    ]

    data[
      condition == "MW",
      Q := equation * (
        (1 - get(zhat)) * get(Dcol) * get(Zcol) -
          get(zhat) * get(Dcol) * (1 - get(Zcol))
      ) +
        (1 - equation) * (
          get(zhat) * (1 - get(Dcol)) * (1 - get(Zcol)) -
            (1 - get(zhat)) * (1 - get(Dcol)) * get(Zcol)
        )
    ]
  }

  # ---------------- KR ----------------

  if (any(is_KR)) {
    stopifnot("dval" %in% names(data))
    stopifnot("yval" %in% names(data))
    stopifnot("Avals" %in% names(data))
    stopifnot(!is.null(Dcol), Dcol %in% names(data))
    stopifnot(Ycol %in% names(data))

    dmin <- min(data[[Dcol]], na.rm = TRUE)

    data[
      condition == "KR" & dval == dmin,
      Q := -as.numeric(get(Ycol) %in% Avals[[1L]] & get(Dcol) == dval),
      by = .(dval, yval)
    ]

    data[
      condition == "KR" & dval > dmin,
      Q := as.numeric(get(Dcol) == dval) -
        as.numeric(get(Ycol) %in% Avals[[1L]] & get(Dcol) == dval),
      by = .(dval, yval)
    ]
  }

  # ---------------- sanity check ----------------

  if (any(is.na(data$Q))) {
    bad_Q <- unique(data[
      is.na(Q),
      intersect(
        c("condition", "linear", "equation", "zmargin", "dval", "yval"),
        names(data)
      ),
      with = FALSE
    ])

    print(bad_Q)

    stop("Some rows still have missing Q.", call. = FALSE)
  }

  ## ---------------- pre-check: drop cells with no variation in Q ----------------
  ## Checked here, before Q.hat is ever fit -- not just later, on the
  ## residual Q - Q.hat (see the bad_Q block further below). A cell with a
  ## completely constant Q carries no information for any downstream test,
  ## parametric or not, but under parametric = TRUE specifically feols()
  ## hard-errors trying to fit Q ~ X against a constant dependent variable,
  ## instead of the (still uninformative, but non-crashing) fit a forest
  ## would produce under parametric = FALSE. Screening it out up front makes
  ## both paths behave the same way instead of one of them crashing. Mirrors
  ## the later bad_Q block's own convention exactly (same "at least one
  ## sample part" grouping/wording, same drop-the-whole-margin-cell
  ## behavior), just applied to Q itself instead of its residual, before
  ## nuisance fitting instead of after.
  byvars_preQ <- c("sample", margins)
  data[, nQ_preQ := data.table::uniqueN(Q, na.rm = TRUE), by = byvars_preQ]

  if (length(margins) == 0L) {
    drop_all_preQ <- data[, any(nQ_preQ < 2L, na.rm = TRUE)]

    if (isTRUE(drop_all_preQ)) {
      message("Dropping all rows because at least one sample part has no usable variation in Q, before nuisance fitting.")
      data <- data[0]

      if (exists("margin_index")) {
        margin_index <- margin_index[0]
      }
    }

  } else {
    bad_preQ <- unique(data[nQ_preQ < 2L, ..margins])

    if (nrow(bad_preQ) > 0L) {
      if (length(margins) == 1L) {
        msg <- paste(bad_preQ[[margins]], collapse = ", ")
        message(
          "Dropping ", nrow(bad_preQ),
          " margin cell(s) because at least one sample part has no usable variation in Q, before nuisance fitting: ",
          msg
        )
      } else {
        bad_labels <- apply(bad_preQ, 1, function(r) {
          paste(paste(names(r), r, sep = "="), collapse = ", ")
        })
        message(
          "Dropping ", nrow(bad_preQ),
          " margin cell(s) because at least one sample part has no usable variation in Q, before nuisance fitting:\n",
          paste("  -", bad_labels, collapse = "\n")
        )
      }

      bad_preQ[, drop_preQ__ := TRUE]

      data <- bad_preQ[data, on = margins]
      data <- data[is.na(drop_preQ__)][, drop_preQ__ := NULL]

      if (exists("margin_index")) {
        margin_index <- bad_preQ[margin_index, on = margins]
        margin_index <- margin_index[is.na(drop_preQ__)][, drop_preQ__ := NULL]
      }
    }
  }

  data[, nQ_preQ := NULL]

  if (nrow(data) == 0L) {
    stop(
      "No remaining variation in Q for any margins. ",
      "Likely identification issue -- does X perfectly explain Z or D?",
      call. = FALSE
    )
  }

  time=rbind(time,"Stack data across margins"=proc.time())

  ## Estimate Q.hat after Q has been constructed
  data[, Q.hat := NA_real_]

  ## MW rows get no Q.hat fit at all -- MW's own score is assigned directly
  ## from Z.hat (see the MW Q block above) and never reads Q.hat, so fitting
  ## one would be pure wasted computation, and it would wrongly cost MW's
  ## test statistic degrees of freedom via x_rank_vars for a fit it never
  ## uses. Q.hat stays NA_real_ for MW rows (never consumed downstream:
  ## Q_needs_resid excludes MW from the sd_resQ check below, and MW's own
  ## scoring doesn't touch Q.hat).
  i_base <- which(!data$condition %in% c("MW", "AHS"))
  i_yres <- which(data$condition %in% c("AHS"))

  ## Rows with ordinary Q nuisance: Q ~ fml.Q | FE
  if (length(i_base)) {
    Q_hat_base <- estimate_conditional_mean(
      DT = data,
      y_name = "Q",
      x_expr = X_expr_Q,
      fe_expr = FE_expr,
      out_hat = "Q.hat",

      by = margins,
      sample_var = "sample",

      weight = weight,
      cluster = cluster,
      parametric = parametric,
      foldname = foldname,
      crossfit = crossfit,
      crossfit_label = "Q",

      forest_opts = Qparameters,
      fixest_opts = Qparameters,

      ## Reuse X_forest_info's FE-residualized columns when X_expr_Q is the
      ## same formula (and fixest_opts) as the main one -- see
      ## X_names_reuse_for_Q above. Otherwise let the wrapper build its own.
      x_names = X_names_reuse_for_Q,
      x_prefix = "__xq",

      keep_x = TRUE,
      return_residual = FALSE,
      partial_out_y_fe = TRUE,
      i = i_base
    )
  }

  if (length(i_yres)) {
    X_expr_Q_yres <- append_to_x_expr(X_expr_Q, y_name_rhs)

    ## If parametric = FALSE and you already have X_Q, you can pass
    ## c(X_Q, y_name_rhs). If not, let wrapper build from X_expr_Q_yres.
    x_names_yres <- NULL
    if (!isTRUE(parametric) && exists("X_Q", inherits = FALSE) && !is.null(X_Q)) {
      x_names_yres <- null_if_empty(c(X_Q, y_name_rhs))
    }

    Q_hat_yres <- estimate_conditional_mean(
      DT = data,
      y_name = "Q",
      x_expr = X_expr_Q_yres,
      fe_expr = FE_expr,
      out_hat = "Q.hat",

      by = margins,
      sample_var = "sample",

      weight = weight,
      cluster = cluster,
      parametric = parametric,
      foldname = foldname,
      crossfit = crossfit,
      crossfit_label = "Q",

      forest_opts = Qparameters,
      fixest_opts = Qparameters,

      x_names = x_names_yres,
      x_prefix = "__xq_yres",

      keep_x = TRUE,
      return_residual = FALSE,
      partial_out_y_fe = TRUE,
      i = i_yres
    )
  }

  ## Check for usable (residual) variation in Q
  byvars <- c("sample", margins)

  # common summaries
  data[, n_group := .N, by = byvars]
  data[, nQ      := data.table::uniqueN(Q, na.rm = TRUE), by = byvars]
  data[, sdQ     := stats::sd(Q, na.rm = TRUE), by = byvars]

  # family/type flags
  data[, Q_continuous :=
         condition == "MW" |
         (condition %in% c("simple", "AHS") & linearD)]
  data[, Q_needs_resid :=
         condition != "MW"]

  # residual variation for rows that need Q.hat
  data[, sd_resQ := NA_real_]
  if ("Q.hat" %in% names(data)) {
    data[Q_needs_resid == TRUE,
         sd_resQ := stats::sd(Q - `Q.hat`, na.rm = TRUE),
         by = byvars]
  }

  # discrete-support helper only where relevant
  data[, min_cell_Q := NA_integer_]
  data[Q_continuous == FALSE,
       min_cell_Q := {
         qtab <- table(Q, useNA = "no")
         if (length(qtab) == 0L) NA_integer_ else min(qtab)
       },
       by = byvars]

  # initialize
  data[, bad_Q := FALSE]

  # MW: continuous, no conditioning on Q.hat
  data[condition == "MW",
       bad_Q := (
         n_group < min_n |
           is.na(sdQ) | sdQ == 0
       ),
       by = byvars]

  # linear-D / linear-DZ: continuous and residualized
  data[Q_continuous == TRUE & Q_needs_resid == TRUE,
       bad_Q := (
         n_group < min_n |
           is.na(sd_resQ) | sd_resQ == 0
       ),
       by = byvars]

  # KR / simple / AHS / simple_linearZ: discrete and residualized
  data[Q_continuous == FALSE,
       bad_Q := (
         n_group < min_n |
           nQ < 2 |
           is.na(min_cell_Q) | min_cell_Q < min_n |
           is.na(sd_resQ) | sd_resQ == 0
       ),
       by = byvars]

  ## ------------------------------------------------------------
  ## Drop entire margin cells if any sample part / condition cell fails
  ## ------------------------------------------------------------

  if (length(margins) == 0L) {
    drop_all <- data[, any(bad_Q, na.rm = TRUE)]

    if (isTRUE(drop_all)) {
      message("Dropping all rows because at least one sample part has no usable variation in Q.")
      data <- data[0]

      if (exists("margin_index")) {
        margin_index <- margin_index[0]
      }
    }

  } else {
    bad_margins <- unique(
      data[bad_Q == TRUE, ..margins]
    )

    if (nrow(bad_margins) > 0L) {
      if (length(margins) == 1L) {
        msg <- paste(bad_margins[[margins]], collapse = ", ")
        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part has no usable variation in Q: ",
          msg
        )
      } else {
        bad_labels <- apply(bad_margins, 1, function(r) {
          paste(paste(names(r), r, sep = "="), collapse = ", ")
        })
        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part has no usable variation in Q:\n",
          paste("  -", bad_labels, collapse = "\n")
        )
      }

      bad_margins[, drop_bad_Q__ := TRUE]

      data <- bad_margins[data, on = margins]
      data <- data[is.na(drop_bad_Q__)][, drop_bad_Q__ := NULL]

      if (exists("margin_index")) {
        margin_index <- bad_margins[margin_index, on = margins]
        margin_index <- margin_index[is.na(drop_bad_Q__)][, drop_bad_Q__ := NULL]
      }
    }
  }

  if (nrow(data) == 0L) {
    stop(
      "No remaining residual variation for any margins. ",
      "Likely identification issue - does X perfectly explain Z or D?",
      call. = FALSE
    )
  }
  ## ------------------------------------------------------------
  ## Cleanup helper columns
  ## ------------------------------------------------------------

  helper_cols_Q <- intersect(
    c("n_group", "nQ", "sdQ", "Q_continuous", "Q_needs_resid",
      "sd_resQ", "min_cell_Q", "bad_Q"),
    names(data)
  )

  if (length(helper_cols_Q) > 0L) {
    data[, (helper_cols_Q) := NULL]
  }

  time=rbind(time,"Estimate nuisance for outcomes Q"=proc.time())

  ########## ESTIMATE ALL CAUSAL/REGRESSION/IV FORESTS AND  predict in/out of sample ##########
  if (!"C" %in% crossfit) foldname=NULL #Do not crossfit causal forest, just the nuisances - use OOB for forest.


  ### SIMPLE, KR ###

  i_simple_kr <- which(data$condition %in% c("simple", "KR"))

  if (length(i_simple_kr)) {
    fit_models(
      DT = data,
      i = i_simple_kr,

      forest_type = "causal",
      y_name = "Q",
      x_names = X_forest,

      margins = margins,

      w_name = Z,

      folds = foldname,
      weight_name = weight,
      cluster_name = cluster,

      forest_opts = Cparameters,
      aipw.clip = aipw.clip,
      shrink = (shrink > 0),
      verbose = FALSE,

      doubly.robust = doubly.robust,
      z_linear_score_name = "z_use_linear_score"
    )
  }

  ##AHS

  i_ahs <- which(data$condition == "AHS")

  if (length(i_ahs)) {
    fit_models(
      DT = data,
      i = i_ahs,

      forest_type = "causal",
      y_name = "Q",
      x_names = null_if_empty(c(X_forest, y_name_rhs)),

      margins = margins,

      w_name = Z,

      folds = foldname,
      weight_name = weight,
      cluster_name = cluster,

      forest_opts = Cparameters,
      aipw.clip = aipw.clip,
      shrink = (shrink > 0),
      verbose = FALSE,

      doubly.robust = doubly.robust,
      z_linear_score_name = "z_use_linear_score"
    )
  }

  ## ------------------------------------------------------------
  ## MW: no target; conditional mean E[Q | X]
  ## ------------------------------------------------------------
  if (any(data$condition == "MW")) {
    i_mw <- which(data$condition == "MW")
    data[i_mw, scores := Q]

    if (testtype == "forest") {
      fit_models(
        DT = data,
        i = i_mw,
        forest_type = "regression",
        y_name = "Q",
        x_names = null_if_empty(c(X_forest, y_name_rhs)),
        folds = foldname,
        margins = margins,
        weight_name = weight,
        cluster_name = cluster,
        forest_opts = Cparameters,
        shrink = (shrink > 0)
      )
    }
  }

  ## ------------------------------------------------------------
  ## Target-weighting: `target` no longer changes the score formula (see
  ## make_scores_vec()) -- it only changes the weight fed into the final
  ## forest_test()/CART_test() moment. target == "overlap" reweights by the
  ## per-row conditional variance `scores_v` (attached by fit_models() for
  ## simple/KR/AHS); target == "all" leaves the weight untouched. MW is
  ## deliberately excluded -- it was never touched by `target` even before
  ## this refactor (its own plug-in moment doesn't go through
  ## make_scores_vec() at all).
  ## ------------------------------------------------------------

  data[, w_eff := if (is.null(weight)) 1 else get(weight)]

  if (identical(target, "overlap") && "scores_v" %in% names(data)) {
    data[
      condition %in% c("simple", "KR", "AHS"),
      w_eff := w_eff * scores_v
    ]
  }

  ## ------------------------------------------------------------
  ## Test-side centering: for need_ols_v rows only (doubly.robust=FALSE &
  ## target=="overlap"), the through-origin score's ratio structure only
  ## reproduces classical FWL/OLS at the level `Z.hat` was normalized at --
  ## an adaptively-found leaf nested inside that level isn't guaranteed to
  ## have zero mean(Z-Z.hat), so forest_test()/CART_test()'s final held-out
  ## test uses crv1_mean(..., center=TRUE) instead, re-fitting a with-
  ## intercept coefficient fresh on each candidate leaf's own rows. See
  ## crv1_mean()'s `center` argument and montest.R's `need_ols_v`. Applies
  ## only to the final test-side moment, never to the search/screening
  ## phase; and only to conditions that actually route through
  ## make_scores_vec() (excludes MW, matching the w_eff gate above), and
  ## only when those conditions' raw residual columns were populated.
  ## ------------------------------------------------------------

  resid_cols_exist <- all(c(".resid_treat", ".resid_outcome") %in% names(data))

  data[, use_centered_test := FALSE]
  if (need_ols_v && resid_cols_exist) {
    data[
      condition %in% c("simple", "KR", "AHS"),
      use_centered_test := TRUE
    ]
  }

  ## Only pass the centering column names through to forest_test()/
  ## CART_test() when they actually exist (e.g. not for an MW-only run,
  ## which never populates `.resid_treat`/`.resid_outcome`) -- otherwise
  ## `use_centered_test` is all-FALSE and these arguments would be dead
  ## weight, but forest_test_core()/CART_test() validate any non-NULL
  ## `resid_treat`/`resid_outcome` column name actually exists in `data`.
  center_arg <- if (resid_cols_exist) "use_centered_test" else NULL
  resid_treat_arg <- if (resid_cols_exist) ".resid_treat" else NULL
  resid_outcome_arg <- if (resid_cols_exist) ".resid_outcome" else NULL

  ## `stabilize.scores`'s narrower test-side fix (see its own docs): reuses
  ## the same `.resid_treat`/`.resid_outcome` columns above, plus `scores_v`
  ## (`v`, always populated alongside them by the same make_scores_vec()
  ## call) and, for CART_test() only, `pred` as the tau source (forest_test()
  ## already defaults its own `pred` argument to "pred", so no separate
  ## argument is needed there). `recenter_binary_arg` picks which of the two
  ## corrections gets redone locally, matching whichever one is actually in
  ## use for these rows -- reuses `z_use_lin_any` (see `need_v_hat` above),
  ## since the continuous/binary representation choice is a single
  ## montest()-call-level constant, never a per-row one (see
  ## `stabilize.scores`'s docs).
  v_arg <- if (resid_cols_exist) "scores_v" else NULL
  tau_arg <- if (resid_cols_exist) "pred" else NULL
  recenter_propensity_arg <- isTRUE(stabilize.scores)
  recenter_binary_arg <- !z_use_lin_any

  ###EMPIRICAL BAYES SHRINKAGE IF SHRINK>0 #######

  if (shrink>0&testtype=="forest") {
    shrink_te_crossfit(
    data        = data,
    pred        = "pred",
    pred_var     = "pred_var",
    pred_out    = "pred_o",
    pred_out_var = "pred_o_var",
    margins     = margins,
    sample      = "sample",
    gamma       = shrink
  )
  }


  time=rbind(time,"Estimate causal forests"=proc.time())

  ######################################## FIND OPTIMAL SUBSET TO TEST AND TEST IN OPPOSITE SAMPLE #####################
  poolmargins=pool[pool %in% c(margins,"sample")]
  selectmargins=select[select %in% c(margins,"sample")]

  if ("forest" == testtype) res=forest_test(data,cluster=cluster,weight="w_eff",minsize=minsize,x_names=X_forest,pool=poolmargins,select=selectmargins,gridpoints=gridpoints,margins=margins,screen=screen,alpha=alpha,fe_expr=FE_expr,fe_rank_adj=fe_rank_adj,fe_rank_conservative = fe_rank_conservative,x_rank_vars=x_rank_vars,center=center_arg,resid_treat=resid_treat_arg,resid_outcome=resid_outcome_arg,sample_weight=weight,recenter_propensity=recenter_propensity_arg,recenter_binary=recenter_binary_arg,v=v_arg)
  if ("CART" == testtype) res=CART_test(data, x_names=X_forest,margins=margins,weight="w_eff",cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,screen=screen,cluster=cluster,select=selectmargins,rpart_options=Rparameters,fe_expr=FE_expr,fe_rank_adj=fe_rank_adj,x_rank_vars=x_rank_vars,center=center_arg,resid_treat=resid_treat_arg,resid_outcome=resid_outcome_arg,sample_weight=weight,recenter_propensity=recenter_propensity_arg,recenter_binary=recenter_binary_arg,tau=tau_arg,v=v_arg)

  ## Xmeans/Xmeans_all/XSD (when present) are keyed on the internal
  ## `__xf_raw_*`/`__xf_res_*` forest-feature columns from
  ## make_X_residualized_from_FE() -- FE-residualized whenever `fml` has an
  ## FE part -- rather than the original covariate names. Rename them back
  ## to the clean (model-matrix) covariate names for output, so downstream
  ## consumers like montestplot() show the real names. Guard against a
  ## clean name colliding with the pool/select vocabulary or with another
  ## covariate's cleaned name, either of which would make the rename
  ## ambiguous to name-based lookups; such columns keep their internal name.
  if (!is.null(X_forest_info) && !is.null(X_forest_info$x_names)) {
    reserved <- c("zmargin", "dval", "yval", "condition", "equation", "sample")
    old_nm <- X_forest_info$x_names
    new_nm <- X_forest_info$clean_names
    bad <- new_nm %in% reserved | duplicated(new_nm) | duplicated(new_nm, fromLast = TRUE)
    if (any(bad)) {
      warning(
        "Could not restore original covariate name(s) for: ",
        paste(unique(new_nm[bad]), collapse = ", "),
        " (collides with a pool/select dimension name or another covariate's ",
        "cleaned name); keeping the internal column name in Xmeans/Xmeans_all/XSD instead.",
        call. = FALSE
      )
      old_nm <- old_nm[!bad]
      new_nm <- new_nm[!bad]
    }
    for (xtbl in c("Xmeans", "Xmeans_all", "XSD")) {
      if (!is.null(res[[xtbl]])) {
        keep <- old_nm %in% names(res[[xtbl]])
        if (any(keep)) data.table::setnames(res[[xtbl]], old_nm[keep], new_nm[keep])
      }
    }
  }

  time=rbind(time,"Find promising subset and test"=proc.time())


  ################ 7: Multiple hypothesis testing and output #####################
    if (nrow(res$results[train==FALSE])==1) {
      res$minwhere=NA
      res$minp=rep(res$results[train==FALSE,p.raw],6)
      names(res$minp)=paste0("p.",c("raw","holm","hochberg","BH","BY","CCT"))
    } else {
      for (m in c("holm","hochberg","BH","BY")) {
        res$results[train==FALSE,paste0("p.",m):=p.adjust(replace(p.raw, is.na(p.raw), 1),method=m)]
      }
      byv=c("sample",margins)[!c("sample",margins) %in% pool]
      res$minwhere=res$results[train == FALSE & is.finite(p.raw)][which.min(p.raw), ..byv]
      res$minp=apply(res$results[train==FALSE,c("p.raw","p.holm","p.hochberg","p.BH","p.BY")],2,min)
      res$minp=c(res$minp,p.CCT=cct_pvalue(replace(res$results[train==FALSE,p.raw],is.na(res$results[train==FALSE,p.raw]),1)))
    }

    res$global[,p.raw:=pnorm(t)]
    if (nrow(res$global)>1) {
      for (m in c("holm","hochberg","BH","BY")) {
        res$global[,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }
    }


  time=rbind(time,"Correct for multiple hypothesis testing"=proc.time())
  time = time[-1, , drop = FALSE] - time[-nrow(time), , drop = FALSE]
  time=time[,1:3]
  time=rbind(time,"Total"=colSums(time))

  ## `options`: the *resolved* value of every montest() argument (after all
  ## match.arg()/default-resolution logic above has run), keyed by pulling
  ## each formal's current binding out of this call frame -- so it stays in
  ## sync automatically as arguments are added/changed, without hand-
  ## maintaining a list. `data` is excluded: by this point it is the fully
  ## expanded/mutated internal working table, not the caller's original
  ## data, and would needlessly bloat the returned object.
  options <- mget(setdiff(names(formals(montest)), "data"), envir = environment())

  out <- c(res, list(
    call = mc,
    options = options,
    time = time,
    obs = obs,
    dropped = c(fe_singletons = n_dropped_singletons, novar_Z = n_dropped_novar_Z),
    margins = margin_index[, setdiff(names(margin_index), "Avals"), with = FALSE]
  ))
  return(out)
}
