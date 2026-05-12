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
#' are expanded across margins and discretized into bins before estimation. See Details.
#'
#' @param data A \code{data.frame} or \code{data.table} containing the analysis sample.
#'   Observations with missing values in any variables used by the call are dropped.
#' @param fml A fixest-style formula on the form Y~X|fe|D~Z that specified the IV call to be tested. Importantly, the instrument Z should be coded so that higher values of Z imply weakly higher values of D for everyone - positive monotonicity.
#' The formula may be one-sided and omit Y if testing only the simple first stage condition. Note that the exact functional form does not matter in the default case when \code{parametric=FALSE} because the command uses semiparametric methods.
#' @param parametric A boolean indicating whether nuisances should be estimated using the parametric functional form specified or (the default) using semiparametric methods. In the latter case,
#' all fixed effects are one-hot encoded, while the functional form in the main part of the formula is ignored and determined by the corresponding regression forests for the nuisance parameters.
#' @param one.hot boolean indicating whether fixed effects should be one-hot encoded or added as a continuous variable. Default is TRUE.
#' @param condition Character vector selecting which tests to run. Allowed values are any combination of
#'   \code{"simple"}, \code{"KR"} (Kwan-Roth conditions), \code{"MW"} (Mourifie and Wan conditions), \code{"AHS"} (Andresen-Huber-Sloczynski), or \code{"all"}.
#'   If \code{Y} is omitted, only \code{"simple"} is allowed.
#'  @param target a scalar equal to either "all" or  "overlap" to determine whether the function should target the average treatment effect (all) in the testing subsample, or the overlap-adjusted / weighted ATE ("overlap"). For linear Z versions, overlap targets the residualized slope, while all targets the average partial effect.
#' @param inner.folds Optional integer giving the number of within-sample folds used for
#'   cross-fitting nuisance functions and, optionally, forest predictions. Set to
#'   \code{NULL} to disable the inner split. Defaults to NULL - nuisances and predictions from causal forests are fit out-of-bag. See option crossfit, which decides which parts this applies to.
#' @param crossfit Character vector of what parts of the procedure to cross-fit. Accepts "Z","Q","Y","C". If e.g. "Z" appears in crossfit, nuissances for Z are cross fit, either across outer sample part (if inner.folds==NULL), or within outer sample part across inner folds. If "Z" does not appear, OOB predictions are used. "C" is for the causal forest fit.
#' @param normalize.Z Logical, default TRUE; if \code{TRUE}, estimated instrument propensity scores are
#'   normalized after estimation.
#' @param aipw.clip Positive scalar in \code{(0,1)}, default 1e-3, used to trim estimated propensity
#'   scores when augmented inverse-probability weighted scores are constructed and when normalizing propensity scores.
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
#' @param gridpoints Optional integer controlling the number of candidate cutoffs searched
#'   by the forest-based test. If \code{NULL}, all eligible cutoffs are considered.
#' @param min_n Integer minimum number of treated and untreated instrument observations
#'   required within each sample half and margin cell considered for residual variation.
#' @param pool Character vector controlling which dimensions are pooled when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dval"},
#'   \code{"yval"}, \code{"condition"}, \code{"equation"},
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. No margin can appear in both pool and select, but "sample" can, implying adaptive pooling. Relevant margins that appear in neither are all tested, and tests are corrected for multiple hypothesis testing.
#' @param select Character vector controlling which dimensions are selected over when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dval"},
#'   \code{"yval"}, \code{"condition"}, \code{"equation"},
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. No margin can appear in both pool and select, but "sample" can, implying adaptive pooling. Relevant margins that appear in neither are all tested, and tests are corrected for multiple testing.
#' @param screen Screening rule for deciding what determines a "promising" leaf or cell to carry forward to testing. May be "minimum","negative","nonpositive","stepdown","fg_relevant","none". Defaults to stepdown, described below.
#' @param cp,maxrankcp,alpha,prune Tuning parameters for the CART-based search
#'   routine. See Details.
#' @param Zparameters,Yparameters,Qparameters,Cparameters,Rparameters Named lists of
#'   additional arguments passed to the underlying estimation routines for different
#'   nuisance or target models. See regression_forest, causal_forest, feols and rpart for for details.
#'   @param joint specifies that all Kwan-Roth conditions should be included in the test, not only those for which the subset A contains only one outcome value. Defaults to TRUE.
#'
#' @details
#' The rough steps of the montest algorithm is as follows
#'
#' \enumerate{
#'   \item Data is restricted to complete cases of the variables used, and the sample is split into two halves (optionally respecting clusters).
#'   \item Continuous or multivalued instruments, treatments and outcomes discretized into bins.
#'   \item The data is stacked across instrument margins, treatment margins, outcomes,
#'   equations, and test conditions. The outcome variable Q is defined, depending on condition and margins.
#'   \item Nuisance functions for as \code{Z.hat}, and the outcome \code{Q.hat} are estimated and outcomes are residualized
#'   using regression forests within margins.
#'   \item Separate causal forests of the outcome \code{Q}, on
#'   the instrument \code{Z} using features \code{X} (and optionally \code{Y} for MW and AHS conditions),
#'   treatment effects are predicted in and out of sample and scores constructed
#'   \item Each sample part (optionally within margins, depending on the options in \code{pool}) is sorted
#'   according to treatment effects, and the mean of scores is estimated numerically for all possible cutoffs
#'   in predicted treatment effects. Select the cutoff with the smallest t-statistic on the mean of scores
#'   Alternatively, subset selection can be done using a CART algorithm.
#'   \item Depending on the choices in \code{select}, promising subsets (using selection rules in \code{screen}) are evaluated  in the opposite sample half. If performing multiple tests,
#'   depending on the option \code{pool}, p-values are adjusted for multiple testing.
#' }
#'
#' \code{montest} supports weights and clustering, and allow for multivalued treatments and instruments
#' by binarizing instruments and treatments into quantile or equisized bins. The command also tests for
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
#'   \item \code{condition = "MW} tests the sharp conditions from Mourifie and Wan (2017), which tests monotonicity
#'   and exclusion conditional on instrument validity. This is only allowed for a binary treatment. There are a total of 2K such conditions.
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

montest=function(fml,data,condition=NULL,inner.folds=NULL,crossfit=NULL,
                 normalize.Z=FALSE,aipw.clip=1e-3,weight=NULL,cluster=NULL,seed=10101,minsize=50L,
                 gridtypeY="equidistant",gridtypeD="equisized",gridtypeZ="equisized",stratify=TRUE,joint=TRUE,
                 Ysubsets = 4L, Dsubsets = 4L,Zsubsets=4L,Y.res=TRUE,testtype="forest",
                 gridpoints=NULL,min_n=1L,pool=NULL,select=NULL,shrink=0,linear="none",target="all",one.hot=TRUE,
                 cp=0,maxrankcp=10L,Rparameters=list(),alpha=0.05,prune=TRUE,screen="stepdown",parametric=FALSE,
                 Zparameters=list(),Yparameters=list(),Qparameters=list(),Dparameters=list(),Cparameters=list()
){

  time=rbind(start=proc.time())
  if (!is.null(seed)) set.seed(seed)


  ################### 1 CHECK INPUT #####################
  data <- data.table::as.data.table(data.table::copy(data))
  target=match.arg(target,c("all","overlap"))
  ##Check formula and validate
  v=validate_iv(fml,data)
  X=v$X
  Y=v$Y
  D=v$D
  Z=v$Z
  z_col=Z
  d_col=D
  y_col=Y
  FE=v$FE
  FE_expr=v$FE_expr
  fml_rf=v$reduced_form_rhs

  ##One-hot encode FE or add to X
  if (!is.null(FE_expr)) {
    if (one.hot) {
      fe_out <- fixest_fe_to_mm(
        DT = data,
        fe_expr = FE_expr,
        prefix = "fe",
        drop_first = TRUE
      )

      X <- null_if_empty(unique(c(X, fe_out$new_cols)))

    } else {
      X <- null_if_empty(unique(c(X, all.vars(FE_expr))))
    }
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

  linear=match.arg(linear,c("none","Z","D","DZ","all"),several.ok=TRUE)
  if ("all" %in% linear) {
    linear <- c("none", "Z", "D", "DZ")
  }
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
  else if (sum(pool=="all")==1) pool=c("zmargin","dval","yval","condition","equation","sample","linear")
  else if (sum(pool=="none")==1) pool=c()
  else if (is.null(pool)==FALSE) pool <- match.arg(
    pool,
    c("zmargin", "dval", "yval", "condition", "equation", "sample", "linear"),
    several.ok = TRUE
  )
  else pool=c("zmargin","dval","yval","sample","linear")

  if ((sum(select=="none")==1)&(sum(select=="all")==1)) stop("Do not specify both none and all in select().")
  else if (sum(select=="all")==1) select=c("zmargin","dval","yval","condition","equation","sample","linear")
  else if (sum(select=="none")==1) select=c()
  else if (is.null(select)==FALSE) select=match.arg(select,c("zmargin","dval","yval","condition","equation","sample","linear"),several.ok=TRUE)
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

  ##Downgrade linear options if D or Z is binary:
  ## Number of support points in original D and Z
  Zname=Z
  Dname=D

  J <- length(unique(data[,get(d_col)]))
  K <- length(unique(data[,get(z_col)]))

  normalize_linear <- function(linear, J, K) {
    linear <- unique(linear)

    out <- character()

    for (ll in linear) {
      if (ll == "none") {
        out <- c(out, "none")

      } else if (ll == "Z") {
        if (K > 2L) {
          out <- c(out, "Z")
        } else {
          out <- c(out, "none")
        }

      } else if (ll == "D") {
        if (J > 2L) {
          out <- c(out, "D")
        } else {
          out <- c(out, "none")
        }

      } else if (ll == "DZ") {
        z_lin <- K > 2L
        d_lin <- J > 2L

        if (z_lin && d_lin) {
          out <- c(out, "DZ")
        } else if (z_lin && !d_lin) {
          out <- c(out, "Z")
        } else if (!z_lin && d_lin) {
          out <- c(out, "D")
        } else {
          out <- c(out, "none")
        }

      } else {
        stop("Unknown linear option: ", ll)
      }
    }

    unique(out)
  }

  linear <- normalize_linear(linear, J = J, K = K)

  ## Conditions for which linear has an effect
  linear_conditions <- c("simple", "AHS")

  has_linear_conditions <- any(condition %in% linear_conditions)
  has_other_conditions  <- any(!condition %in% linear_conditions)

  ## If any non-linear-aware condition is requested, ordinary margins must also be present.
  ## Otherwise KR/MW/etc. would have no valid non-linear block to use.
  if (has_other_conditions && !"none" %in% linear) {
    stop(
      "When condition includes conditions other than ",
      paste(linear_conditions, collapse = " or "),
      ", linear must include 'none'. ",
      "The non-linear-aware conditions ignore linearized variants and require ordinary margin rows. ",
      "Use e.g. linear = c('none', ",
      paste(sprintf("'%s'", setdiff(linear, "none")), collapse = ", "),
      ")."
    )
  }

  ## Inform the user that linear variants only apply to simple/AHS.
  if (has_other_conditions && any(linear != "none")) {
    msg_other <- paste(setdiff(condition, linear_conditions), collapse = ", ")
    msg_linear <- paste(setdiff(linear, "none"), collapse = ", ")

    message(
      "Note: linear option(s) ",
      msg_linear,
      " only apply to condition(s) ",
      paste(linear_conditions, collapse = ", "),
      ". For the other requested condition(s) ",
      msg_other,
      ", only linear = 'none' will be used."
    )
  }

  ## Conditions that always need margin/binary versions
  nonlinear_conditions <- setdiff(condition, linear_conditions)

  ## Need binarized Z if:
  ## - any non-linear-aware condition is requested, or
  ## - simple/AHS are requested with linear = none or D
  need_binarized_Z <-
    length(nonlinear_conditions) > 0L ||
    any(condition %in% linear_conditions) && any(linear %in% c("none", "D"))

  ## Need original/linear Z if:
  ## - simple/AHS are requested with linear = Z or DZ
  need_linear_Z <-
    any(condition %in% linear_conditions) && any(linear %in% c("Z", "DZ"))

  ## Need margin conditions if:
  ## - ordinary simple/AHS rows are requested, linear = "none"
  ## - linear-D simple/AHS rows are requested, linear = "D"
  ## - KR/MW are requested
  ##
  ## These rows require the one-sided monotonicity pre-check, which itself
  ## needs the binarized treatment margin variable D.bin.
  has_margin_conditions <-
    any(c("none", "D") %in% linear) ||
    any(c("KR", "MW") %in% condition)

  ## Need binarized D for Q construction if:
  ## - any non-linear-aware condition is requested, or
  ## - simple/AHS are requested with linear = none or Z
  need_binarized_D_for_Q <-
    length(nonlinear_conditions) > 0L ||
    (any(condition %in% linear_conditions) && any(linear %in% c("none", "Z")))

  ## Need binarized D for the one-sided monotonicity screen.
  need_binarized_D_for_os <- has_margin_conditions

  ## Actually create D.bin if either Q construction or the one-sided screen needs it.
  need_binarized_D <- need_binarized_D_for_Q || need_binarized_D_for_os

  ## Need original/linear D if:
  ## - simple/AHS are requested with linear = D or DZ
  need_linear_D <- any(condition %in% linear_conditions) && any(linear %in% c("D", "DZ"))

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

    allowed_select <- c("zmargin","dval","yval","condition","equation","sample","linear")
    bad_select <- setdiff(select, allowed_select)
    if (length(bad_select) > 0L) {
      stop(
        "Invalid `select` for testtype = \"CART\". ",
        "For CART, `select` may only contain \"zmargin\",\"dval\",\"yval\",\"condition\",\"equation\",\"sample\",\"linear\". ",
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


  allvars <- unique(c(X, Y, D, Z, weight, cluster,FE))
  allvars <- allvars[!is.na(allvars)]
  data <- data[, ..allvars]
  dropped=sum(!complete.cases(data))
  if (dropped>0) message(paste("Note: dropped ",dropped," observations with missing data on one or more input variables."))
  data=data[complete.cases(data)]
  n=nrow(data)
  if (is.null(cluster)==FALSE) {
    G <- uniqueN(data[[cluster]])
    if (G<=2*minsize) stop("Number of clusters is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
  } else {
    if (n<=2*minsize) stop("Number of observations is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
    G=NA
    }
  obs=c(N=n,G=G)
  data[, id_ := .I]


  ##OUTER SPLIT
  stopifnot(is.logical(stratify), length(stratify) == 1L, !is.na(stratify))

  strat <- NULL

  if (isTRUE(stratify) && is.null(cluster)) {
    strat <- Z
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

  if (J>2&("MW" %in% condition)) stop("Multivalued treatment not supported with condition MW.")
  if (J==2&K==2&is.null(X)==TRUE&!any(condition %in% c("AHS","MW","KR"))) {
    stop("Nothing to test with a binary treatment, a binary instrument, the simple first stage condition and no variables in X.")
  }

  time=rbind(time,"Binarize treatment, instrument and outcome"=proc.time())

  ######################## 6a STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
  margins=c()

  if (K > 2L && (need_binarized_Z || need_linear_Z)) {

    zmap_list <- list()

    ## Margin-specific expansion.
    ##
    ## - lowest Z gets second-lowest margin, so Zbin >= zmargin gives 0
    ## - highest Z gets highest margin, so Zbin >= zmargin gives 1
    ## - intermediate Z gets current and next-higher margins, giving one 1 and one 0
    if (need_binarized_Z) {
      stopifnot(!is.null(Zbincol), Zbincol %in% names(data))

      zvals <- sort(unique(data[[Zbincol]]))
      zvals <- zvals[!is.na(zvals)]

      if (length(zvals) < 2L) {
        stop("Need at least two non-missing Z values to construct zmargin.")
      }

      minZ <- zvals[1L]
      maxZ <- zvals[length(zvals)]

      zmap_list[["margin"]] <- data[, {
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

        .(
          rowid = .I,
          zmargin = zm,
          z_is_linear = FALSE
        )
      }, by = id_]
    }

    ## Linear-Z copy.
    ##
    ## One copy per original row, with zmargin = NA.
    ## Later, these rows should only match linear %in% c("Z", "DZ")
    ## in margin_index.
    if (need_linear_Z) {
      zmap_list[["linear"]] <- data[
        ,
        .(
          rowid = .I,
          zmargin = NA_real_,
          z_is_linear = TRUE
        ),
        by = id_
      ]
    }

    zmap <- data.table::rbindlist(zmap_list, use.names = TRUE)

    data <- data[zmap$rowid]
    data[, zmargin := zmap$zmargin]
    data[, z_is_linear := zmap$z_is_linear]

    margins <- unique(c(margins, "zmargin"))

    ## Assign Z on the expanded rows.
    ##
    ## Margin rows use binarized Z at the local zmargin.
    ## Linear rows keep the original Z.
    data[z_is_linear == FALSE, (Z) := as.integer(get(Zbincol) >= zmargin)]
    data[z_is_linear == TRUE,  (Z) := get(Zname)]

    if (!is.null(Zbincol) && Zbincol %in% names(data)) {
      data[, (Zbincol) := NULL]
    }

  } else if (need_linear_Z) {

    ## No Z-margin stacking needed, but linear-Z designs need original Z.
    data[, (Z) := get(Zname)]
    data[, zmargin := NA_real_]
    data[, z_is_linear := TRUE]

    if (!is.null(Zbincol) && Zbincol %in% names(data)) {
      data[, (Zbincol) := NULL]
    }

  } else if (need_binarized_Z) {

    ## Binary-Z case or otherwise no stacking needed.
    ## If Zbincol exists, use it as the working binary Z.
    if (!is.null(Zbincol) && Zbincol %in% names(data)) {
      data[, (Z) := get(Zbincol)]
      data[, (Zbincol) := NULL]
    }

    data[, zmargin := NA_real_]
    data[, z_is_linear := FALSE]
  }

  ## Safety fallback if none of the branches above ran.
  if (!("zmargin" %in% names(data))) {
    data[, zmargin := NA_real_]
  }

  if (!("z_is_linear" %in% names(data))) {
    data[, z_is_linear := FALSE]
  }

  ## Conditional variance of Z is only needed for continuous/linear-Z rows
  ## when the global target is "all".
  need_z_var_global <- identical(target, "all") && any(data$z_is_linear == TRUE, na.rm = TRUE)
  need_z_var_rows <- if (need_z_var_global) which(data$z_is_linear == TRUE) else integer()

  ##estimate Z.hat for each margin in stacked data
  zhat <- paste0(Z, ".hat")
  if (!is.null(X)) {
    if (parametric==FALSE) {
      crossfit_hat(
      data,
      y_name = Z,
      x_names = X,
      margins = margins,
      weight_name = weight,
      forest_opts = Zparameters,
      folds=foldname,
      mode=ifelse(("Z" %in% crossfit & is.null(foldname)==TRUE),"across","within")
    )
    } else {
      parametric_hat(
              fml_rf,
              data=data,
             outcome_name=Z,
             margins = margins,
             weight_name = weight,
             fixest_opts = Zparameters)
    }
  } else { ##leave-cluster out mean
    leave_cluster_out_mean(
      DT = data,
      y_name = Z,
      cluster_name = cluster,
      by = c("sample", margins),
      out_name = zhat
    )
  }
  if (isTRUE(normalize.Z)) {
    zhat_col <- paste0(Z, ".hat")

    ## Normalize fitted Z so residualized instrument has mean zero:
    ## mean(Z - Z.hat) = 0 within sample ?? margin cell.
    data[
      ,
      (zhat_col) := get(zhat_col) + mean(get(z_col) - get(zhat_col), na.rm = TRUE),
      by = c("sample", margins)
    ]

    ## Binary Z: keep propensities valid after recentering.
    data[
      z_is_linear != TRUE,
      (zhat_col) := pmin(pmax(get(zhat_col), aipw.clip), 1 - aipw.clip)
    ]
  }

  ##Optional: Estimate conditional variance of Z, if required for target="all" & linear Z
  Zvarhat <- NULL

  if (need_z_var_global) {
    i_zvar <- need_z_var_rows

    data[i_zvar, Zvar := (get(z_col) - get(zhat))^2]
    Zvarhat <- "Zvar.hat"

    if (!is.null(X)) {
      if (!parametric) {
        crossfit_hat(
          data,
          i = i_zvar,
          y_name = "Zvar",
          x_names = X,
          margins = margins,
          weight_name = weight,
          forest_opts = Zparameters,
          folds = foldname,
          mode = if ("Z" %in% crossfit && is.null(foldname)) "across" else "within"
        )
      } else {
        parametric_hat(
          fml_rf,
          i = i_zvar,
          data = data,
          outcome_name = "Zvar",
          margins = margins,
          weight_name = weight,
          fixest_opts = Zparameters
        )
      }
    } else {
      leave_cluster_out_mean(
        DT = data,
        y_name = "Zvar",
        cluster_name = cluster,
        by = c("sample", margins),
        out_name = Zvarhat
      )
    }

    data[i_zvar, (Zvarhat) := pmax(get(Zvarhat), 1e-6)]
  }


  ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
  if (any(condition %in% c("MW","AHS"))&Y.res==TRUE)  {
    if (!is.null(X)) {
    crossfit_hat(
      data,
      y_name = Y,
      x_names = X,
      folds=foldname,
      margins = margins,
      weight_name = weight,
      mode=ifelse(("Y" %in% crossfit & is.null(foldname)==TRUE),"across","within"),
      forest_opts = Yparameters
    )
    } else {
        leave_cluster_out_mean(
          DT = data,
          y_name = Y,
          cluster_name = cluster,
          by = c("sample", margins),
          out_name = paste0(Y,".hat"),
          weight_name = weight
        )
    }
    data[, paste0(y_col, ".res") := get(y_col) - get(paste0(y_col, ".hat"))]
    data[, c(y_col, paste0(y_col, ".hat")) := NULL]
    y_name_rhs <- paste0(Y, ".res")
  } else y_name_rhs <- Y

  ##Check for onesided noncompliance
  if (has_margin_conditions) {
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

    if (
      nrow(os_res$threshold[one_sided == FALSE]) == 0 &&
      !any(condition %in% c("KR", "MW"))
    ) {
      stop(
        "One-sided noncompliance for all margins of Z and D - ",
        "all specified conditions are trivially satisfied."
      )
    }
  }
  ##Check for residual variation in Z
  byvars <- c("sample", margins)

  # Common summaries
  data[, n_group := .N, by = byvars]
  data[, nZ := uniqueN(get(z_col), na.rm = TRUE), by = byvars]
  data[, sdZ := stats::sd(get(z_col), na.rm = TRUE), by = byvars]
  data[, sd_res := stats::sd(get(z_col) - get(zhat), na.rm = TRUE), by = byvars]

  # Initialize helper
  if (!("min_cell_Z" %in% names(data))) data[, min_cell_Z := NA_integer_]
  data[, bad_Z := FALSE]
  if (need_linear_Z) {
    # Continuous Z: require group size, variation in Z, and variation in residualized Z
    data[z_is_linear == TRUE, bad_Z := (
      n_group < min_n |
        is.na(sdZ) | sdZ == 0 |
        is.na(sd_res) | sd_res == 0
    ), by = byvars]
  }
  if (has_margin_conditions) {
    # Discrete Z: require at least 2 support points, enough observations in each Z cell,
    # and variation in residualized Z
    data[z_is_linear == FALSE, min_cell_Z := {
      ztab <- table(get(z_col), useNA = "no")
      if (length(ztab) == 0L) NA_integer_ else min(ztab)
    }, by = byvars]

    data[z_is_linear == FALSE, bad_Z := (
      n_group < min_n |
        nZ < 2 |
        is.na(min_cell_Z) | min_cell_Z < min_n |
        is.na(sd_res) | sd_res == 0
    ), by = byvars]
  }

  if (length(margins) == 0L) {
    # No margins: if any sample part fails, drop everything
    drop_all <- data[, any(bad_Z, na.rm = TRUE)]

    if (isTRUE(drop_all)) {
      stop("At least one sample part has no residual variation in Z.")
      data <- data[0]
    }

  } else {
    # Identify margin cells where at least one sample part fails
    bad_margins <- unique(
      data[bad_Z == TRUE, ..margins]
    )

    if (nrow(bad_margins) > 0L) {
      # Build a readable message
      if (length(margins) == 1L) {
        msg <- paste(bad_margins[[margins]], collapse = ", ")
        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part has no usable variation in Z: ",
          msg
        )
      } else {
        bad_labels <- apply(bad_margins, 1, function(r) {
          paste(paste(names(r), r, sep = "="), collapse = ", ")
        })
        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part has no usable variation in Z:\n",
          paste("  -", bad_labels, collapse = "\n")
        )
      }

      # Anti-join on margins only: removes both sample parts
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

  ## helper: all subsets of Y support
  all_subsets <- function(vals, min_size = 1L, max_size = length(vals)) {
    out <- vector("list", 0L)
    for (k in min_size:max_size) {
      out <- c(out, combn(vals, k, simplify = FALSE))
    }
    out
  }

  # --------------------------------------------------
  # Build one unified condition index
  #   Generic columns:
  #     condition, linear, equation, zmargin, dval, yval, Avals
  # --------------------------------------------------

  idx_blocks <- list()

  ## helper for simple / AHS
  make_linear_blocks <- function(cond_name) {
    blocks <- list()

    ## ordinary margins: binarized Z and binarized D
    if ("none" %in% linear) {
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

      tmp[, condition := cond_name]
      tmp[, linear := "none"]

      blocks[["none"]] <- tmp
    }

    ## linear in D, margin in Z
    if ("D" %in% linear) {
      if (K > 2L) {
        tmp <- data.table::data.table(zmargin = Zsup[-1L])
      } else {
        tmp <- data.table::data.table()
      }

      tmp[, condition := cond_name]
      tmp[, linear := "D"]

      blocks[["D"]] <- tmp
    }

    ## linear in Z, margin in D
    if ("Z" %in% linear) {
      if (J > 2L) {
        tmp <- data.table::data.table(dval = Dsup[-1L])
      } else {
        tmp <- data.table::data.table()
      }

      tmp[, condition := cond_name]
      tmp[, linear := "Z"]

      blocks[["Z"]] <- tmp
    }

    ## linear in both Z and D
    if ("DZ" %in% linear) {
      tmp <- data.table::data.table()

      tmp[, condition := cond_name]
      tmp[, linear := "DZ"]

      blocks[["DZ"]] <- tmp
    }

    data.table::rbindlist(blocks, use.names = TRUE, fill = TRUE)
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

    idx_blocks[["MW"]] <- tmp
  }

  # KR
  if ("KR" %in% condition) {
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

    tmp[, condition := "KR"]

    A_specs <- list()

    for (dv in Dsup) {
      ## admissible A sets at each exact treatment value:
      ## bottom d: only singleton sets
      ## interior d: all nonempty subsets
      ## top d: all nonempty proper subsets
      if (joint == FALSE) {
        A_list <- all_subsets(Ysup, min_size = 1L, max_size = 1L)
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

  has_linear_col <- "linear" %in% names(data)

  linear_requested__ <- unique(as.character(linear))
  if (is.null(linear_requested__) || length(linear_requested__) == 0L) {
    linear_requested__ <- "none"
  }

  if (!is.null(Dcol) && Dcol %in% names(data)) {
    Dvals__ <- sort(unique(stats::na.omit(data[[Dcol]])))
    D_is_binary__ <- length(Dvals__) <= 2L && all(Dvals__ %in% c(0, 1))
  } else {
    D_is_binary__ <- FALSE
  }


  # ---------------- simple / AHS: linear-D Q ----------------
  # If there is a linear column, only rows with linear D/DZ get raw D.
  # If there is no linear column, the whole simple/AHS block follows the requested linear argument.

  linear_D_rows <- if (has_linear_col) {
    is_linear_condition & data[["linear"]] %in% c("D", "DZ")
  } else {
    is_linear_condition & any(linear_requested__ %in% c("D", "DZ"))
  }

  if (any(linear_D_rows)) {
    stopifnot(!is.null(Dname), Dname %in% names(data))

    data[
      linear_D_rows,
      Q := as.numeric(get(Dname))
    ]
  }


  # ---------------- simple / AHS: binarized-D Q ----------------
  # These are rows where Q should be 1(D >= dval).
  # If there is no linear column and linear = "D", this is FALSE, as desired.

  binarized_rows <- if (has_linear_col) {
    is_linear_condition & data[["linear"]] %in% c("none", "Z")
  } else {
    is_linear_condition & any(linear_requested__ %in% c("none", "Z"))
  }

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
        "ordinary simple/AHS rows or linear='none'/'Z' rows.",
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
      Q := as.numeric(get(Dcol)) -
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

  time=rbind(time,"Stack data across margins"=proc.time())

  ##Estimate Q.hat in stacked data
  if (!is.null(X)) {
    if (any(condition %in% c("simple", "KR"))) {
      i <- if (length(condition) > 1L) {
        which(data$condition %in% c("simple", "KR"))
      } else {
        NULL
      }

      if (!parametric) {
        crossfit_hat(
          data,
          i = i,
          y_name = "Q",
          x_names = X,
          folds = foldname,
          margins = margins,
          weight_name = weight,
          mode = if ("Q" %in% crossfit && is.null(foldname)) "across" else "within",
          forest_opts = Qparameters
        )
      } else {
        parametric_hat(
          fml_rf,
          data = data,
          outcome_name = "Q",
          i = i,
          margins = margins,
          weight_name = weight,
          fixest_opts = Qparameters
        )
      }
    }
  } else {
    i <- which(data$condition != "AHS")

    if (length(i)) {
      leave_cluster_out_mean(
        DT = data,
        i = i,
        y_name = "Q",
        cluster_name = cluster,
        by = c("sample", margins),
        out_name = "Q.hat",
        weight_name = weight
      )
    }
  }

  if ("AHS" %in% condition) {
    i <- if (length(condition) > 1L) {
      which(data$condition == "AHS")
    } else {
      NULL
    }

    if (!parametric) {
      crossfit_hat(
        data,
        i = i,
        y_name = "Q",
        folds = foldname,
        x_names = null_if_empty(c(X, y_name_rhs)),
        margins = margins,
        weight_name = weight,
        mode = if ("Q" %in% crossfit && is.null(foldname)) "across" else "within",
        forest_opts = Qparameters
      )
    } else {
      parametric_hat(
        append_to_main_rhs(fml_rf, y_name_rhs),
        data = data,
        outcome_name = "Q",
        i = i,
        margins = margins,
        weight_name = weight,
        fixest_opts = Qparameters
      )
    }
  }

  ## Check for usable (residual) variation in Q
  byvars <- c("sample", margins)

  # common summaries
  data[, n_group := .N, by = byvars]
  data[, nQ      := uniqueN(Q, na.rm = TRUE), by = byvars]
  data[, sdQ     := stats::sd(Q, na.rm = TRUE), by = byvars]

  # family/type flags
  data[, Q_continuous :=
         condition == "MW" |
         (
           condition %in% c("simple", "AHS") &
             "linear" %in% names(data) &
             linear %in% c("D", "DZ")
         )]
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
        x_names = null_if_empty(c(X, y_name_rhs)),
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
  ## FWL / residual scores for all linear-Z rows
  ## Q.hat and Z.hat already exist.
  ## ------------------------------------------------------------
  if (identical(target, "overlap") && any(data$z_is_linear == TRUE)) {
    i_fwl <- which(data$z_is_linear == TRUE)

    make_scores(
      DT = data,
      y_name = "Q",
      z_name = Z,
      y_hat_name = "Q.hat",
      z_hat_name = zhat,
      tau_name = NULL,
      target = "overlap",
      zvar_name = NULL,
      score_name = "scores",
      i = i_fwl,
      z_is_linear_name = "z_is_linear",
      clip = aipw.clip,
      var_floor = 1e-6
    )
    if (testtype=="forest") {
      fit_models(
        DT = data,
        i = i_fwl,
        forest_type = "regression",
        y_name = "scores",
        x_names = X,
        w_name = NULL,
        zvar_name = NULL,
        folds = foldname,
        margins = margins,
        weight_name = weight,
        cluster_name = cluster,
        forest_opts = Cparameters,
        aipw.clip = aipw.clip,
        shrink = (shrink > 0),
        target = "overlap"
      )

    }
  }


  ## ------------------------------------------------------------
  ## AIPW scores for non-linear / binary-Z rows + linear Z rows when target="all"
  ## ------------------------------------------------------------
  if (any(data$condition %in% c("simple", "KR"))) {
    i_simple_kr <- which(
      data$condition %in% c("simple", "KR") &
        (
          data$z_is_linear != TRUE |
            identical(target, "all")
        )
    )

    if (length(i_simple_kr)) {
      fit_models(
        DT = data,
        i = i_simple_kr,
        forest_type = "causal",
        y_name = "Q",
        x_names = X,
        w_name = Z,
        zvar_name = Zvarhat,
        folds = foldname,
        margins = margins,
        weight_name = weight,
        cluster_name = cluster,
        forest_opts = Cparameters,
        aipw.clip = aipw.clip,
        shrink = (shrink > 0),
        target = target
      )
    }
  }


  ## ------------------------------------------------------------
  ## AHS: AIPW scores for non-linear / binary-Z rows,
  ## with RHS outcome controls included in X
  ## ------------------------------------------------------------
  if (any(data$condition == "AHS")) {
    i_ahs <- which(
      data$condition == "AHS" &
        (
          data$z_is_linear != TRUE |
            identical(target, "all")
        )
    )

    if (length(i_ahs)) {
      fit_models(
        DT = data,
        i = i_ahs,
        forest_type = "causal",
        y_name = "Q",
        x_names = null_if_empty(c(X, y_name_rhs)),
        w_name = Z,
        zvar_name = Zvarhat,
        folds = foldname,
        margins = margins,
        weight_name = weight,
        cluster_name = cluster,
        forest_opts = Cparameters,
        aipw.clip = aipw.clip,
        shrink = (shrink > 0),
        target = target
      )
    }
  }

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

  if ("forest" == testtype) res=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,x_names=X,pool=poolmargins,select=selectmargins,gridpoints=gridpoints,margins=margins,screen=screen,alpha=alpha)
  if ("CART" == testtype) res=CART_test(data, x_names=X,margins=margins,weight=weight,cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,screen=screen,cluster=cluster,select=selectmargins,rpart_options=Rparameters)


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
  out <- c(res, list(
    time = time,
    obs = obs,
    margins = margin_index[, setdiff(names(margin_index), "Avals"), with = FALSE]
  ))
  return(out)
}
