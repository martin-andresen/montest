#' Monotonicity and LATE assumptions tests using sample splitting and machine learning
#'
#' \code{montest()} searches for violations of monotonicity and LATE assumptions in
#' data-adaptive subsets of the sample. It combines sample splitting,
#' cross-fitting, and generalized random forest (or CART-based subset search) to identify
#' regions of the covariate space or margins of the instrument or treatment where test
#' statistics are most negative, and then evaluates those regions in the held-out sample.
#' The function supports several testable conditions, including a simple test of whether
#' the first stage is negative \code{"simple"}, the Balke and Pearl (1997) condition \code{"BP"},
#' the Mourifie and Wan (2017) conditions  \code{"MW"}, and the first stage conditional on Y-test from
#' Andresen, Huber and Sloczynski (2026) \code{"AHS"}. Multivalued instruments, treatments and outcomes
#' are expanded across margins and discretized into bins before estimation. See Details.
#'
#' @param data A \code{data.frame} or \code{data.table} containing the analysis sample.
#'   Observations with missing values in any variables used by the call are dropped.
#' @param D Character scalar giving the name of the treatment variable.
#' @param Z Character scalar giving the name of the instrumental variable. \bold{Z should be coded so that higher values weakly increases treatment.}
#' @param X Optional character vector of covariate names used for subset discovery and
#'   nuisance estimation.
#' @param Y Optional character scalar for the outcome variable. Required for tests other
#'   than \code{"simple"}.
#' @param condition Character vector selecting which tests to run. Allowed values are any combination of
#'   \code{"simple"}, \code{"KR"} (Kwan-Roth conditions), \code{"MW"} (Mourifi?? and Wan conditions), \code{"AHS"} (Andresen-Huber-Sloczynski), or \code{"all"}.
#'   If \code{Y} is omitted, only \code{"simple"} is allowed.
#' @param inner.folds Optional integer giving the number of within-sample folds used for
#'   cross-fitting nuisance functions and, optionally, forest predictions. Set to
#'   \code{NULL} to disable the inner split. Defaults to NULL - nuissances and predictions from causal forests are fit out-of-bag. See option crossfit, which decides which parts this applies to.
#' @param crossfit Character vector of what parts of the procedure to cross-fit. Accepts "Z","Q","Y","C". If e.g. "Z" appears in crossfit, nuissances for Z are cross fit, either across outer sample part (if inner.folds==NULL), or within outer sample part across inner folds. If "Z" does not appear, OOB predictions are used. "C" is for the causal forest fit.
#' @param normalize.Z Logical, default TRUE; if \code{TRUE}, estimated instrument propensity scores are
#'   normalized after estimation.
#' @param aipw.clip Positive scalar in \code{(0,1)}, dfeault 1e-3, used to trim estimated propensity
#'   scores when augmented inverse-probability weighted scores are constructed.
#' @param weight Optional character scalar naming a nonnegative weight variable.
#' @param cluster Optional character scalar naming a cluster identifier. Cluster-robust
#'   inference is used in forest-based testing. CART testing canot be combined with cluster.
#' @param num.trees Integer, default 2000 giving the number of trees for the main forest fits.
#' @param seed Integer random seed, default 10101, set to NULL to disable setting seed
#' @param minsize Integer minimum effective sample size or minimum cluster count required
#'   for subset search and testing. Default 50.
#' @param shrink Shrink predicted treatment effects using empirical bayes before sorting. Default 0: No shrinkage. 1: Full shrinkage
#' @param gridtypeY,gridtypeD,gridtypeZ Character strings controlling how continuous
#'   variables are discretized before stacking. Must be one of \code{"equisized"} or
#'   \code{"equidistant"}.
#' @param sim Logical; for development and testing
#' @param Ysubsets,Dsubsets,Zsubsets Integers giving the number of bins used when
#'   discretizing \code{Y}, \code{D}, and \code{Z}, respectively. Dsubsets and Zsubsets may be set to 0L for linear models.
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
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. No margin can appear in both pool and select. Relevant margins that appear in neither are all tested, and tests are corrected for multiple hypothesis testing.
#' @param select Character vector controlling which dimensions are selected over when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dval"},
#'   \code{"yval"}, \code{"condition"}, \code{"equation"},
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. No margin can appear in both pool and select. Relevant margins that appear in neither are all tested, and tests are corrected for multiple testing.
#' @param cp,maxrankcp,alpha,prune,preselect Tuning parameters for the CART-based search
#'   routine. See Details.
#' @param Zparameters,Yparameters,Qparameters,Dparameters,Cparameters Named lists of
#'   additional arguments passed to the underlying GRF estimation routines for different
#'   nuisance or target models.See regression_forest and causal_forest for details.
#'
#' @details
#' The procedure works in several stages:
#'
#' \enumerate{
#'   \item The data are restricted to the variables used in the call, complete cases are
#'   kept, and the sample is split into two halves (optionally respecting clusters).
#'   \item Continuous or multivalued instruments, treatments and outcomes discretized into bins.
#'   \item The data is stacked across instrument margins, treatment margins, outcomes,
#'   equations, and test conditions.
#'   \item Nuisance functions such as \code{Z.hat}, \code{Q.hat} are estimated and outcomes are residualized
#'   using regression forests within margins and and cross-fitted within sample halves.
#'   \item Separate causal forests of the outcome \code{Q}, which depends on the condition being tested, on
#'   the instrument \code{Z} using features \code{X} (and optionally \code{Y} for MW and AHS conditions),
#'   treatment effects are predicted in and out of sample and scores constructed
#'   \item Each sample part (optionally within margins, depending on the options in code{pool}) is sorted
#'   according to treatment effects, and the mean of scores is estimated numerically for all possible cutoffs
#'   in predicted treatment effects. Select the cutoff with the smallest t-statistic on the mean of scores
#'   Alternatively, subset selection can be done using a CART algorithm.
#'   \item Depending on the choices in \code{select}, promising subsets are evaluated in the opposite sample half. If performing multiple tests,
#'   depending on the option \code{pool}, p-values are adjusted for multiple testing.
#' }
#'
#' \code{montest} supports weights and clustering, and allow for multivalued treatments and instruments
#' by binarizing instruments and treatments into quantile or equisized bins. The command also tests for
#' one-sided monotonicity (within margins of the treatment and instrument) and if found, warns the user
#' and skips testing any trivially satisified conditions.
#'
#'
#' @return
#' A named list:
#'
#' \describe{
#'   \item{\code{results}}{A Matrix of train and test results by sample and margin cell.}
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
#' @section CART tuning parameters:
#' When \code{testtype = "CART"}, subset search is based on regression trees fit to the
#' score variable. The argument \code{cp} controls tree complexity, \code{maxrankcp}
#' limits the pruning table rank considered, \code{alpha} enters the preselection rule,
#' \code{prune} toggles pruning, and \code{preselect} determines which candidate leaves
#' are retained before testing. Supported preselection rules are \code{"none"},
#' \code{"minimum"}, \code{"negative"}, and \code{"nonpositive"}.
#'
#' @section Reserved names:
#' The function uses internal variable names such as \code{Q}, \code{*.hat},
#' \code{sample}, \code{condition}, \code{equation}, \code{ybin}, \code{zmargin}, and
#' \code{dmargin}. Some of these are created internally, and names conflicting with the
#' \code{*.hat} conventions are rejected.
#'
#' @examples
#' \dontrun{
#'
#' #Generate data - Simulation DGP from Farbmacher et. al.
#' data=fct_datasim(setup="A",dgp=2,n=3000)
#'
#' # Simple monotonicity-style test
#' out <- montest(
#'   data = data,
#'   D = "D",
#'   Z = "Z",
#'   X = c("Xvar1", "Xvar2", "Xvar3"),
#'   condition = "simple",
#'   testtype = "forest")
#'
#' # Test multiple conditions, pooling evidence
#' out2 <- montest(
#'   data = data,
#'   D = "D",
#'   Z = "Z",
#'   Y="Y",
#'   X = c("Xvar1", "Xvar2", "Xvar3"),
#'   condition = c("simple","KR","MW"))
#' }
#'
#' @seealso montestplot LATEtest
#' @export

montest=function(data,D,Z,X=NULL,Y=NULL,condition=NULL,inner.folds=NULL,crossfit=NULL,
                 normalize.Z=TRUE,aipw.clip=0,weight=NULL,cluster=NULL,seed=10101,minsize=50L,
                 gridtypeY=NULL,gridtypeD=NULL,gridtypeZ=NULL,stratify=NULL,joint=TRUE,
                 Ysubsets = 4L, Dsubsets = 4L,Zsubsets=4L,Y.res=TRUE,testtype="forest",
                 gridpoints=NULL,min_n=1L,pool="all",select="none",shrink=0, ##forest opts
                 cp=0,maxrankcp=10L,rpart_options=NULL,alpha=0.05,prune=TRUE,preselect="fgk_relevant", ##CART opts
                 Zparameters=list(),Yparameters=list(),Qparameters=list(),Dparameters=list(),Cparameters=list()
                 #tune.Qparameters="none",tune.Zparameters="none",tune.Cparameters="none",tune.Yparameters="none",tune.Dparameters="none",
                 #tune.num.trees=200,tune.num.reps=50,tune.num.draws=1000,tunetype="one" ##tuning options
){

  time=proc.time()
  if (is.null(seed)==FALSE) set.seed(seed)


  ################### 1 CHECK INPUT #####################
  gridtypeY=match.arg(gridtypeY,c("equidistant","equisized"))
  gridtypeD=match.arg(gridtypeD,c("equidistant","equisized"))
  gridtypeZ=match.arg(gridtypeZ,c("equidistant","equisized"))
  if (is.null(crossfit)==FALSE) match.arg(crossfit,c("Z","Q","C","Y"),several.ok=TRUE)
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

  if (is.integer(minsize)==FALSE|minsize<=0) stop("Minsize must be an integer >0.")

  if (is.null(condition)==TRUE) {
    if (is.null(Y)==TRUE) {
      condition="simple"
    } else condition="all"
  }

  condition=match.arg(condition,c("simple","KR","MW","AHS","all"),several.ok=TRUE) ##took out support for Kitagawa
  if (length(Y)>1) stop("More than one outcome currently not supported.")
  if ("all" %in% condition) {
    condition=c("simple","KR","MW","AHS")
  }
  if (is.null(Y)==TRUE) {
    if (sum(!condition %in% "simple")>0) {
      stop("Other conditions than simple may not be used when Y is not specified. Specify the Y argument or use condition=simple")
    }
  }

  if (Ysubsets<=1|is.integer(Ysubsets)==FALSE) stop("Ysubsets must be an integers larger than 1")
  if (Dsubsets<0|Dsubsets==1|is.integer(Dsubsets)==FALSE) stop("Dsubsets must be an integer equal to 0 (linear model) or larger than 1. Use #L notation for integer.")
  if (Zsubsets<0|Zsubsets==1|is.integer(Zsubsets)==FALSE) stop("Zsubsets must be an integer equal to 0 (linear model) or larger than 1. Use #L notation for integer.")

  if (Dsubsets==0&"KR" %in% condition) stop("Dsubsets=0 (linear model) incompatible with condition=KR, use only for condition=simple or condition=AHS.")
  if (Zsubsets==0&"KR" %in% condition) stop("Zsubsets=0 (linear model) incompatible with condition=KR, use only for condition=simple or condition=AHS.")

  if ((sum(pool=="none")==1)&(sum(pool=="all")==1)) stop("Do not specify both none and all in pool().")
  else if (sum(pool=="all")==1) pool=c("zmargin","dval","yval","condition","equation","sample")
  else if (sum(pool=="none")==1) pool=c()
  else if (is.null(pool)==FALSE) pool=match.arg(pool,c("zmargin","dval","yval","condition","sample"),several.ok=TRUE)

  if ((sum(select=="none")==1)&(sum(select=="all")==1)) stop("Do not specify both none and all in select().")
  else if (sum(select=="all")==1) select=c("zmargin","dval","yval","condition","equation","sample")
  else if (sum(select=="none")==1) select=c()
  else if (is.null(select)==FALSE) select=match.arg(select,c("zmargin","dval","yval","condition","equation","sample"),several.ok=TRUE)

  overlap <- intersect(pool, select)
  if (length(overlap) > 0L) {
    stop("Error: The following names appear in both pool and select: ",
         paste(overlap, collapse = ", "))
  }

  if (sum(sum(grepl("Z.hat",colnames(data))))) stop("Variable name beginning with Z.hat discovered, reserved for internal use. Please rename.")
  if (sum(sum(grepl("D.hat",colnames(data))))) stop("Variable name beginning with D.hat discovered, reserved for internal use. Please rename.")
  if (sum(sum(grepl("Q.hat",colnames(data))))) stop("Variable name beginning with Q.hat discovered, reserved for internal use. Please rename.")
  if ("Q" %in% colnames(data)) stop("Data contains variable named Q, which is reserved for internal use. Please rename.")

  gridtypeZ=match.arg(gridtypeZ,c("equidistant","equisized"))
  gridtypeY=match.arg(gridtypeY,c("equidistant","equisized"))
  gridtypeD=match.arg(gridtypeD,c("equidistant","equisized"))

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

  if (sum(condition %in% c("AHS","MW","KR"))==0&is.null(Y)==FALSE) Y=NULL
  time=rbind("Check input"=proc.time())

  ###################### 2 Prepare data #########################3
  XW=X

  if (sum(c("MW","AHS") %in% condition)>0) {XWY=c(X,Y)} else {XWY=XW}

  data=data.table(data)
  allvars=c(X,Y,D,Z,weight,cluster)
  data=data[,..allvars]
  dropped=sum(!complete.cases(data))
  if (dropped>0) message(paste("Note: dropped ",dropped," observations with missing data on one or more input variables."))
  data=data[complete.cases(data)]
  n=nrow(data)
  if (is.null(cluster)==FALSE) {
    G=length(unique(data[,get(..cluster)]))
    if (G<=2*minsize) stop("Number of clusters is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
  } else {
    if (n<=2*minsize) stop("Number of observations is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
    G=NA
    }
  obs=c(N=n,G=G)
  data[,id_:=(1:n)]


  ##OUTER SPLIT
  if (is.null(stratify)==TRUE) {
    if (is.null(cluster)==TRUE&Zsubsets>0) strat=Z else strat=NULL
  }
  make_group_folds(data,K = 2,cluster_name = cluster, fold_col = "sample",verbose = FALSE,diag_prefix=NULL,strat_col=strat)

  ##OPTIONAL INNER SPLIT
  if (is.null(inner.folds)==FALSE) {
    make_group_folds(data,K = inner.folds,cluster_name = cluster,fold_col = "cf_fold",verbose = FALSE,by_col="sample",diag_prefix=NULL,strat_col=stratify)
    foldname="cf_fold"
  } else {
    foldname=NULL
  }

  if (is.null(cluster)==TRUE) cluster="id_"

  time=rbind(time,"Prepare data"=proc.time())

  ############### 3 Discretize Z, D and Y into subsets ###############


  if (is.null(Dsubsets)==FALSE&Dsubsets>0) { ##bin treatment
    data <- binarize_var(
      data    = data,
      var     = D,
      ngroups = Dsubsets,
      gridtype = gridtypeD,
      wvar    = wvar
    )
  }

  if (is.null(Zsubsets)==FALSE&Zsubsets>0) { ##bin instrument
    data <- binarize_var(
      data    = data,
      var     = Z,
      ngroups = Zsubsets,
      gridtype = gridtypeZ,
      wvar    = wvar
    )
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
    ybincol=paste0(Y,".bin")
    Ysup=sort(unique(data[,get(..ybincol)]));L=length(Ysup)
  }

  n=nrow(data);
  Dsup=sort(unique(data[,get(..D)]));J=length(Dsup)
  Zsup=sort(unique(data[,get(..Z)]));K=length(Zsup)

  if (J>2&("MW" %in% condition)) stop("Multivalued treatment not supported with condition MW.")
  if (J==2&K==2&is.null(X)==TRUE&all(condition=="simple")) {
    stop("Nothing to test with a binary treatment, a binary instrument, the simple first stage condition and no variables in X.")
  }

  time=rbind(time,"Binarize treatment, instrument and outcome"=proc.time())

  ######################## 6a STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
  margins=c()

  ##Stack for margins of Z
  if (K>2) { ##Expand to all margins of Z
    minZ=min(data[,..Z]);maxZ=max(data[,..Z]);zvals=sort(unique(data[,get(..Z)]))
    times = 1L + 1L * as.integer(data[[Z]] > minZ & data[[Z]] < maxZ)
    data = data[rep.int(seq_len(nrow(data)), times)]

    zname=Z
    data[get(zname)==minZ,zmargin:=zvals[2]]
    data[get(zname)==maxZ,zmargin:=zvals[length(zvals)]]
    data[is.na(zmargin),
         zmargin := {
           z0 <- .SD[[1L]][1L]
           k  <- match(z0, zvals)
           stopifnot(.N == 2L, !is.na(k), k < length(zvals))
           c(zvals[k], zvals[k + 1L])
         },
         by = id_,
         .SDcols = zname]

    data[,(Z):=as.integer(get(..Z)>=zmargin)]
    margins=c(margins,"zmargin")
  } else if (Zsubsets>0) { ##MAKE SURE Z is dummy!
    data[, (Z) := as.integer(get(..Z) == max(get(..Z)))]
  }


  ##estimate Z.hat for each margin in stacked data
  if (is.null(X)==FALSE) {
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
  } else { ##leave-cluster out mean
    zhat <- paste0(Z, ".hat")
    by_cl <- c("sample", margins, cluster)

    # cluster totals within each sample x margins x cluster
    cl <- data[, .(
      z_sum = sum(get(Z)),
      n_cl  = .N
    ), by = by_cl]

    # totals within each sample x margins cell
    tot <- cl[, .(
      z_tot = sum(z_sum),
      n_tot = sum(n_cl),
      G     = .N
    ), by = c("sample", margins)]

    # leave-cluster-out mean
    cl <- tot[cl, on = c("sample", margins)]
    cl[, (zhat) := fifelse(
      n_tot > n_cl,
      (z_tot - z_sum) / (n_tot - n_cl),
      NA_real_
    )]

    # merge back to observation level
    data <- cl[, c(by_cl, zhat), with = FALSE][data, on = by_cl]
  }
  if (normalize.Z==TRUE&Zsubsets>0) { #Normalize propensity scores
    data[, (paste0(Z, ".hat")) :=
           get(paste0(..Z, ".hat")) * .N / sum(get(..Z) / get(paste0(..Z, ".hat"))),
         by = c("sample", margins)]
  }


  ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
  if ((sum(condition %in% c("MW","AHS"))>0)&Y.res==TRUE)  {
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
    data[,paste0(Y,".res"):=get(..Y)-get(paste0(..Y,".hat"))]
    data[, c(Y, paste0(Y, ".hat")) := NULL]
  }

  ##Check for onesided noncompliance
  if (Dsubsets>0&Zsubsets>0) {
    os_res <- test_one_sided_noncompliance(data = data, D = D, Z = Z,zmargin_var=margins)
    osmargins=margins
    if (nrow(os_res$threshold[one_sided==TRUE])>0) os_res
    if (nrow(os_res$threshold[one_sided==FALSE])==0&sum(condition %in% c("KR","MW"))==0) {
      stop("One-sided noncompliance for all margins of Z and D - all specified conditions are trivially satisfied.")
    }
  }

  ##Check for residual variation in Z
  byvars <- c("sample", margins)
  Zhat <- paste0(Z, ".hat")
  Zcol <- Z

  # Common summaries
  data[, n_group := .N, by = byvars]
  data[, nZ      := uniqueN(get(Zcol), na.rm = TRUE), by = byvars]
  data[, sdZ     := stats::sd(get(Zcol), na.rm = TRUE), by = byvars]
  data[, sd_res  := stats::sd(get(Zcol) - get(Zhat), na.rm = TRUE), by = byvars]

  # Initialize helper
  if (!("min_cell_Z" %in% names(data))) data[, min_cell_Z := NA_integer_]

  if (Zsubsets == 0) {
    # Continuous Z: require group size, variation in Z, and variation in residualized Z
    data[, bad_Z := (
      n_group < min_n |
        is.na(sdZ) | sdZ == 0 |
        is.na(sd_res) | sd_res == 0
    ), by = "sample"]

  } else {
    # Discrete Z: require at least 2 support points, enough observations in each Z cell,
    # and variation in residualized Z
    data[, min_cell_Z := {
      ztab <- table(get(Zcol), useNA = "no")
      if (length(ztab) == 0L) NA_integer_ else min(ztab)
    }, by = byvars]

    data[, bad_Z := (
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
  ##CONSTRUCT INDEX MATRIX OF MARGINS
  # helper: all subsets of Y support
  all_subsets <- function(vals, min_size = 1L, max_size = length(vals)) {
    out <- vector("list", 0L)
    for (k in min_size:max_size) {
      out <- c(out, combn(vals, k, simplify = FALSE))
    }
    out
  }

  # --------------------------------------------------
  # Build one unified condition index
  #    Generic columns:
  #      condition, equation, zmargin, dval, yval
  # --------------------------------------------------

  idx_blocks <- list()

  # simple
  if ("simple" %in% condition) {
    tmp <- os_res$threshold[one_sided == FALSE]
    keep <- intersect(c("zmargin", "dmargin"), names(tmp))
    tmp <- tmp[, ..keep]
    if ("dmargin" %in% names(tmp)) setnames(tmp, "dmargin", "dval")
    tmp[, condition := "simple"]
    idx_blocks[["simple"]] <- tmp
  }

  # AHS
  if ("AHS" %in% condition) {
    tmp <- os_res$threshold[one_sided == FALSE]
    keep <- intersect(c("zmargin", "dmargin"), names(tmp))
    tmp <- tmp[, ..keep]
    if ("dmargin" %in% names(tmp)) setnames(tmp, "dmargin", "dval")
    tmp[, condition := "AHS"]
    idx_blocks[["AHS"]] <- tmp
  }
  # MW
  if ("MW" %in% condition) {
    tmp <- os_res$threshold[one_sided == FALSE]
    keep <- intersect(c("zmargin", "dmargin","direction"), names(tmp))
    tmp <- tmp[, ..keep]
    if ("dmargin" %in% names(tmp)) setnames(tmp, "dmargin", "dval")
    tmp[, condition := "MW"]
    eq_idx <- data.table(equation = 0:1, dummy__ = 1L)
    tmp[, dummy__ := 1L]
    tmp <- eq_idx[tmp, on = "dummy__", allow.cartesian = TRUE][, dummy__ := NULL]

    if (nrow(os_res$threshold[one_sided==TRUE])>0) {
      tmp <- tmp[!((direction == "Never-takers only"  & equation == 1L) | (direction == "Always-takers only" & equation == 0L))]
    }
    tmp[, direction := NULL]
    idx_blocks[["MW"]] <- tmp
  }
  # KR
  if ("KR" %in% condition) {
    tmp <- os_res$exact[trivial_exact == FALSE]
    keep <- intersect(c("zmargin", "dval"), names(tmp))
    tmp <- tmp[, ..keep]
    if ("dmargin" %in% names(tmp)) setnames(tmp, "dmargin", "dval")
    tmp[, condition := "KR"]

    A_specs <- list()

    for (dv in Dsup) {
      # admissible A sets at each exact treatment value:
      # bottom d: only singleton sets (joint adds nothing)
      # interior d: all nonempty subsets
      # top d: all nonempty proper subsets
      if (joint==FALSE) {
        A_list <- all_subsets(Ysup, min_size = 1L, max_size = 1L)
      } else if (dv == min(Dsup)) {
        A_list <- all_subsets(Ysup, min_size = 1L, max_size = 1L)
      } else if (dv == max(Dsup)) {
        A_list <- if (L >= 2L) all_subsets(Ysup, min_size = 1L, max_size = L - 1L) else list()
      } else {
        A_list <- all_subsets(Ysup, min_size = 1L, max_size = L)
      }

      if (length(A_list)) {
        A_specs[[as.character(dv)]] <- rbindlist(lapply(A_list, function(a) {
          lbl <- paste0("{", paste(a, collapse = ","), "}")
          data.table(
            dval = dv,
            yval = lbl,
            Avals = list(a)
          )
        }))
      }
    }

    A_idx <- rbindlist(A_specs, use.names = TRUE, fill = TRUE)
    tmp <- A_idx[tmp, on = "dval", allow.cartesian = TRUE]
    idx_blocks[["KR"]] <- tmp
  }

  margin_index <- rbindlist(idx_blocks, use.names = TRUE, fill = TRUE)
  margins=setdiff(names(margin_index),"Avals")

  # --------------------------------------------------
  # 2) Cross with zmargin-stacked data
  # --------------------------------------------------

  if (!"zmargin" %in% names(data)) {
    data[, zmargin := NA_integer_]
    margin_index[, zmargin := NA_integer_]
  }
  data <- margin_index[data, on = "zmargin", allow.cartesian = TRUE]
  if ("zmargin" %in% names(margin_index) && all(is.na(margin_index[["zmargin"]]))) {
    margin_index[, zmargin := NULL]
  }

  # --------------------------------------------------
  #  Create Q by condition
  # --------------------------------------------------

  Dcol   <- D
  Ycol   <- paste0(Y,".bin")
  Zcol   <- Z
  Zhat   <- paste0(Zcol, ".hat")

  data[, Q := NA_real_]


  if ("dval" %in% names(data)) {
    data[, D_bin := fifelse(
           is.na(dval),
           as.numeric(get(Dcol)),
           as.numeric(get(Dcol) >= dval)
         )]
  } else {
    data[,D_bin := as.numeric(get(Dcol))]
  }

  # simple / AHS
  if (any(c("AHS","simple") %in% condition)) {
    data[condition %in% c("simple","AHS"),Q:=D_bin]
  }
  # MW
  if ("MW" %in% condition) {
  data[condition == "MW",
       Q := equation * (
         (1 - get(Zhat)) * D_bin * get(Zcol) -
           get(Zhat) * D_bin * (1 - get(Zcol))
       ) +
         (1 - equation) * (
           get(Zhat) * (1 - D_bin) * (1 - get(Zcol)) -
             (1 - get(Zhat)) * (1 - D_bin) * get(Zcol)
         )]
  }


  # KR (cellwise + joint)
  if ("KR" %in% condition) {
    dmin <- min(data[[Dcol]], na.rm = TRUE)

    data[condition == "KR" & dval == dmin,
         Q := -as.integer(get(Ycol) %in% Avals[[1L]] & get(Dcol) == dval),
         by = .(dval, yval)]

    data[condition == "KR" & dval > dmin,
         Q := D_bin -
           as.integer(get(Ycol) %in% Avals[[1L]] & get(Dcol) == dval),
         by = .(dval, yval)]
  }

  time=rbind(time,"Stack data across margins"=proc.time())

  ##Estimate Q.hat in stacked data
  if (sum(condition %in% c("simple","KR"))>0) {
    if (length(condition)>1) i=which(data$condition %in%  c("simple","KR")) else i=NULL
    crossfit_hat(
      data,
      i=i,
      y_name = "Q",
      x_names = X,
      folds=foldname,
      margins = margins,
      weight_name = weight,
      mode=ifelse(("Q" %in% crossfit & is.null(foldname)==TRUE),"across","within"),
      forest_opts = Qparameters
    )
  }

  if (("AHS" %in% condition)>0) {
    if (length(condition)>1) i=which(data$condition=="AHS") else i=NULL
    yname_lhs= if (Y.res==TRUE) paste0(Y,".res") else Y
    crossfit_hat(
      data,
      i=i,
      y_name = "Q",
      folds=foldname,
      x_names = c(X,yname_lhs),
      margins = margins,
      weight_name = weight,
      mode=if ("Q" %in% crossfit & is.null(foldname)==TRUE) "across" else "within",
      forest_opts = Qparameters
    )
  }

  ##Check for usable (residual) variation in Q
  byvars <- c("sample", margins)

  # common summaries
  data[, n_group := .N, by = byvars]
  data[, nQ := uniqueN(Q, na.rm = TRUE), by = byvars]
  data[, sdQ := stats::sd(Q, na.rm = TRUE), by = byvars]

  # family/type flags
  data[, Q_continuous :=
         condition == "MW" |
         (condition %in% c("simple", "AHS") & Dsubsets == 0)]

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

  # simple / AHS continuous case
  data[condition %in% c("simple", "AHS") & Q_continuous == TRUE,
       bad_Q := (
         n_group < min_n |
           is.na(sd_resQ) | sd_resQ == 0
       ),
       by = byvars]

  # simple / AHS discrete case
  data[condition %in% c("simple", "AHS") & Q_continuous == FALSE,
       bad_Q := (
         n_group < min_n |
           nQ < 2 |
           is.na(min_cell_Q) | min_cell_Q < min_n |
           is.na(sd_resQ) | sd_resQ == 0
       ),
       by = byvars]

  # KR always discrete
  data[condition == "KR",
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
          " margin cell(s) because at least one sample part / condition cell has no usable variation in Q: ",
          msg
        )
      } else {
        bad_labels <- apply(bad_margins, 1, function(r) {
          paste(paste(names(r), r, sep = "="), collapse = ", ")
        })
        message(
          "Dropping ", nrow(bad_margins),
          " margin cell(s) because at least one sample part / condition cell has no usable variation in Q:\n",
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
  if (!"C" %in% crossfit) foldname=NULL #Do not crossfit causal forest, just the nuissances - use OOB for forest.

  if (sum(condition %in% c("simple","KR"))>0) {
    if (length(condition)>1) i=which(data$condition %in% c("simple","KR")) else i=NULL
    fit_models(data,
               i=i,
               forest_type = "causal",
               y_name="Q",
               x_names=X,
               w_name=Z,
               folds=foldname,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = Cparameters,
               aipw.clip=aipw.clip,
               shrink=(shrink>0))
  }


  if (sum(condition == "AHS")>0) {
    if (length(condition)>1) i=which(data$condition=="AHS") else i=NULL
    yname_lhs= if (Y.res==TRUE) paste0(Y,".res") else Y
    fit_models(data,
               forest_type = "causal",
               i=i,
               y_name="Q",
               x_names=c(X,yname_lhs),
               w_name=Z,
               folds=foldname,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = Cparameters,
               aipw.clip=aipw.clip,
               shrink=(shrink>0))
  }

  if (sum(condition == "MW")>0) {
    if (length(condition)>1) i=which(data$condition=="MW") else i=NULL
    yname_lhs= if (Y.res==TRUE) paste0(Y,".res") else Y
    fit_models(data,
               forest_type = "regression",
               i=i,
               y_name="Q",
               x_names=c(X,yname_lhs),
               folds=foldname,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = Cparameters,
               aipw.clip=aipw.clip,
               shrink=(shrink>0))
  }

  ###EMPIRICAL BAYES SHRINKAGE IF SHRINK>0 #######


  if (shrink>0) {
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

  ##res=list()
  ##if (sim==TRUE) {
  ##  poollist=list(character(0),margins[!margins %in% "condition"],c(margins[!margins %in% "condition"],"sample"),c(margins,"sample"),margins)
  ##  treelist=c("CART","forest")
  ##
  ##  for (p in 1:5) {
  ##    res[[paste0(c("forest",p),collapse="_")]]=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,x_names=X,pool=poollist[[p]],gridpoints=gridpoints,margins=margins)
  ##    res[[paste0(c("CART",p),collapse="_")]]=CART_test(data, x_names=X,margins=margins,weight=weight,cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,preselect=preselect,cluster=cluster,pool=poollist[[p]])
  ##  }
  ##} else {   ##}

  if ("forest" == testtype) res=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,x_names=X,pool=poolmargins,select=selectmargins,gridpoints=gridpoints,margins=margins,preselect=preselect,alpha=alpha)
  if ("CART" == testtype) res=CART_test(data, x_names=X,margins=margins,weight=weight,cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,preselect=preselect,cluster=cluster,select=selectmargins,pool=poolmargins,rpart_options=rpart_options)


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
  return=c(res,list(time=time,obs=obs,margins=margin_index[, setdiff(names(margin_index), "Avals"), with = FALSE]))
  return(return)
}
