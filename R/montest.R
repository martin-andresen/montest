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
#' @param Z Character scalar giving the name of the instrumental variable.
#' @param X Optional character vector of covariate names used for subset discovery and
#'   nuisance estimation.
#' @param Y Optional character vector of outcome variable names. Required for tests other
#'   than \code{"simple"}.
#' @param test Character vector selecting which tests to run. Allowed values are any combination of
#'   \code{"simple"}, \code{"BP"}, \code{"MW"}, \code{"AHS"}, or \code{"all"}.
#'   If \code{Y} is omitted, only \code{"simple"} is allowed.
#' @param inner.folds Optional integer giving the number of within-sample folds used for
#'   cross-fitting nuisance functions and, optionally, forest predictions. Set to
#'   \code{NULL} to disable the inner split. Defaults to 5.
#' @param crossfit Character vector of what parts of the procedure to cross-fit. Accepts "Z","Q,"Y","C","D". If e.g. "Z" appears in crossfit, nuissances for Z are cross fit, either across outer sample part (if inner.folds==NULL), or within outer sample part across inner folds. If "Z" does not appear, OOB predictions are used. "C" is for the causal forest fit.
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
#'   discretizing \code{Y}, \code{D}, and \code{Z}, respectively.
#' @param Y.res Logical; if \code{TRUE}, outcomes are residualized from X before tests that use
#'   outcome on the right hand side \code{MW,AHS}.
#' @param testtype Character string selecting the subset-search routine. Must be
#'   \code{"forest"} or \code{"CART"}. Default: Forest.
#' @param gridpoints Optional integer controlling the number of candidate cutoffs searched
#'   by the forest-based test. If \code{NULL}, all eligible cutoffs are considered.
#' @param min_n Integer minimum number of treated and untreated instrument observations
#'   required within each sample half and margin cell.
#' @param pool Character vector controlling which dimensions are pooled when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dmargin"},
#'   \code{"ybin"}, \code{"condition"}, \code{"equation"}, \code{"outcome"},
#'   \code{"sample"}, \code{"all"}, and \code{"none"}. No margin can appear in both pool and select. Relevant margins that appear in neither are all tested, and tests are corrected for multiple testing.
#' @param select Character vector controlling which dimensions are selected over when finding testing subsets
#'   testing subsets. Allowed values are \code{"zmargin"}, \code{"dmargin"},
#'   \code{"ybin"}, \code{"condition"}, \code{"equation"}, \code{"outcome"},
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
#'   in predicted treatment effects. Select the cutoff with the smallest t-statistic on themean of scores
#'   Alterantively, subset selection can be done using a CART algorithm.
#'   \item The corresponding subset is evaluated in the opposite sample half. If performing multiple tests,
#'   depending on the option  \cite{pool}, p-values are optionally adjusted for multiple testing.
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
#'   \item{\code{results}}{A Matrix of ain train/test results by sample and margin cell.}
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
#'   test = "simple",
#'   testtype = "forest")
#'
#' # Test multiple conditions, pooling evidence
#' out2 <- montest(
#'   data = data,
#'   D = "D",
#'   Z = "Z",
#'   Y="Y",
#'   X = c("Xvar1", "Xvar2", "Xvar3"),
#'   test = c("simple","BP","MW"))
#' }
#'
#' @seealso montestplot LATEtest
#' @export

montest=function(data,D,Z,X=NULL,Y=NULL,test=NULL,inner.folds=5,crossfit=c("Z","Q","forest","Y"),
                 normalize.Z=TRUE,aipw.clip=1e-3,weight=NULL,cluster=NULL,num.trees=2000,seed=10101,minsize=50,
                 gridtypeY="equisized",gridtypeD="equisized",gridtypeZ="equisized",
                 Ysubsets = 4, Dsubsets = 4,Zsubsets=4,Y.res=TRUE,testtype="forest",
                 gridpoints=NULL,min_n=1L,pool="all",select="none",shrink=0, ##forest opts
                 cp=0,maxrankcp=10L,rpart_options=NULL,alpha=0.05,prune=TRUE,preselect="fgk_relevant", ##CART opts
                 Zparameters=list(),Yparameters=list(),Qparameters=list(),Dparameters=list(),Cparameters=list(),
                 tune.Qparameters="none",tune.Zparameters="none",tune.Cparameters="none",tune.Yparameters="none",tune.Dparameters="none",
                 tune.num.trees=200,tune.num.reps=50,tune.num.draws=1000,tunetype="one" ##tuning options
){

  time=proc.time()
  if (is.null(seed)==FALSE) set.seed(seed)


  ################### 1 CHECK INPUT #####################
  crossfit=match.arg(crossfit,c("Z","Q","C","Y"),several.ok=TRUE)
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

  if ((minsize!=floor(minsize))|minsize<0) stop("Minsize must be an integer >0.")

  if (is.null(test)==TRUE) {
    if (is.null(Y)==TRUE) {
      test="simple"
    } else test="all"
  }


  test=match.arg(test,c("simple","BP","MW","AHS","all"),several.ok=TRUE) ##took out support for Kitagawa
  if ("all" %in% test) {
    test=c("simple","BP","MW","AHS")
  }
  if (length(test)==1) {condition=test} else condition=NULL

  if (is.null(Y)==TRUE) {
    if (sum(!test %in% "simple")>0) {
      stop("Other tests than simple may not be used when Y is not specified. Specify the Y argument or use test=simple")
    }
  }

  if ((Ysubsets<=1)|(Dsubsets<=1)|(Zsubsets<=1)) stop("Ysubsets, Dsubsets and Zsubsets must be integers larger than 1")

  if ((sum(pool=="none")==1)&(sum(pool=="all")==1)) stop("Do not specify both none and all in pool().")
  else if (sum(pool=="all")==1) pool=c("zmargin","dmargin","ybin","condition","equation","outcome","sample")
  else if (sum(pool=="none")==1) pool=c()
  else if (is.null(pool)==FALSE) pool=match.arg(pool,c("zmargin","dmargin","ybin","condition","equation","outcome","sample"),several.ok=TRUE)

  if ((sum(select=="none")==1)&(sum(select=="all")==1)) stop("Do not specify both none and all in select().")
  else if (sum(select=="all")==1) select=c("zmargin","dmargin","ybin","condition","equation","outcome","sample")
  else if (sum(select=="none")==1) select=c()
  else if (is.null(select)==FALSE) select=match.arg(select,c("zmargin","dmargin","ybin","condition","equation","outcome","sample"),several.ok=TRUE)
  #if (testtype=="CART") pool=c()

  overlap <- intersect(pool, select)
  if (length(overlap) > 0L) {
    stop("Error: The following names appear in both pool and select: ",
         paste(overlap, collapse = ", "))
  }

  if (sum(sum(grepl("Z.hat",colnames(data))))) stop("Variable name beginning with Z.hat discovered, reserved for internal use. Please rename.")
  if (sum(sum(grepl("D.hat",colnames(data))))) stop("Variable name beginning with D.hat discovered, reserved for internal use. Please rename.")
  if (sum(sum(grepl("Q.hat",colnames(data))))) stop("Variable name beginning with Q.hat discovered, reserved for internal use. Please rename.")
  if ("Q" %in% colnames(data)) stop("Data contains variable named Q, which is reserved for internal use. Please rename.")

  tunetype=match.arg(tunetype,c("one","all"))
  if (tune.Qparameters!="none") {
    tune.Qparameters=match.arg(tune.Qparameters,several.ok=TRUE,c("all","sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty"))
  }
  if (tune.Cparameters!="none") {
    tune.Cparameters=match.arg(tune.Cparameters,several.ok=TRUE,c("all","sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty"))
  }
  if (tune.Zparameters!="none") {
    tune.Zparameters=match.arg(tune.Zparameters,several.ok=TRUE,c("all","sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty"))
  }
  if (tune.Yparameters!="none") {
    tune.Yparameters=match.arg(tune.Yparameters,several.ok=TRUE,c("all","sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty"))
  }
  if (tune.Dparameters!="none") {
    tune.Yparameters=match.arg(tune.Yparameters,several.ok=TRUE,c("all","sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty"))
  }
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

  if (sum(test %in% c("AHS","MW","BP"))==0&is.null(Y)==FALSE) Y=NULL
  time=rbind(start=time,checks=proc.time())

  ###################### 2 Prepare data #########################3
  XW=X

  if (sum(c("MW","AHS") %in% test)>0) {XWY=c(X,Y)} else {XWY=XW}

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
  make_group_folds(data,K = 2,cluster_name = cluster, fold_col = "sample",verbose = FALSE,diag_prefix=NULL)

  ##OPTIONAL INNER SPLIT
  if (is.null(inner.folds)==FALSE) {
    make_group_folds(data,K = inner.folds,cluster_name = cluster,fold_col = "cf_fold",verbose = FALSE,by_col="sample",diag_prefix=NULL)
    foldname="cf_fold"
  } else {
    foldname=NULL
  }

  if ((is.null(cluster)==TRUE)&(is.null(stack)==FALSE)) cluster="id_"

  time=rbind(time,prepare=proc.time())

  ############### 3 Discretize Z, D and Y into subsets ###############33

  if (is.null(Dsubsets)==FALSE) { ##bin treatment
    data <- binarize_var(
      data    = data,
      var     = D,
      ngroups = Dsubsets,
      gridtype = gridtypeD,
      wvar    = wvar
    )
  }

  if (is.null(Zsubsets)==FALSE) { ##instrument
    data <- binarize_var(
      data    = data,
      var     = Z,
      ngroups = Zsubsets,
      gridtype = gridtypeZ,
      wvar    = wvar
    )
  }


  if (sum(test %in% c("BP","K"))>0) { ##bin outcome(s)
    maxlevs=c();Ybin=c()
    for (Yno in Y) {
        data <- binarize_var(
        data    = data,
        var     = Yno,
        ngroups = Ysubsets,
        gridtype = gridtypeY,
        wvar    = wvar,
        newvar = paste0(Yno, ".bin")
      )
      col <- paste0(Yno, ".bin")
      maxlevs <- c(maxlevs, max(data[[col]], na.rm = TRUE))
      Ybin=c(Ybin,paste0(Yno,".bin"))
    }
  }

  n=nrow(data); J=length(unique(data[,get(..D)]))-1; K=length(unique(data[,get(..Z)]))-1;L=length(Y)

  results=c();tunable.Cparams=c()

  if (J>1&(sum(test %in% "simple")!=length(test))) stop("Multivalued treatments and testable conditions other than simple not supported")

  if (J==1&K==1&is.null(X)==TRUE&is.null(Y)==TRUE) {
    stop("Nothing to test with a binary treatment, a binary instrument and no other variables in data.")
  }

  XW=X
  XWY=c(X,Y)

  ##group common tree argments
  forest_opts=list(num.trees=max(50,num.trees/4),tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps)
  time=rbind(time,discretize=proc.time())

  ######################## 6a STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
  margins=c()

  ##HANDLE Z STACKS
  if (K>1) { ##Expand to all margins of Z
    minZ=min(data[,..Z]);maxZ=max(data[,..Z]);zvals=sort(unique(data[,get(..Z)]))
    times = 1L + 1L * as.integer(data[[Z]] > minZ & data[[Z]] < maxZ)
    data = data[rep.int(seq_len(nrow(data)), times)]

    data[, zmargin := {
      z <- get(Z)
      k <- match(z[1L], zvals)

      if (k == 1L || k == length(zvals)) {
        rep(z[1L], .N)
      } else {
        c(zvals[k], zvals[k + 1L])[seq_len(.N)]
      }
    }, by = "id_"]
    data[,(Z):=get(..Z)>=zmargin]

    margins=c(margins,"zmargin")
  } else { ##MAKE SURE Z is dummy!
    data[, (Z) := as.integer(get(..Z) == max(get(..Z), na.rm = TRUE))]
  }



  ##estimate Z.hat for each margin in stacked data (also if K==1!)
  if (is.null(X)==FALSE) {
    crossfit_hat(
      data,
      y_name = Z,
      x_names = X,
      margins = margins,
      weight_name = weight,
      forest_opts = c(forest_opts,Zparameters),
      mode=if ("Z" %in% crossfit) "across" else "within"
    )
  } else {
    data[,(paste0(Z,".hat")):=(mean(get(..Z))*.N-get(..Z))/(.N-1),by=c("sample",margins)] ##leave one out mean
  }
  if (normalize.Z==TRUE) { #Normalize propensity scores
    data[, (paste0(Z, ".hat")) :=
           get(paste0(..Z, ".hat")) * .N / sum(get(..Z) / get(paste0(..Z, ".hat"))),
         by = c("sample", margins)]

  }


  ##STACK MULTIPLE OUTCOMES
  if (L>1) {
    data=melt(data, measure = list(paste0(Y,rep(".bin",length(Y))),Y),value.name = c("Ystack.bin", "Ystack"),variable.name="outcome")
    data[,outcome:=Y,by=c("_id",margins)]
    data[,maxlevsY:=maxlevs,by=c("id_",margins)]
    margins=c(margins,"outcome")
    Y="Ystack"
  } else if (sum(test %in% c("BP","K"))>0) maxlevsY=maxlevs

  ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
  if ((sum(test %in% c("MW","AHS"))>0)&Y.res==TRUE)  {
    crossfit_hat(
      data,
      y_name = Y,
      x_names = X,
      folds=foldname,
      margins = margins,
      weight_name = weight,
      mode=if ("Y" %in% crossfit) "across" else "within",
      forest_opts = c(forest_opts,Yparameters)
    )
    data[,paste0(Y,".res"):=get(..Y)-get(paste0(..Y,".hat"))]
  }

  ##STACK MARGINS OF D
  if (J>1) { ##Stack data for all margins of D
    rows=rep(seq_len(nrow(data)),each=J)
    data= data[rows]    # replicate each row J times
    dvals=sort(unique(data[,get(..D)]))
    data[,dmargin:=..dvals[-1],by=c("id_",margins)]
    data[,D:=D>=dmargin,env=list(D=D)]
    margins=c(margins,"dmargin")
  }

  ##Estimate D.hat if test includes "K"
  if ("K" %in% test) {
    if (length(test)>1) i=which(data$condition=="K") else i=NULL
    crossfit_hat(
      data,
      i=i,
      y_name = D,
      x_names = X,
      folds=foldname,
      margins = margins,
      weight_name = weight,
      mode=if ("D" %in% crossfit) "across" else "within",
      forest_opts = c(forest_opts,Dparameters)
    )
  }

  ##Check for onesided noncompliance

  os_res <- test_one_sided_noncompliance(data = data, D = D, Z = Z, margins = margins)
  osmargins=margins

  ##Expand multiple conditions for testing (except K + BP, which has same def of Q - expand later)
  if (length(test)>1) {
    if (sum(test %in% c("BP","K"))==2) testexpand=c(test[!test %in% c("BP","K")],"BPK")
    else testexpand=test
    if (length(testexpand)>1) {
      data=data[rep(seq(.N),length(testexpand))]
      data[,condition:=testexpand,by=c("id_",margins)]
    } else {
      data[,condition:=ifelse(length(test)==2,"BPK",test)]
    }
    margins=c(margins,"condition")
    condition="condition"
  } else {condition=test}

  ## DROP TRIVIALLY SATISFIED ROWS DUE TO ONESIDED NONCOMPLIANCE,
  ## drop rows with condition="FS" or "AHS" in margins with onesided noncompliance (trivial)
  ## EXPAND BP,MW rows to equation=0:1 for nontrivially satisfied conditions
  ## (or keep just the nontrivially satisfied one)
  ## Expand to equation=0,1 only when BOTH equations are nontrivial

  if (length(osmargins) == 0L) {

    # global one-sided result (os_res should be 1 row)
    os_one <- isTRUE(os_res[["one_sided"]][1])
    os_eq  <- os_res[["trivial_equation"]][1]

    data[, `:=`(os_one_sided = os_res[["one_sided"]][1],
                os_triv_eq   = os_res[["trivial_equation"]][1])]

  } else {

    data[os_res, on = osmargins,
         `:=`(os_one_sided = i.one_sided,
              os_triv_eq   = i.trivial_equation)]

    # unmatched cells treated as not one-sided
    data[is.na(os_one_sided), os_one_sided := FALSE]
  }

  # ------------------------------------------------------------
  # 2) Drop trivially satisfied SIMPLE/AHS in one-sided cells
  # ------------------------------------------------------------

  data <- data[!(condition %chin% c("simple","AHS") & os_one_sided)]

  # ------------------------------------------------------------
  # EXPAND to EQ=0:1
  # ------------------------------------------------------------

  if (sum(test %in% c("BP","K","BPK","MW"))>0) {
    data[,expand_n:=1]
    data[(condition %chin% c("BP","K","BPK","MW")) & os_one_sided==FALSE,expand_n:=2]
    data=data[rep(seq(.N),expand_n)]
    data[(condition %chin% c("BP","K","BPK","MW")) & !os_one_sided,equation := seq_len(.N) - 1L,by = c("id_", margins)]
    data[(condition %chin% c("BP","K","BPK","MW")) & os_one_sided==TRUE,equation:=1-os_triv_eq]
    margins=c(margins,"equation")
    data[,expand_n:=NULL]
    }

  data[, c("os_one_sided", "os_triv_eq") := NULL]

  ##Expand to all groups of Ybin for BP, K conditions
  if (sum(test %in% c("BP","K"))>0) {
    data=data[rep(seq(.N), 1+maxlevsY*(condition %in% c("BPK","K","BP")))]
    data[condition %in% c("BP","K","BPK"),ybin:=(0:maxlevsY),by=c("id_",margins)]
    margins=c(margins,"ybin")
  }

  ##Create outcome variable Q in stacked data
  if (sum(test %in% c("simple","AHS"))>0) data[condition %in% c("simple","AHS"),Q:=D,env=list(D=D)]
  if (sum(test %in% c("BPK","K","BP"))>0) data[condition %in% c("BPK","BP","K"),Q:=equation*(D*(name==ybin))-(1-equation)*(1-D)*(name==ybin),env=list(D=D,name=paste0(Y,".bin"))]
  if ("MW" %in% test) data[condition=="MW",Q:=equation*((1-get(paste0(..Z,".hat")))*D*Z-get(paste0(..Z,".hat"))*D*(1-Z))+(1-equation)*(get(paste0(..Z,".hat"))*(1-D)*(1-Z)-(1-get(paste0(..Z,".hat")))*(1-D)*Z),env=list(Z=Z,D=D)] ##ERROR HERE!!! get(paste0(...))
  if ("K" %in% test) {
    data[condition %in% c("BPK","K"),D.hat:=D.hat*equation+(1-D.hat)*(1-equation)]
    data[condition %in% c("BPK","K"),D:=D*equation+(1-D)*(1-equation),env=list(D=D)]
    ##reverse treatment indicator for equation==0 and condition=="K"
  }

  ##TEST FOR SUPPORT and VARIATION IN Z AND Q within each group

  data[,nQ:=uniqueN(Q), by=c("sample",margins)]
  data[,nZ:=uniqueN(get(..Z)), by=c("sample",margins)]
  data[,n0:=sum(get(..Z)==0), by=c("sample",margins)]
  data[,n1:=sum(get(..Z)==1), by=c("sample",margins)]
  data[,sd_res:=sd(get(..Z)-get(paste0(..Z,".hat"))),by=c("sample",margins)]
  data[,bad:= (nQ<2| nZ<2 | n0<min_n | n1 <min_n | sd_res==0)]
  data=data[bad==FALSE]

  ##ALSO DROP OTHER SAMPLE PART WITHIN A MARGIN IF ONE PART HAS BEEN DROPPED
  data[,sdsample:=sd(sample),by=margins]
  data=data[sdsample>0]
  data[,sdsample:=NULL]

  if (nrow(data)==0) {
    stop("No observations with variation in Z and Q and overlap remains - identification issue.")
  }

  ##Estimate Q.hat in stacked data
  if (sum(test %in% c("simple","BP","K","BPK"))>0) {
    if (length(test)>1) i=which(data$condition %in%  c("simple","BP","K","BPK")) else i=NULL
    crossfit_hat(
      data,
      i=i,
      y_name = "Q",
      x_names = X,
      folds=foldname,
      margins = margins,
      weight_name = weight,
      mode=if ("Q" %in% crossfit) "across" else "within",
      forest_opts = c(forest_opts,Qparameters)
    )
  }

  if ("AHS" %in% test) {
    if (length(test)>1) i=which(data$condition=="AHS") else i=NULL
    crossfit_hat(
      data,
      i=i,
      y_name = "Q",
      folds=foldname,
      x_names = c(X,paste0(Y,".res")),
      margins = margins,
      weight_name = weight,
      mode=if ("Q" %in% crossfit) "across" else "within",
      forest_opts = Qparameters
    )
  }

  ##Split BP and K conditions if doing both
  if ((sum(test %in% c("BP","K"))>0)) {
    if (sum(test %in% c("BP","K"))==2) data=data[rep(seq(.N),1+1*(condition=="BPK"))]
    data[condition=="BPK",condition:=test[test %in% c("BP","K")],by=c("id_",margins)]
  }

  time=rbind(time,stack_nuisance=proc.time())

  ########## ESTIMATE ALL CAUSAL/REGRESSION/IV FORESTS AND  predict in/out of sample ##########
  forest_opts=list(num.trees=max(50,num.trees),tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps)
  if ("C" %in% crossfit) foldname=NULL #Do not crossfit causal forest, just the nuissances

  if (sum(test %in% c("simple","BP"))>0) {
    if (length(test)>1) i=which(data$condition %in% c("simple","BP")) else i=NULL
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
               forest_opts = c(forest_opts,Cparameters),
               aipw.clip=aipw.clip,
               shrink=(shrink>0))
  }

  if (sum(test == "AHS")>0) {
    if (length(test)>1) i=which(data$condition=="AHS") else i=NULL
    fit_models(data,
               forest_type = "causal",
               i=i,
               y_name="Q",
               x_names=c(X,paste0(Y,".res")),
               w_name=Z,
               folds=foldname,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = c(forest_opts,Cparameters),
               aipw.clip=aipw.clip,
               shrink=(shrink>0))
  }

  if (sum(test == "MW")>0) {
    if (length(test)>1) i=which(data$condition=="MW") else i=NULL
    fit_models(data,
               forest_type = "regression",
               i=i,
               y_name="Q",
               x_names=c(X,paste0(Y,".res")),
               folds=foldname,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = c(forest_opts,Cparameters),
               aipw.clip=aipw.clip,
               shrink=(shrink>0))
  }

  if (sum(test == "K")>0) {
    if (length(test)>1) i=which(data$condition=="K") else i=NULL
    fit_models(data,
               forest_type = "instrumental",
               i=i,
               y_name="Q",
               x_names=X,
               w_name=D,
               folds=foldname,
               z_name=Z,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = c(forest_opts,Cparameters),
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


  time=rbind(time,causal_forest=proc.time())

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

  if ("forest" == testtype) res=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,x_names=X,pool=poolmargins,select=selectmargins,gridpoints=gridpoints,margins=margins)
  if ("CART" == testtype) res=CART_test(data, x_names=X,margins=margins,weight=weight,cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,preselect=preselect,cluster=cluster,select=selectmargins,rpart_options=rpart_options)


  time=rbind(time,find_and_test=proc.time())


  ################ 7: Multiple hypothesis testing and output #####################
    if (nrow(res$results[train==FALSE])==1) {
      res$minwhere=NA
      res$minp=rep(res$results[train==FALSE,p.raw],6)
      names(res$minp)=paste0("p.",c("raw","holm","hochberg","BH","BY","CCT"))
    } else {
      for (m in c("holm","hochberg","BH","BY")) {
        res$results[train==FALSE&is.na(t)==FALSE,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }
      byv=c("sample",margins)[!c("sample",margins) %in% pool]
      res$minwhere=res$results[train == FALSE & is.finite(p.raw)][which.min(p.raw), ..byv]
      res$minp=apply(res$results[train==FALSE&is.na(t)==FALSE,c("p.raw","p.holm","p.hochberg","p.BH","p.BY")],2,min)
      res$minp=c(res$minp,p.CCT=cct_pvalue(res$results[train==FALSE,p.raw]))
    }

    res$global[,p.raw:=pnorm(t)]
    if (nrow(res$global)>1) {
      for (m in c("holm","hochberg","BH","BY")) {
        res$global[,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }
    }

  time=rbind(time,finalize=proc.time())
  time = time[-1, , drop = FALSE] - time[-nrow(time), , drop = FALSE]
  time=time[,1:3]
  time=rbind(time,total=colSums(time))
  return=c(res,list(time=time),obs=obs)
  return(return)
}
