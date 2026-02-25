#' @param data A data.table containing the analysis sample.
#' @param ... Other arguments (document these properly).
#' @return A list with components \code{results}, \code{meta}, and optional fitted objects.
#' @export

montest=function(data,D,Z,X=NULL,Y=NULL,W=NULL,test=NULL,inner.folds=5,crossfit.forest=TRUE,
                 normalize.Z=TRUE,aipw_clip=1e-3,weight=NULL,cluster=NULL,num.trees=2000,seed=10101,minsize=50,
                 gridtypeY="equisized",gridtypeD="equisized",gridtypeZ="equisized",sim=FALSE,
                 Ysubsets = 4, Dsubsets = 4,Zsubsets=4,Y.res=TRUE,testtype="forest",
                 gridpoints=NULL,min_n=1L,pool="all", ##forest opts
                 cp=0,maxrankcp=10L,alpha=0.05,prune=TRUE,preselect="negative", ##CART opts
                 Zparameters=list(),Yparameters=list(),Qparameters=list(),Dparameters=list(),Cparameters=list(),
                 tune.Qparameters="none",tune.Zparameters="none",tune.Cparameters="none",tune.Yparameters="none",tune.Dparameters="none",
                 tune.num.trees=200,tune.num.reps=50,tune.num.draws=1000,tunetype="one" ##tuning options
){

  time=proc.time()
  set.seed(seed)


  ################### 1 CHECK INPUT #####################
  if ((is.null(cluster)==FALSE)&("CART" %in% testtype)) stop("Clustering not supported with testtype = CART.")
  testtype=match.arg(testtype,c("forest","CART"),several.ok=TRUE)
  if (!is.null(aipw_clip)) {
    stopifnot(
      is.numeric(aipw_clip),
      length(aipw_clip) == 1L,
      is.finite(aipw_clip),
      aipw_clip > 0,
      aipw_clip < 1
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
  if (sum(pool=="all")==1) pool=c("zmargin","dmargin","ybin","condition","equation","outcome","sample")
  if (sum(pool=="none")==1) pool=c()
  pool=match.arg(pool,c("zmargin","dmargin","ybin","condition","equation","outcome","sample"),several.ok=TRUE)
  #if (testtype=="CART") pool=c()

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
  XW=c(X,W)

  if (sum(c("MW","AHS") %in% test)>0) {XWY=c(X,W,Y)} else {XWY=XW}

  data=data.table(data)
  allvars=c(X,W,Y,D,Z,weight,cluster)
  data=data[,..allvars]
  dropped=sum(!complete.cases(data))
  if (dropped>0) message(paste("Note: dropped ",dropped," observations with missing data on one or more input variables."))
  data=data[complete.cases(data)]
  n=nrow(data)
  if (is.null(cluster)==FALSE) {
    G=length(unique(data[,get(..cluster)]))
    if (G<=2*minsize) stop("Number of clusters is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")
  } else if (n<=2*minsize) stop("Number of observations is smaller than 2x minsize. There is not enough data to split the sample and test in a large enough sample. Reconsider specification or reduce minsize.")

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
    if (length(unique(data[,get(..D)]))>Dsubsets){
      if (gridtypeD=="equidistant") {
        data[,(D):=as.integer(cut(get(..D), breaks = seq(from = min(get(..D)) - 0.001, to = max(get(..D)) + 0.001, length.out = Dsubsets + 1)))-1]
      }
      else {
        data[,(D):={
          w <- if (is.na(..wvar)) NULL else get(..wvar)
          as.integer(cut(get(..D),breaks=c(-Inf,unique(weighted_quantile(get(..D),w=w,probs=(1:(Dsubsets-1))/Dsubsets)),Inf)))-1
        }]

      }
    }
  }

  if (is.null(Zsubsets)==FALSE) { ##bin instrument
    if (length(unique(data[,get(..Z)]))>Zsubsets){
      if (gridtypeZ=="equidistant") {
        data[,(Z):=as.integer(cut(get(..Z), breaks = seq(from = min(get(..Z)) - 0.001, to = max(get(..Z)) + 0.001, length.out = Zsubsets + 1)))-1]
      }
      else {
        data[,(Z):={
          w <- if (is.na(..wvar)) NULL else get(..wvar)
          as.integer(cut(get(..Z),breaks=c(-Inf,unique(weighted_quantile(get(..Z),w=w,probs=(1:(Zsubsets-1))/Zsubsets)),Inf)))-1
        }]
      }
    }
  }

  if (sum(test %in% c("BP","K"))>0) { ##bin outcome(s)
    maxlevs=c();Ybin=c()
    for (Yno in Y) {
      if (length(unique(data[,get(..Yno)]))>Ysubsets) {
        if (gridtypeY=="equidistant") {
          data[,(paste0(Yno,".bin")):=as.integer(cut(get(..Yno), breaks = seq(from = min(get(..Yno)) - 0.001, to = max(get(..Yno)) + 0.001, length.out = Ysubsets + 1)))]
        }
        else {
          data[,(paste0(Yno,".bin")):={
            w <- if (is.na(..wvar)) NULL else get(..wvar)
            as.integer(cut(get(..Yno),breaks=c(-Inf,unique(weighted_quantile(get(..Yno),w=w,probs=(1:(Ysubsets-1))/Ysubsets)),Inf)))
          }]
        }
      } else {
        data[,(paste0(Yno,".bin")):=as.numeric(cut(get(..Yno),breaks=length(unique(get(..Yno)))))]
      }
      maxlevs=c(maxlevs,data[, max(get(paste0(Yno, ".bin")))])
      Ybin=c(Ybin,paste0(Yno,".bin"))
    }
  }

  n=nrow(data); J=length(unique(data[,get(..D)]))-1; K=length(unique(data[,get(..Z)]))-1;L=length(Y)

  results=c();tunable.Cparams=c()

  if (J>1&(sum(test %in% "simple")!=length(test))) stop("Multivalued treatments and testable conditions other than simple not supported")

  if (J==1&K==1&is.null(X)==TRUE&is.null(W)==TRUE&is.null(Y)==TRUE) {
    stop("Nothing to test with a binary treatment, a binary instrument and no other variables in data.")
  }

  XW=c(X,W)
  XWY=c(X,W,Y)

  ##group common tree argments
  forest_opts=list(num.trees=max(50,num.trees/4),tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps)
  time=rbind(time,discretize=proc.time())

  ######################## 6a STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
  margins=c()

  ##HANDLE Z STACKS
  if (K>1) { ##Expand to all margins of Z
    times = 1L + 1L * as.integer(data[[Z]] > 0 & data[[Z]] < K)
    data = data[rep.int(seq_len(nrow(data)), times)]
    data[,zmargin:=seq(.N)+get(..Z)-1*(get(..Z)>0),by=id_]
    data[,(Z):=get(..Z)>=zmargin]
    margins=c(margins,"zmargin")
  }

  ##estimate Z.hat for each margin in stacked data (also if K==1!)
  if (is.null(X)==FALSE) {
    crossfit_hat(
      data,
      y_name = Z,
      x_names = X,
      margins = margins,
      weight_name = weight,
      forest_opts = c(forest_opts,Zparameters)
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
      forest_opts = c(forest_opts,Yparameters)
    )
    data[,paste0(Y,".res"):=get(..Y)-get(paste0(..Y,".hat"))]
  }

  ##STACK MARGINS OF D
  if (J>1) { ##Stack data for all margins of D
    rows=rep(seq_len(nrow(data)),each=J)
    data= data[rows]    # replicate each row J times
    data[,dmargin:=seq_len(J),by=c("id_",margins)]
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
    data=data[rep(seq(.N), 1+(maxlevsY-1)*(condition %in% c("BPK","K","BP")))]
    data[condition %in% c("BP","K","BPK"),ybin:=(1:maxlevsY),by=c("id_",margins)]
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
  if (crossfit.forest==FALSE) foldname=NULL #Do not crossfit causal forest, just the nuissances

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
               aipw_clip=aipw_clip)
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
               aipw_clip=aipw_clip)
  }

  if (sum(test == "MW")>0) {
    if (length(test)>1) i=which(data$condition=="MW") else i=NULL
    fit_models(data,
               forest_type = "regression",
               i=i,
               y_name="Q",
               x_names=X,
               folds=foldname,
               margins = margins,
               weight_name = weight,
               cluster_name = cluster,
               forest_opts = c(forest_opts,Cparameters),
               aipw_clip=aipw_clip)
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
               aipw_clip=aipw_clip)
  }


  time=rbind(time,causal_forest=proc.time())

  ######################################## FIND OPTIMAL SUBSET TO TEST AND TEST IN OPPOSITE SAMPLE #####################
  poolmargins=pool[pool %in% c(margins,"sample")]
  if (is.null(poolmargins)==TRUE) poolmargins=character(0)

  res=list()
  if (sim==TRUE) {
    poollist=list(margins[!margins %in% "condition"],c(margins[!margins %in% "condition"],"sample"),c(margins,"sample"),margins)
    treelist=c("CART","forest")

    for (p in 1:4) {
      res[[paste0(c("forest",p),collapse="_")]]=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,x_names=X,pool=poollist[[p]],gridpoints=gridpoints,margins=margins)
      res[[paste0(c("CART",p),collapse="_")]]=CART_test(data, x_names=X,margins=margins,weight=weight,cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,preselect=preselect,cluster=cluster,pool=poollist[[p]])
    }
  } else {
  if ("forest" %in% testtype) res[["forest"]]=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,x_names=X,pool=poolmargins,gridpoints=gridpoints,margins=margins)
  if ("CART" %in% testtype) res[["CART"]]=CART_test(data, x_names=X,margins=margins,weight=weight,cp = cp,maxrankcp = maxrankcp,alpha = alpha,prune = prune,  minsize = minsize,preselect=preselect,cluster=cluster,pool=poolmargins)
  }
  time=rbind(time,find_and_test=proc.time())


  ################ 7: Multiple hypothesis testing and output #####################
  for (r in 1:length(res)) {
    if (nrow(res[[r]]$results[train==FALSE])==1) {
      minwhere=NA
      minp=res[[r]]$results[train==FALSE,p.raw]
    } else {
      for (m in c("holm","hochberg","BH","BY")) {
        res[[r]]$results[train==FALSE&is.na(t)==FALSE,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }
      byv=c("sample",margins)[!c("sample",margins) %in% pool]
      minwhere=res[[r]]$results[train == FALSE & is.finite(p.raw)][which.min(p.raw), ..byv]
      res[[r]]$minp=apply(res[[r]]$results[train==FALSE&is.na(t)==FALSE,c("p.raw","p.holm","p.hochberg","p.BH","p.BY")],2,min)
      res[[r]]$minp=c(res[[r]]$minp,p.CCT=cct_pvalue(res[[r]]$results[train==FALSE,p.raw]))
    }

    res[[r]]$global[,p.raw:=pnorm(t)]
    if (nrow(res[[r]]$global)>1) {
      for (m in c("holm","hochberg","BH","BY")) {
        res[[r]]$global[,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }
    }

  }

  time=rbind(time,finalize=proc.time())
  time = time[-1, , drop = FALSE] - time[-nrow(time), , drop = FALSE]
  time=time[,1:3]
  time=rbind(time,total=colSums(time))
  return=c(res,list(time=time))
  return(return)
}
