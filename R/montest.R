  #' @param infile Path to the input file
  #' @return A matrix of the infile
  #' @export

  montest=function(data,D,Z,X=NULL,Y=NULL,W=NULL,treefraction=0.5,test=NULL,
                   normalize=TRUE,stack=NULL,treetype="forest",
                   gridtypeY="equidistant",gridtypeD="equidistant",gridtypeZ="equidistant",
                   Ysubsets = 4, Dsubsets = 4,Zsubsets=4,Y.res=TRUE,saveforest=F,
                   weight=NULL,cluster=NULL,num.trees=2000,seed=10101,minsize=50,
                   maxrankcp=5,prune=TRUE,cp=0,alpha=0.05,preselect="nonpositive", ##CART options
                   Zparameters=NULL,Yparameters=NULL,Qparameters=NULL,Dparameters=NULL,Cparameters=NULL,
                   tune.Qparameters="none",tune.Zparameters="none",tune.Cparameters="none",tune.Yparameters="none",tune.Dparameters="none",
                   tune.num.trees=200,tune.num.reps=50,tune.num.draws=1000,tunetype="one" ##tuning options
                   ){

    time=proc.time()
    set.seed(seed)


    ################### 1 CHECK INPUT #####################
    if (is.null(test)==TRUE) {
      if (is.null(Y)==TRUE) {
        test="simple"
      } else test="all"
    }

    test=match.arg(test,c("simple","BP","MW","AHS","all"),several.ok=TRUE) ##took out support for Kitagawa
    if ("all" %in% test) {
      test=c("simple","BP","MW","AHS")
    }

    if (is.null(Y)==TRUE) {
      if (sum(!test %in% "simple")>0) {
      stop("Other tests than simple may not be used when Y is not specified. Specify the Y argument or use test=simple")
      }
    }

    if ((Ysubsets<=1)|(Dsubsets<=1)|(Zsubsets<=1)) stop("Ysubsets, Dsubsets and Zsubsets must be integers larger than 1")

    if (is.null(stack)==TRUE&treetype=="forest") stack=TRUE
    if (is.null(stack)==TRUE&treetype=="CART") stack=FALSE
    if (is.logical(stack)==FALSE) stop("stack must be a logical (TRUE/FALSE)")

    if (treetype=="CART"&stack==TRUE) stop("Cannot stack margins with the CART algorithm because it doesn't support clustering.")
    if (treetype=="CART"&is.null(cluster)==FALSE) stop("Cannot combine the CART algorithm with clustering.")

    if (sum(sum(grepl("Z.hat",colnames(data))))) stop("Variable name beginning with Z.hat discovered, reserved for internal use. Please rename.")
    if (sum(sum(grepl("D.hat",colnames(data))))) stop("Variable name beginning with D.hat discovered, reserved for internal use. Please rename.")
    if (sum(sum(grepl("Q.hat",colnames(data))))) stop("Variable name beginning with Q.hat discovered, reserved for internal use. Please rename.")
    if ("Q" %in% colnames(data)) stop("Data contains variable named Q, which is reserved for internal use. Please rename.")

    treetype=match.arg(treetype,c("forest","CART"))
    preselect=match.arg(preselect,c("minimum","none","negative","nonpositive"))
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
    compute.oob.predictions=fifelse(treetype=="forest",TRUE,FALSE)
    gridtypeZ=match.arg(gridtypeZ,c("equidistant","equisized"))
    gridtypeY=match.arg(gridtypeY,c("equidistant","equisized"))
    gridtypeD=match.arg(gridtypeD,c("equidistant","equisized"))

    if ((gridtypeZ=="equisized"|gridtypeY=="equisized"|gridtypeD=="equisized")&is.null(weight)==FALSE) stop("Cannot combine weights and quantile-based splits of Y, D or Z.")

    if ((length(D)!=1)|(!(D %in% colnames(data)))) {
      stop("Argument D must be the name of a single column in data")
    }

    if (is.null(weight)==FALSE) {
      if ((length(weight)!=1)|(!(weight %in% colnames(data)))) {
        stop("Argument weight must be the name of a single column in data")
      }
    }

    if (is.null(cluster)==FALSE) {
      if ((length(cluster)!=1)|(!(cluster %in% colnames(data)))) {
        stop("Argument cluster must be the name of a single column in data")
      }
    }
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

    if (sum(test %in% c("AHS","MW"))==0&is.null(Y)==FALSE) Y=NULL

    ###################### 2 Prepare data #########################3
    XW=c(X,W)
    if (sum(c("MW","AHS") %in% test)>0) {XWY=c(X,W,Y)} else {XWY=XW}

    data=data.table(data)
    allvars=c(X,W,Y,D,Z,weight,cluster)
    data=data[,..allvars]
    data=data[complete.cases(data)]
    n=nrow(data)

    data[,id:=(1:n)]
    if (is.null(cluster)==TRUE) {
      data[,sample:=1+1*(id %in% sample(n,size=floor(treefraction*n), replace = FALSE))]
    } else {
      s2=sample(as.matrix(unique(data[,cluster,with=F])),size=floor(treefraction*nrow(unique(data[,cluster,with=F]))),replace=TRUE)
      data[,sample:=1+1*(cluster %in% s2)]
      rm(s2)
    }
    if ((is.null(cluster)==TRUE)&(is.null(stack)==FALSE)) cluster="id"

    ############### 3 Discretize Z, D and Y into subsets ###############33
    ##TODO: ALLOW FOR WEIGHTED QUANTILES

    if (is.null(Dsubsets)==FALSE) { ##bin treatment
      unique(data[,D,env=list(D=D)])
      if (length(unique(data[,D,env=list(D=D)]))>Dsubsets){
        if (gridtypeD=="equidistant") {
          data[,D:=as.numeric(cut(D, breaks = seq(from = min(D) - 0.001, to = max(D) + 0.001, length.out = Dsubsets + 1)))-1,env=list(D=D)]
        }
        else {
          data[,D:=as.numeric(cut(D, breaks = c(-Inf, quantile(D, seq(1/Dsubsets, 1 - 1/Dsubsets, 1/Dsubsets)), Inf)))-1,env=list(D=D)]
        }
      }
    }

    if (is.null(Zsubsets)==FALSE) { ##bin instrument
      if (length(unique(data[,Z,env=list(Z=Z)]))>Zsubsets){
        if (gridtypeZ=="equidistant") {
          data[,Z:=as.numeric(cut(Z, breaks = seq(from = min(Z) - 0.001, to = max(Z) + 0.001, length.out = Zsubsets + 1)))-1,env=list(Z=Z)]
        }
        else {
          data[,Z:=as.numeric(cut(Z, breaks = c(-Inf, quantile(Z, seq(1/Zsubsets, 1 - 1/Zsubsets, 1/Zsubsets)), Inf)))-1,env=list(Z=Z)]
        }
      }
    }

    if (sum(test %in% c("BP","K"))>0) { ##bin outcome(s)
      maxlevs=c();Ybin=c()
      for (Yno in Y) {
      if (length(unique(data[,Yno,env=list(Yno=Yno)]))>Ysubsets) {
        if (gridtypeY=="equidistant") {
          data[,name:=as.numeric(cut(Yno, breaks = seq(from = min(Yno) - 0.001, to = max(Yno) + 0.001, length.out = Ysubsets + 1))),env=list(Yno=Yno,name=paste0(Yno,".bin"))]
        }
        else {
          data[,name:=as.numeric(cut(Yno, breaks = c(-Inf, quantile(Yno, seq(1/Ysubsets, 1 - 1/Ysubsets, 1/Ysubsets)), Inf))),env=list(Yno=Yno,name=paste0(Yno,".bin"))]
        }
      } else {
        data[,name:=as.numeric(cut(Yno,breaks=length(unique(Yno))))-1,env=list(Yno=Yno,name=paste0(Yno,".bin"))]
        }
      maxlevs=c(maxlevs,max(data[,paste0(Yno,".bin"),with=F]))
      Ybin=c(Ybin,paste0(Yno,".bin"))
      }
    }

    n=nrow(data); J=dim(unique(data[,D,with=F]))[1]-1; K=dim(unique(data[,Z,with=F]))[1]-1;L=length(Y)

    results=c();tunable.Cparams=c()

    if (J>1&(sum(test %in% "simple")!=length(test))) stop("Multivalued treatments and testable conditions other than simple not supported")

    if (J==1&K==1&is.null(X)==TRUE&is.null(W)==TRUE&is.null(Y)==TRUE) {
      stop("Nothing to test with a binary treatment, a binary instrument and no other variables in data.")
    }

    XW=c(X,W)
    XWY=c(X,W,Y)

    ############ GROUP COMMON FOREST OPTIONS #############

    ##group common tree argments
    forest_opts=list(num.trees=max(50,num.trees/4),tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps)

    ######################## 6 STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
    if (stack==TRUE) {
      margins=c()

      ##HANDLE Z STACKS
        if (K>1) { ##Expand to all margins of Z
          data=data[rep(seq(.N),1+1*(Z>0&Z<K)),env=list(Z=Z)]
          data[,zmargin:=seq(.N)+Z-1*(Z>0),by=id,env=list(Z=Z)]
          data[,Z:=Z>=zmargin,env=list(Z=Z)]
          margins=c(margins,"zmargin")
        }

        ##estimate Z.hat for each margin in stacked data (also if K==1!)
        if (is.null(X)==FALSE) {
        data[,paste0(Z,".hat"):=do.call(regression_forest, materialize_args(.SD,model="rf",
            y_name=..Z, x_names=..XW,forest_opts=c(forest_opts,Zparameters),
            weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
            cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
            $prediction,by=c("sample",margins),env=list(Z=Z)]
        if (normalize==TRUE) { #Normalize propensity scores
          data[,Z.hat:=.N*Z.hat/sum(Z/Z.hat),by=c("sample",margins),env=list(Z=Z)]
        }
      } else {
      data[,Z.hat:=(mean(Z)*.N-Z)/(.N-1),by=c("sample",margins),env=list(Z=Z)] ##leave one out mean
      }


      ##STACK MULTIPLE OUTCOMES
      if (L>1) {
        data=melt(data, measure = list(paste0(Y,rep(".bin",length(Y))),Y),value.name = c("Ystack.bin", "Ystack"),variable.name="outcome")
        data[,outcome:=Y,by=c("id",margins)]
        data[,maxlevsY:=maxlevs,by=c("id",margins)]
        margins=c(margins,"outcome")
        Y="Ystack"
      } else maxlevsY=maxlevs

      ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
      if ((sum(test %in% c("MW","AHS"))>0)&Y.res==TRUE)  {
        data[,paste0(Y,".res"):=Y-do.call(regression_forest, materialize_args(.SD,model="rf",
             y_name=..Y, x_names=..XW,forest_opts=c(forest_opts,Zparameters),
              weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
             cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
             $prediction,by=c("sample",margins),env=list(Y=Y)]
        }

      ##STACK MARGINS OF D
      if (J>1) { ##Stack data for all margins of D
            data=cbind(CJ(dmargin=(1:J),id=data$id),data[,!"id"])
            data[,D:=D>=dmargin,env=list(D=D)]
            margins=c(margins,"dmargin")
          }

      ##Estimate D.hat if test includes "K"
      if ("K" %in% test) {
        data[,paste0(D,".hat"):=do.call(regression_forest, materialize_args(.SD,model="rf",
            y_name=..D, x_names=..XW,forest_opts=c(forest_opts,Zparameters),
             weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
            cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
             $prediction,by=c("sample",margins),env=list(D=D)]
        }

        ##Expand multiple conditions for testing (except K + BP, which has same def of Q - expand later)
          if (length(test)>1) {
              if (sum(test %in% c("BP","K"))==2) testexpand=c(test[!test %in% c("BP","K")],"BPK")
              else testexpand=test
              if (length(testexpand)>1) {
                data=data[rep(seq(.N),length(testexpand))]
                data[,condition:=testexpand,by=c("id",margins)]
              } else {
                data[,condition:=ifelse(length(test)==2,"BPK",test)]
              }
            margins=c(margins,"condition")
          } else {condition=test}

        ##Expand to a=0,1 for BP, K and MW conditions
          if (sum(test %in% c("MW","BP","K"))>0) {
            data=data[rep(seq(.N),1+condition %in% c("BPK","MW","BP","K"))]
            data[condition %in% c("BPK","MW","BP","K"),above:=seq(.N)-1*(condition %in% c("BPK","MW","BP","K")),by=c("id",margins)]
            margins=c(margins,"above")
            }

        ##Expand to all groups of Ybin for BP, K conditions
           if (sum(test %in% c("BP","K"))>0) {
              data=data[rep(seq(.N), 1+(maxlevsY-1)*(condition %in% c("BPK","K","BP")))]
              data[condition %in% c("BP","K","BPK"),ybin:=(1:maxlevsY),by=c("id",margins)]
              margins=c(margins,"ybin")
           }

        ##Create outcome variable Q in stacked data
        if (sum(test %in% c("simple","AHS"))>0) data[condition %in% c("simple","AHS"),Q:=D,env=list(D=D)]
        if (sum(test %in% c("BPK","K","BP"))>0) data[condition %in% c("BPK","BP","K"),Q:=above*(D*(name==ybin))-(1-above)*(1-D)*(name==ybin),env=list(D=D,name=paste0(Y,".bin"))]
        if ("MW" %in% test) data[condition=="MW",Q:=above*((1-Z.hat)*D*Z-Z.hat*D*(1-Z))+(1-above)*(Z.hat*(1-D)*(1-Z)-(1-Z.hat)*(1-D)*Z),env=list(Z=Z,D=D)]
        if ("K" %in% test) {
            data[condition %in% c("BPK","K"),D.hat:=D.hat*above+(1-D.hat)*(1-above)]
            data[condition %in% c("BPK","K"),D:=D*above+(1-D)*(1-above),env=list(D=D)]
            ##reverse treatment indicator for above==0 and condition=="K"
        }

        ##drop groups/margins/conditions/outcomes with constant Q
        data[,rQ:=max(Q)-min(Q),by=c("sample",margins)]
        data=data[,rQ:=min(rQ),by=margins]
        data[,rQ:=NULL]

        ##Estimate Q.hat in stacked data
        if (sum(test %in% c("simple","BP","K","BPK"))>0) {
          data[condition %in% c("simple","BPK","BP","K"),Q.hat:=do.call(regression_forest, materialize_args(.SD,model="rf",
               y_name="Q", x_names=..XW,forest_opts=c(forest_opts,Zparameters),
               weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
               cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))$prediction,by=c("sample",margins)]
        }
        if ("AHS" %in% test) {
          data[condition=="AHS",Q.hat:=do.call(regression_forest, materialize_args(.SD,model="rf",
             y_name="Q", x_names=..XWY,forest_opts=c(forest_opts,Zparameters),
             weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
             cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))$prediction,by=c("sample",margins)]
          }

        ##Split BP and K conditions if doing both
        if ((sum(test %in% c("BP","K"))>0)) {
            if (sum(test %in% c("BP","K"))==2) data=data[rep(seq(.N),1+1*(condition=="BPK"))]
            data[condition=="BPK",condition:=test[test %in% c("BP","K")],by=c("id",margins)]
        }

    ########## ESTIMATE ALL CAUSAL/REGRESSION/IV FORESTS AND  predict in/out of sample ##########

        driver_name <- "condition"

        model_spec <- list(
          "simple" = list(model = "cf", xvars = c(X, W)),          # CF: Q ~ Z ; cov: X,W
          "MW"     = list(model = "rf", xvars = c(X, W, Y)),        # RF: Q ~ X,W,Y
          "BP"     = list(model = "cf", xvars = c(X, W)),           # CF: Q ~ Z ; cov: X,W
          "AHS"    = list(model = "cf", xvars = c(X, W, Y)),        # CF: Q ~ Z ; cov: X,W,Y
          "K"      = list(model = "iv", xvars = c(X, W))            # IV: Q ~ D (Z) ; cov: X,W
        )

        # Housekeeping cols
        data[, idx__ := .I]
        data[, grp_id := .GRP, by = margins]

        # Helper: fetch model spec per key and validate X columns
        get_spec <- function(key) {
          sp <- model_spec[[as.character(key)]]
          if (is.null(sp)) stop(sprintf("No model_spec entry for '%s'", key))
          miss <- setdiff(sp$xvars, names(data))
          if (length(miss)) stop(sprintf("Missing X columns for '%s': %s", key, paste(miss, collapse = ", ")))
          sp
        }

        # Optional: column names for weights/clusters if present in parent frame
        wcol <- get0("weight",  inherits = TRUE, ifnotfound = NULL)   # e.g. "w"
        ccol <- get0("cluster", inherits = TRUE, ifnotfound = NULL)   # e.g. "id"

        # ---------- Single pass (fit CF/RF/IVF + write back) ----------
        data[, {
          # determine condition for this group
          key <- if (is.character(condition) && length(condition) == 1L) {
            condition
          } else {
            get(driver_name)[1L]
          }
          sp  <- get_spec(key)

          idx <- idx__                # row indices in the original data
          gid <- grp_id[1L]
          s0  <- sample[1L]           # this group's sample (e.g., 1 or 2)

          # build args with your helper (adds hats for cf/iv if present)
          if (sp$model == "rf") {
            args <- materialize_args(
              .sd        = .SD,
              model      = "rf",
              y_name     = "Q",
              x_names    = sp$xvars,
              forest_opts= c(get0("forest_opts", inherits = TRUE, ifnotfound = list()),
                             get0("Zparameters", inherits = TRUE, ifnotfound = list())),
              weight_col = wcol,
              cluster_col= ccol
            )
            fit <- do.call(grf::regression_forest, args)
            pred_in   <- as.numeric(predict(fit)$predictions)
            scores_in <- .SD[,Q]   # MW: scores = Q

          } else if (sp$model == "cf") {
            args <- materialize_args(
              .sd        = .SD,
              model      = "cf",
              y_name     = "Q",
              w_name     = ..Z,                   # treatment = Z
              x_names    = sp$xvars,              # covariates
              forest_opts= c(get0("forest_opts", inherits = TRUE, ifnotfound = list()),
                             get0("Zparameters", inherits = TRUE, ifnotfound = list())),
              weight_col = wcol,
              cluster_col= ccol
            )
            fit <- do.call(grf::causal_forest, args)
            pred_in   <- as.numeric(predict(fit)$predictions)
            scores_in <- as.numeric(get_scores(fit))

          } else if (sp$model == "iv") {
            args <- materialize_args(
              .sd        = .SD,
              model      = "iv",
              y_name     = "Q",
              w_name     = ..D,                   # treatment = D
              z_name     = ..Z,                   # instrument = Z
              x_names    = sp$xvars,              # covariates
              forest_opts= c(get0("forest_opts", inherits = TRUE, ifnotfound = list()),
                             get0("Zparameters", inherits = TRUE, ifnotfound = list())),
              weight_col = wcol,
              cluster_col= ccol
            )
            fit <- do.call(grf::instrumental_forest, args)
            pred_in   <- as.numeric(predict(fit)$predictions)
            scores_in <- as.numeric(get_scores(fit))

          } else {
            stop(sprintf("Unknown model '%s'", sp$model))
          }

          # write back in-sample
          set(data, i = idx, j = "pred",   value = pred_in)
          set(data, i = idx, j = "scores", value = scores_in)

          # cross-sample predict: use the model fit on this group to predict the "other" sample within the same margins
          other_idx <- data[grp_id == gid & sample != s0, which = TRUE]
          if (length(other_idx)) {
            X_other <- as.matrix(data[other_idx, sp$xvars, with = FALSE])
            pred_o  <- as.numeric(predict(fit, X_other)$predictions)
            set(data, i = other_idx, j = "pred_o", value = pred_o)
          }

          NULL
        }, by = .(grp_id, sample)]

        # cleanup
        data[, c("idx__","grp_id") := NULL]

        ######################################## FIND OPTIMAL SUBSET TO TEST AND TEST IN OPPOSITE SAMPLE #####################


        setorderv(data,cols=c("sample","pred"))

        if (is.null(weight)) weight=1
        data[,N:=seq_len(.N),by=sample]
        data[, `:=`(a = weight * Y, b = weight), env=list(weight=weight)] ## Per-row basics

        ## Within-cluster running totals (in current order)
        data[, `:=`(WgY = cumsum(a),Wg  = cumsum(b)), by = c("sample",cluster)]

        ## Global running totals
        data[, `:=`(SW  = cumsum(b),SWY = cumsum(a)), by=sample]
        data[,m:= SWY / SW,by=sample]

        ## Per-row deltas for the three cross-cluster aggregates
        ## Using (x_new^2 - x_old^2) and (x_new*y_new - x_old*y_old) to avoid shifts
        data[, `:=`(dTA2 =  WgY^2 - (WgY - a)^2,dTB2 =  Wg^2  - (Wg  - b)^2,dTAB =  WgY*Wg - (WgY - a)*(Wg - b))]

        ## Cumulative (global) aggregates across rows
        data[, `:=`(TA2 = cumsum(dTA2),TB2 = cumsum(dTB2),TAB = cumsum(dTAB)),by=sample]

        ## Number of unique clusters seen so far
        data[, G := cumsum(!duplicated(cluster)),by=sample,env=list(cluster=cluster)]

        ## Cluster-robust SE for the weighted mean (CRV1 with small-sample adj.)
        ## sumS2 = (TA2 - 2*m*TAB + m^2*TB2)/SW^2
        data[, sumS2 := (TA2 - 2*m*TAB + (m^2)*TB2) / (SW^2)]
        data[, se := sqrt( (G / pmax.int(G - 1L, 1L)) * sumS2 )]
        data[G < 2, se := NA_real_]

        ## t-stat for H0: mean = 0
        data[, t := m / se]

        data[, c("a","b","WgY","Wg","dTA2","dTB2","dTAB","TA2","TB2","TAB","sumS2") := NULL]

        ##FIND DUAL-constraint TAU CUTOFF that ensures minsize clusters in both samples
        tau <- data[G==minsize,min(pred),by=sample]
        setorderv(data,cols=c("sample","pred_o"))

        tau <- cbind(tau,data[cumsum(!duplicated(cluster))==minsize,min(pred_o),by=sample,env=list(cluster=cluster)])

        tau=c(max(tau[1,2],tau[2,4]),max(tau[2,2],tau[1,4]))
        data[,tau:=fifelse(sample==1,tau[1],tau[2])]

        ###TRAIN SAMPLE RESULTS
        res=data[pred>tau, .SD[which.min(t)], by = sample,.SDcols=c("G","N","m","se","t","pred")] ##FIND OPTIMAL CUTOFF


        #####PERFORM TEST IN TRAIN SAMPE- CR1 cluster robust inference ############
        clust=as.formula(paste0("~",cluster))
        fe=feols(scores~1,data=data,vcov=clust,split=~sample)
        GN=cbind(data[pred_o<=tau,uniqueN(cluster),by=sample,env=list(cluster=cluster)][,V1],data[pred_o<=tau,.N,by=sample,env=list(cluster=cluster)][,N])

        res=cbind(res,GN,rbind(fe$`sample.var: sample; sample: 1`$coeftable[1,1:3],fe$`sample.var: sample; sample: 2`$coeftable[1,1:3]))

        colnames(res)=c("sample","G.train","N.train","coef.train","stderr.train","t.train","tau cutoff","G.est","N.est","coef.est","stderr.est","t.est")

        Xmeans <- data[pred_o <= tau,lapply(.SD, mean, na.rm = TRUE),by = sample,.SDcols = get0("XW", inherits = TRUE)]
        Xmeans_all <- data[,lapply(.SD, mean, na.rm = TRUE),by = sample,.SDcols = get0("XW", inherits = TRUE)]
        XSD <- data[,lapply(.SD, sd, na.rm = TRUE),by = sample,.SDcols = get0("XW", inherits = TRUE)]

        marginsmat_all=data[,.N,by=c("sample",margins)]
        marginsmat=data[pred_o<=tau,.N,by=c("sample",margins)]
        marginsmat=merge(marginsmat,marginsmat_all, by=c("sample",margins),all=TRUE)
        marginsmat[,share_all:=N.y/sum(N.y), by=sample]
        marginsmat[,share:=ifelse(is.na(N.x),0,N.x)/sum(N.x,na.rm=TRUE), by=sample]
        marginsmat[,c("N.x","N.y"):=NULL]
    }

    ################ END: Multiple hypothesis testing and output #####################
      res[is.na(t.est)==FALSE,p.raw:=pnorm(t.est)]
      for (m in c("holm","hochberg","BH","BY")) {
        res[is.na(t.est)==FALSE,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }

      minsample=res[which.min(p.raw),sample]
      minp=apply(res[is.na(t.est)==FALSE,c("p.raw","p.holm","p.hochberg","p.BH","p.BY")],2,min)
      time=proc.time()-time
      return=list(results=res,minsample=minsample,minp=minp,time=time,margins=marginsmat,Xmeans_all=Xmeans_all,XSD=XSD,Xmeans=Xmeans)
      return(return)
  }


##################################################### HELPER FUNCTIONS ######################################

# THIS FUNCTION ESTIMATES A TREE USING CART, TEST IN EACH NODE AND DETERMINE relevance
 esttree=function(data,testsample,cp,maxrankcp,alpha,prune,minsize,preselect){
    tree=rpart(scores~.,data=as.data.frame(data[sample==testsample,!c("id","sample")]),method="anova",cp=cp,minbucket=minsize,weights=weight) # run tree with transformed outcome

    maxrankcp=min(maxrankcp, length(tree$cp[, 4]))
    maxcp=tree$cp[maxrankcp, 1]

    if (prune==TRUE){ #prune the tree based on cross-validation
      opcpid=which.min(tree$cp[, 4])
      opcp=tree$cp[opcpid, 1]
      tree=prune(tree, cp = max(maxcp,opcp))
    }      else { # prune the tree based on maximum complexity parameter only (no cross-validation)
      tree=prune(tree, cp = maxcp)
    }
    data[,leaf:=rpart.predict.leaves(tree,newdata=data,type = "where")] # create dummies for leaves (based on predictions))

    ##Drop leaves with no variation in scores
    if (sum(data[sample==1,sd(scores)==0,by=leaf][,V1])>0) data[sample==1&(leaf %in% data[sample==1,sd(scores)==0,by=leaf][V1==TRUE,leaf]),leaf:=NA]
    if (sum(data[sample==2,sd(scores)==0,by=leaf][,V1])>0) data[sample==2&(leaf %in% data[sample==2,sd(scores)==0,by=leaf][V1==TRUE,leaf]),leaf:=NA]

    data[,sample:=ifelse(sample==testsample,"train","est")]
    if (is.null(cluster)==FALSE) {
      data[,G:=cumsum(!duplicated(get(cluster))),by=.(leaf,sample)]
    }else {
      data[,G:=.N,by=.(leaf,sample)]
    }

    res=data[is.na(leaf)==FALSE,data.table(cbind(G=max(G),N=.N,feols(scores~1,data=.SD,vcov=clust,weights=weight)$coeftable)),by=.(leaf,sample)]
    colnames(res)=c("leaf","sample","G","N","coef","stderr","t","p")
    res=res[,p:=NULL]
    res=dcast(res,leaf~sample,value.var=c("G","N","coef","stderr","t"))
    if (preselect!="none"&nrow(res)>1) {
      if (preselect=="nonpositive") {res[,relevant:=ifelse(t_train>=qnorm(1-alpha/nrow(res))&(t_train>min(t_train)),0,1)]}
      else if (preselect=="negative") {res[,relevant:=ifelse(t_train>=qnorm(alpha/nrow(res))&(t_train>min(t_train)),0,1)]}
      else if (preselect=="minimum") {res[,relevant:=ifelse(t_train>min(t_train),0,1)]}
   } else {res[,relevant:=1]}
    res[relevant==0,c("coef_est","stderr_est","t_est"):=NA]
    setcolorder(res, c("leaf","G_train","N_train","coef_train","stderr_train","t_train","relevant","G_est","N_est","coef_est","stderr_est","t_est"))
    return(list(tree=tree,res=res))
 }

 weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75),
                  na.rm = TRUE,
                  interpolation = c("linear", "left", "right", "midpoint")) {
   interpolation <- match.arg(interpolation)

   # basic checks
   if (length(x) != length(w)) stop("x and w must have the same length.")
   if (na.rm) {
     keep <- is.finite(x) & is.finite(w)
     x <- x[keep]; w <- w[keep]
   }
   if (any(w < 0)) stop("weights must be nonnegative.")
   if (!length(x)) return(rep(NA_real_, length(probs)))
   s <- sum(w)
   if (s == 0) return(rep(NA_real_, length(probs)))

   # sort by x
   o <- order(x, na.last = NA)
   x <- x[o]; w <- w[o]

   # cumulative weight share in [0,1]
   cw <- cumsum(w) / s

   # helper to get quantile for a single p
   get_q <- function(p) {
     if (p <= 0) return(x[1L])
     if (p >= 1) return(x[length(x)])

     # first index where cw >= p
     k <- which.max(cw >= p)  # same as min(which(cw >= p)) but faster
     if (cw[k] == p) {
       # exact hit: handle by interpolation rule
       if (interpolation == "left")      return(x[k])
       if (interpolation == "right")     return(x[k])
       if (interpolation == "midpoint") {
         # average of the "mass at p": include any ties on x with same cw threshold
         j <- k
         while (j > 1L && cw[j - 1L] == cw[k]) j <- j - 1L
         i <- k
         while (i < length(x) && cw[i + 1L] == cw[k]) i <- i + 1L
         return((x[j] + x[i]) / 2)
       }
       # "linear" with exact hit returns x[k]
       return(x[k])
     } else {
       # cw[k-1] < p < cw[k]  (k >= 1)
       if (k == 1L) {
         # all mass at the first point
         return(x[1L])
       }
       if (interpolation == "left")  return(x[k - 1L])
       if (interpolation == "right") return(x[k])
       if (interpolation == "midpoint") return((x[k - 1L] + x[k]) / 2)

       # linear interpolation in cumulative weight space
       w_below <- cw[k - 1L]
       w_above <- cw[k]
       t <- (p - w_below) / (w_above - w_below)  # in (0,1)
       return(x[k - 1L] + t * (x[k] - x[k - 1L]))
     }
   }

   # vectorize over probs
   res <- vapply(probs, get_q, numeric(1))
   names(res) <- paste0(probs)
   res
 }


 ######### MAIN FOREST HELPERS ###############

 # ========= Helper: turn names -> GRF args (subset-aligned) =========
 materialize_args <- function(.sd,model="rf",y_name=NULL, x_names=NULL, w_name=NULL, z_name=NULL,
                              forest_opts=list(), weight_col=NULL, cluster_col=NULL) {
   args <- list()
   if (!is.null(x_names)) {
     miss <- setdiff(x_names, names(.sd)); if (length(miss)) stop("Missing X: ", paste(miss, collapse=", "))
     args$X <- as.matrix(.sd[, ..x_names])
   }
   if (!is.null(y_name)) {
     if (!y_name %in% names(.sd)) stop("Missing Y: ", y_name)
     args$Y <- .sd[[y_name]]
     if (model %in% c("cf","iv")) {yhat <- paste0(y_name, ".hat"); if (yhat %in% names(.sd)) args$Y.hat <- .sd[[yhat]]}
   }
   if (!is.null(w_name)) {
     if (!w_name %in% names(.sd)) stop("Missing W: ", w_name)
     args$W <- .sd[[w_name]]
     if (model %in% c("cf","iv")) {what <- paste0(w_name, ".hat"); if (what %in% names(.sd)) args$W.hat <- .sd[[what]]}
   }
   if (!is.null(z_name)) {
     if (!z_name %in% names(.sd)) stop("Missing Z: ", z_name)
     args$Z <- .sd[[z_name]]
     if (model=="iv")  {zhat <- paste0(z_name, ".hat"); if (zhat %in% names(.sd)) args$Z.hat <- .sd[[zhat]]}
   }
   if (!is.null(weight_col))  { if (!weight_col  %in% names(.sd)) stop("Missing weight: ",  weight_col);  args$sample.weights <- .sd[[weight_col]] }
   if (!is.null(cluster_col)) { if (!cluster_col %in% names(.sd)) stop("Missing cluster: ", cluster_col); args$clusters       <- .sd[[cluster_col]] }
   c(args, forest_opts)
 }

