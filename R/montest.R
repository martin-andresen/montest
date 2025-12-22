  #' @param infile Path to the input file
  #' @return A matrix of the infile
  #' @export

  montest=function(data,D,Z,X=NULL,Y=NULL,W=NULL,treefraction=0.5,test=NULL,folds=2,onetest=TRUE,
                   normalize.Z=TRUE,normalize.pred=TRUE,stack=NULL,treetype="forest",
                   gridtypeY="equisized",gridtypeD="equisized",gridtypeZ="equisized",
                   Ysubsets = 4, Dsubsets = 4,Zsubsets=4,Y.res=TRUE,saveforest=F,min_n=1L,
                   weight=NULL,cluster=NULL,num.trees=2000,seed=10101,minsize=50,shrink=FALSE,shrink.alpha=0.3,
                   maxrankcp=5,prune=TRUE,cp=0,alpha=0.05,preselect="nonpositive", ##CART options
                   Zparameters=NULL,Yparameters=NULL,Qparameters=NULL,Dparameters=NULL,Cparameters=NULL,
                   tune.Qparameters="none",tune.Zparameters="none",tune.Cparameters="none",tune.Yparameters="none",tune.Dparameters="none",
                   tune.num.trees=200,tune.num.reps=50,tune.num.draws=1000,tunetype="one" ##tuning options
                   ){

    time=proc.time()
    set.seed(seed)


    ################### 1 CHECK INPUT #####################

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

    if (sum(test %in% c("AHS","MW"))==0&is.null(Y)==FALSE) Y=NULL

    ##Instance specific helper functions
    model_spec <- list(
      "simple" = list(model = "cf", xvars = c(X, W)),          # CF: Q ~ Z ; cov: X,W
      "MW"     = list(model = "rf", xvars = c(X, W, Y)),        # RF: Q ~ X,W,Y
      "BP"     = list(model = "cf", xvars = c(X, W)),           # CF: Q ~ Z ; cov: X,W
      "AHS"    = list(model = "cf", xvars = c(X, W, Y)),        # CF: Q ~ Z ; cov: X,W,Y
      "K"      = list(model = "iv", xvars = c(X, W))            # IV: Q ~ D (Z) ; cov: X,W
    )

    # Helper: fetch model spec per key and validate X columns
    get_spec <- function(key) {
      sp <- model_spec[[as.character(key)]]
      if (is.null(sp)) stop(sprintf("No model_spec entry for '%s'", key))
      miss <- setdiff(sp$xvars, names(data))
      if (length(miss)) stop(sprintf("Missing X columns for '%s': %s", key, paste(miss, collapse = ", ")))
      sp
    }

    ###################### 2 Prepare data #########################3
    XW=c(X,W)

    #
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

    if (is.null(cluster)==TRUE) {
      data[,sample:=1+1*(id_ %in% sample(n,size=floor(treefraction*n), replace = FALSE))]
    } else {
      s=unique(data[[cluster]])
      s2=sample(s,size=floor(0.5*length(s)),replace=FALSE)
      data[,sample:=1L+1L*(get(cluster)) %in% s2]
      rm(s,s2)
    }
    make_group_folds(data, K = folds, cluster_name = cluster)

    if ((is.null(cluster)==TRUE)&(is.null(stack)==FALSE)) cluster="id_"

    ############### 3 Discretize Z, D and Y into subsets ###############33
    ##TODO: ALLOW FOR WEIGHTED QUANTILES

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
            as.integer(cut(get(..Yno),breaks=c(-Inf,unique(weighted_quantile(get(..Yno),w=w,probs=(1:(Ysubsets-1))/Ysubsets)),Inf)))-1
          }]
        }
      } else {
        data[,(paste0(Yno,".bin")):=as.numeric(cut(get(..Yno),breaks=length(unique(get(..Yno)))))-1]
        }
      maxlevs=c(maxlevs,max(data[,get(..Yno)]))
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
    forest_optsC=list(num.trees=max(50,num.trees),tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps)

    ######################## 6a STACK DATA AND ESTIMATE Z.HAT / D.HAT / Q.HAT as early as possible #####
    if (stack==TRUE) {
      margins=c()


      ##HANDLE Z STACKS
        if (K>1) { ##Expand to all margins of Z
          times = 1L + 2L * as.integer(data[[Z]] > 0 & data[[Z]] < K)
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
            cluster_name = cluster,
            fold_col = "cf_fold",
            k_col = "cf_K"
          )
        } else {
          data[,(paste0(Z,".hat")):=(mean(get(..Z))*.N-get(..Z))/(.N-1),by=c("sample",margins)] ##leave one out mean
        }
        if (normalize.Z==TRUE) { #Normalize propensity scores
          data[, (paste0(Z, ".hat")) := get(paste0(Z, ".hat"))*.N / sum(get(Z)/get(paste0(Z,".hat"))), by = c("sample", margins)]
        }


      ##STACK MULTIPLE OUTCOMES
      if (L>1) {
        data=melt(data, measure = list(paste0(Y,rep(".bin",length(Y))),Y),value.name = c("Ystack.bin", "Ystack"),variable.name="outcome")
        data[,outcome:=Y,by=c("id",margins)]
        data[,maxlevsY:=maxlevs,by=c("id",margins)]
        margins=c(margins,"outcome")
        Y="Ystack"
      } else if (sum(test %in% c("BP","K"))>0) maxlevsY=maxlevs

      ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
      if ((sum(test %in% c("MW","AHS"))>0)&Y.res==TRUE)  {

        crossfit_hat(
          data,
          y_name = Y,
          x_names = X,
          margins = margins,
          weight_name = weight,
          cluster_name = cluster,
          fold_col = "cf_fold",
          k_col = "cf_K"
        )
        data[,paste0(Y,".res"):=Y-paste0(Y,".hat"),env=list(Y=Y)]

        }

      ##STACK MARGINS OF D
      if (J>1) { ##Stack data for all margins of D
          rows=rep(seq_len(nrow(data)),each=J)
          data= data[rows]    # replicate each row J times
          data[,dmargin:=seq_len(J),by=c("id_",margins)]
          #data=cbind(CJ(dmargin=(1:J),id=data$id),data[,!"id"])
          data[,D:=D>=dmargin,env=list(D=D)]
          margins=c(margins,"dmargin")
          }

      ##Estimate D.hat if test includes "K"
      if ("K" %in% test) {
        crossfit_hat(
          data,
          y_name = D,
          x_names = X,
          margins = margins,
          weight_name = weight,
          cluster_name = cluster,
          fold_col = "cf_fold",
          k_col = "cf_K"
        )
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
            condition="condition"
          } else {condition=test}

        ##Expand to equation=0,1 for BP, K and MW conditions
          if (sum(test %in% c("MW","BP","K"))>0) {
            data=data[rep(seq(.N),1+condition %in% c("BPK","MW","BP","K"))]
            data[condition %in% c("BPK","MW","BP","K"),equation:=seq(.N)-1*(condition %in% c("BPK","MW","BP","K")),by=c("id",margins)]
            margins=c(margins,"equation")
            }

        ##Expand to all groups of Ybin for BP, K conditions
           if (sum(test %in% c("BP","K"))>0) {
              data=data[rep(seq(.N), 1+(maxlevsY-1)*(condition %in% c("BPK","K","BP")))]
              data[condition %in% c("BP","K","BPK"),ybin:=(1:maxlevsY),by=c("id",margins)]
              margins=c(margins,"ybin")
           }

        ##Create outcome variable Q in stacked data
        if (sum(test %in% c("simple","AHS"))>0) data[condition %in% c("simple","AHS"),Q:=D,env=list(D=D)]
        if (sum(test %in% c("BPK","K","BP"))>0) data[condition %in% c("BPK","BP","K"),Q:=equation*(D*(name==ybin))-(1-equation)*(1-D)*(name==ybin),env=list(D=D,name=paste0(Y,".bin"))]
        if ("MW" %in% test) data[condition=="MW",Q:=equation*((1-Z.hat)*D*Z-Z.hat*D*(1-Z))+(1-equation)*(Z.hat*(1-D)*(1-Z)-(1-Z.hat)*(1-D)*Z),env=list(Z=Z,D=D)]
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
          crossfit_hat(
            data[condition %in% c("simple","BP","K","BPK")],
            y_name = "Q",
            x_names = X,
            margins = margins,
            weight_name = weight,
            cluster_name = cluster,
            fold_col = "cf_fold",
            k_col = "cf_K"
          )

        }

        if ("AHS" %in% test) {
          crossfit_hat(
            data[condition=="AHS"],
            y_name = "Q",
            x_names = c(X,paste0(Y,".res")),
            margins = margins,
            weight_name = weight,
            cluster_name = cluster,
            fold_col = "cf_fold",
            k_col = "cf_K"
          )
          }

        ##Split BP and K conditions if doing both
        if ((sum(test %in% c("BP","K"))>0)) {
            if (sum(test %in% c("BP","K"))==2) data=data[rep(seq(.N),1+1*(condition=="BPK"))]
            data[condition=="BPK",condition:=test[test %in% c("BP","K")],by=c("id",margins)]
        }

    ########## ESTIMATE ALL CAUSAL/REGRESSION/IV FORESTS AND  predict in/out of sample ##########
        # Housekeeping cols
        data[, idx__ := .I]
        data[, grp_id := .GRP, by = margins]

        fit_models(
          data      = data,
          condition = condition,      # name of column in `data`
          get_spec  = get_spec,
          materialize_args = materialize_args,
          variance=FALSE,cleanup=FALSE,
          wcol = weight, ccol = cluster,  Qcol = "Q", Zcol = Z, Dcol = D
          )

        ##Normalize pred, scores, pred_out
        if (normalize.pred==TRUE) {
          data[,scale:=sd(scores),by=c(margins,"sample")]
          data[, (c("pred","pred_o","scores")) := lapply(.SD, function(x) x / sqrt(scale)), .SDcols = c("pred","pred_o","scores") ]
          data[,scale:=NULL]
        }

        ######################################## FIND OPTIMAL SUBSET TO TEST AND TEST IN OPPOSITE SAMPLE #####################
        res=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,onetest=onetest)

        ##Global as comp
        if (is.null(cluster)==FALSE) cl=as.formula(paste0("~",cluster)) else cl=NULL
        if (is.null(weight)==FALSE) wg=as.formula(paste0("~",weight)) else wg=NULL
        global=data[,.(.N,uniqueN(id_)),by=sample]
        global=cbind(global,feols(scores~i(sample)-1,data=data,cluster=cl,weight=wg)$coeftable[,-4])
        colnames(global)=c("sample","N","G","coef","stderr","t")

        w_vec <- if (is.null(weight)) rep(1, nrow(data)) else data[[weight]]
        cols <- get0("XW", inherits = TRUE)

        if (onetest==TRUE) byv=NULL else byv=sample
        Xmeans=data[pred_o <= tau, {
          w <- ..w_vec[.I]
          lapply(.SD, weighted.mean, w = w, na.rm = TRUE)
        }, by = byv, .SDcols = cols]

        Xmeans_all=data[, {
          w <- ..w_vec[.I]
          lapply(.SD, weighted.mean, w = w, na.rm = TRUE)
        }, by = byv, .SDcols = cols]

        XSD=data[pred_o <= tau, {
          w <- ..w_vec[.I]
          lapply(.SD, w_sd, w = w, na.rm = TRUE)
        }, by = byv, .SDcols = cols]

        ###TODO: Y res means if MW,AHS
        if (sum(test %in% c("AHS","MW"))>0) { ##Add Y or Y to Xmeans
          data[pred_o<=tau,mean(Y),by=sample,env=list(Y=Y)]
        }
        ##SHOULD THESE ALSO BE WEIGHTED? Z==0 and Z==K...?
        marginsmat_all=data[,.N,by=c("sample",margins)]
        marginsmat=data[pred_o<=tau,.N,by=c("sample",margins)]
        marginsmat=merge(marginsmat,marginsmat_all, by=c("sample",margins),all=TRUE)
        marginsmat[,share_all:=N.y/sum(N.y), by=sample]
        marginsmat[,share:=ifelse(is.na(N.x),0,N.x)/sum(N.x,na.rm=TRUE), by=sample]
        marginsmat[,c("N.x","N.y"):=NULL]
    }

    ######################## 6b loop over margins and estimate CART or FOREST approach ################
    else {

        ##PLEASE IGNORE THIS BRANCH, NOT DONE!

        ########## PRE-ESTIMATE Z.hat to avoid redoing it for each loop

        for (z in 1:K) {
          data[Z==z|Z==z-1,paste0(Z,".hat",z):=do.call(regression_forest, materialize_args(.SD,model="rf",
             y_name=..Z, x_names=..X,forest_opts=c(forest_opts,Zparameters),
             weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
            cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
            $prediction,by=sample,env=list(Z=Z)]
        }
        ########## FIGURE OUT WHAT TO LOOP OVER ################
        margins=c();list=list();index=data.frame()
        if (K>1) { ##margins of Z
          margins=c(margins,"zmargin")
          list=append(list,list(zmargin=1:K))
        }

        if (J>1) { ##margins of D
          margins=c(margins,"dmargin")
          list=append(list,list(dmargin=1:J))
        }

        if (L>1) { ##Outcomes
          margins=c(margins,"outcome")
          list=append(list,list(outcome=Y))
        }

        if (length(test)>1) { ##testable conditions
          margins=c(margins,"condition")
          list=append(list,list(condition=test))
        }

        if (length(list)>0) index=do.call(CJ,list)

        if (sum(test %in% c("K","BP"))>0) { ##bins of Y
          margins=c(margins,"ybin")
          if (length(list)==0) {
            index=data.table(ybin=0:(maxlevs-1))
          } else {
            index=index[rep(seq(.N),ifelse(condition %in% c("K","BP"),rep(maxlevs,each=length(Y)),1))]
            index[condition %in% c("K","BP"),ybin:=sequence(maxlevs),by=margins[!margins %in% "outcome"]]
          }
        }

        if (sum(test %in% c("MW","K","BP"))>0) { #equation for two-eq. conditions
          margins=c(margins,"equation")
          if (nrow(index)==0) index=data.table(equation=0:1) else {
            if (length(test)>1) {
              index=index[rep(seq(.N),1+(condition %in% c("MW","BP","K")))]
              index[condition %in% c("MW","BP","K"),equation:=seq(.N)-1*(condition %in% c("MW","BP","K")),by=margins]
            } else {
              index=index[rep(seq(.N),2)]
              index[,equation:=seq(.N)-1,by=margins]
            }
          }
        }

        if (is.null(margins)==FALSE) {
          setorderv(index,cols=margins)
          setcolorder(index,margins)
          loopmax=max(nrow(index),1)
        } else loopmax=1

        ######## RUN LOOP ##########

        data[, idx__ := .I]
        data[, grp_id := .GRP]

      for (l in 1:loopmax) {
        if (K>1) z=index[l,zmargin] else z=1
        if (J>1) d=index[l,dmargin] else d=1
        if (sum(c("BP","K","MW") %in% test)>0) eq=index[l,equation] else a=NA
        if (sum(c("BP","K") %in% test)>0) y=index[l,ybin] else y=NA
        if (length(test)>1) cond=index[l,condition] else cond=test
        if (L>1) outcome=index[l,outcome] else outcome=Y

        ##replace Z
        data[,paste0(Z,".hat"):=paste0(Z,".hat",z),env=list(Z=Z)]

        ##define Q
        if (cond %in% c("simple","AHS")) data[Z==z|Z==z-1,Q:=D>=d,env=list(D=D)]
        else if (cond=="MW") data[Z==z|Z==z-1,Q:=a*((1-paste0(Z,".hat"))*(D>=d)*(Z>=z)-paste0(Z,".hat")*(D>=d)*(1-(Z>=z)))+(1-a)*(paste0(Z,".hat")*(1-(D>=d))*(1-(Z>=z))-(1-paste0(Z,".hat"))*(1-(D>=d))*(Z>=z)),env=list(Z=Z,D=D)]
        else data[Z==z|Z==z-1,Q:=a*((D>=d)*(name==y))-(1-a)*(1-(D>=d))*(name==y),env=list(D=D,name=paste0(Y,".bin"))]

        ##Estimate Q.hat
        if (cond %in% c("MW","AHS")) {
          data[Z==z|Z==z-1,Q.hat:=do.call(regression_forest, materialize_args(.SD,model="rf",
            y_name="Q", x_names=..XWY,forest_opts=c(forest_opts,Qparameters),
            weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
            cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
            $prediction,by=sample]
        }
        else if (cond %in% c("BP","K")) {
          data[Z==z|Z==z-1,Q.hat:=do.call(regression_forest, materialize_args(.SD,model="rf",
             y_name="Q", x_names=..XW,forest_opts=c(forest_opts,Qparameters),
             weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
             cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
             $prediction,by=sample]
        }

        if (cond=="K") { ##estimate D.hat
          data[Z==z|Z==z-1,paste0(D,".hat"):=do.call(regression_forest, materialize_args(.SD,model="rf",
            y_name=..D, x_names=..XW,forest_opts=c(forest_opts,Dparameters),
            weight_col=get0("weight",  inherits = TRUE, ifnotfound = NULL),
           cluster_col=get0("cluster",  inherits = TRUE, ifnotfound = NULL)))
               $prediction,by=sample,env=list(D=D)]
        }

        ##Estimate causal forests and predict scores & predicted effects
        if (K>1) idx <- which(which(data[,Z]==z|data[,Z]==z-1)) else idx=NULL
        fit_models(data,
              condition = cond,
              get_spec = get_spec,
              materialize_args = materialize_args,
              wcol = weight, ccol = cluster, var=shrink,
              rows = idx
          )
        }



      }

    ################ 7: Multiple hypothesis testing and output #####################
      if (onetest==TRUE) {
        res[is.na(t)==FALSE,p.raw:=pnorm(t)]
        minsample=NA
        minp=res[is.na(sample)==TRUE,p.raw]
        } else {
          res[is.na(t)==FALSE,p.raw:=pnorm(t)]
          for (m in c("holm","hochberg","BH","BY")) {
            res[train==FALSE&is.na(t)==FALSE,paste0("p.",m):=p.adjust(p.raw,method=m)]
          }
          minsample=res[train==FALSE&which.min(p.raw),sample]
          if (length(minsample)>1) minsample=minsample[1]
          minp=apply(res[train==FALSE&is.na(t)==FALSE,c("p.raw","p.holm","p.hochberg","p.BH","p.BY")],2,min)
      }

      global[,p.raw:=pnorm(t)]
      if (onetest==FALSE) {
        for (m in c("holm","hochberg","BH","BY")) {
        global[,paste0("p.",m):=p.adjust(p.raw,method=m)]
        }
      }

      time=proc.time()-time
      return=list(results=res,global=global,minsample=minsample,minp=minp,time=time,margins=marginsmat,Xmeans_all=Xmeans_all,XSD=XSD,Xmeans=Xmeans)
      return(return)
  }


##################################################### HELPER FUNCTIONS ######################################

# THIS FUNCTION ESTIMATES A TREE USING CART, TEST IN EACH NODE, DETERMINE relevance and perform the test in the other sample
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

 weighted_quantile <- function(x, w = NULL, probs = c(0.25, 0.5, 0.75),
                               na.rm = TRUE,
                               interpolation = c("linear", "left", "right", "midpoint")) {
   interpolation <- match.arg(interpolation)

   # default: no weights -> equal weights
   if (is.null(w)) w <- rep(1, length(x))

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
         j <- k
         while (j > 1L && cw[j - 1L] == cw[k]) j <- j - 1L
         i <- k
         while (i < length(x) && cw[i + 1L] == cw[k]) i <- i + 1L
         return((x[j] + x[i]) / 2)
       }
       return(x[k]) # linear or default on exact hit
     } else {
       # cw[k-1] < p < cw[k]
       if (k == 1L) return(x[1L])
       if (interpolation == "left")      return(x[k - 1L])
       if (interpolation == "right")     return(x[k])
       if (interpolation == "midpoint")  return((x[k - 1L] + x[k]) / 2)

       # linear interpolation in cumulative weight space
       w_below <- cw[k - 1L]
       w_above <- cw[k]
       t <- (p - w_below) / (w_above - w_below)
       return(x[k - 1L] + t * (x[k] - x[k - 1L]))
     }
   }

   # vectorize over probs
   res <- vapply(probs, get_q, numeric(1))
   names(res) <- paste0(probs)
   res
 }

 w_sd <- function(x, w = NULL, na.rm = TRUE, unbiased = FALSE) {
   if (is.null(w)) return(sd(x, na.rm = na.rm))

   if (na.rm) {
     keep <- is.finite(x) & is.finite(w)
     x <- x[keep]; w <- w[keep]
   }
   if (!length(x)) return(NA_real_)

   mu <- weighted.mean(x, w)
   v  <- sum(w * (x - mu)^2) / sum(w)

   if (unbiased && sum(w) > 1)
     v <- v * sum(w) / (sum(w) - 1)

   sqrt(v)
 }




 ######### MAIN FOREST HELPERS ###############

 # ========= Helper: Randomize into K folds =========

 make_group_folds <- function(DT,
                              K = 10L,
                              cluster_name = NULL,
                              fold_col = "cf_fold",
                              diag_prefix = "cf_",
                              verbose = TRUE) {
   stopifnot(is.data.table(DT))
   stopifnot(is.numeric(K) && K >= 2)

   k_col <- paste0(diag_prefix, "K")
   g_col <- paste0(diag_prefix, "G")
   n_col <- paste0(diag_prefix, "n_ok")

   DT[, c(fold_col, k_col, g_col, n_col) := {

     ok <- rep(TRUE, .N)
     n_ok <- sum(ok)

     if (!is.null(cluster_name)) {
       cl <- get(cluster_name)
       ok <- ok & !is.na(cl)
       n_ok <- sum(ok)

       if (n_ok == 0L) {
         return(list(rep(NA_integer_, .N), 0L, 0L, 0L))
       }

       cl_ok <- cl[ok]
       ucl <- unique(cl_ok)
       G <- length(ucl)
       Kg <- as.integer(min(K, G))

       fold <- rep(NA_integer_, .N)
       if (Kg < 2L) {
         fold[ok] <- 1L
         return(list(fold, Kg, G, n_ok))
       }

       fold_cl <- sample(rep_len(seq_len(Kg), G))
       fold[ok] <- as.integer(fold_cl[match(cl_ok, ucl)])

       list(fold, Kg, G, n_ok)
     } else {
       if (n_ok == 0L) {
         return(list(rep(NA_integer_, .N), 0L, NA_integer_, 0L))
       }

       Kg <- as.integer(min(K, n_ok))
       fold <- rep(NA_integer_, .N)

       if (Kg < 2L) {
         fold[ok] <- 1L
         return(list(fold, Kg, NA_integer_, n_ok))
       }

       fold[ok] <- as.integer(sample(rep_len(seq_len(Kg), n_ok)))

       list(fold, Kg, NA_integer_, n_ok)
     }
   }, by = "sample"]

   if (verbose) {
     bad <- DT[get(k_col) < 2L, unique(sample)]
     if (length(bad) > 0L) {
       warning(sprintf(
         "Cross-fitting folds: %d sample(s) have %s < 2 (too few clusters/rows). Using 1 fold in those samples.",
         length(bad), k_col
       ))
     } else {
       message(sprintf(
         "Cross-fitting folds: created %s with up to %d folds within each sample.",
         fold_col, K
       ))
     }
   }

   invisible(DT)
 }



 # ========= Helper: Predict nuissance using regression-forest - cross-fit =========


 crossfit_hat <- function(DT,
                          y_name,
                          x_names,
                          margins = NULL,            # default NULL
                          weight_name = NULL,
                          cluster_name = NULL,
                          fold_col = "cf_fold",
                          k_col = "cf_K",
                          forest_opts = list(),
                          hat_suffix = ".hat",
                          diag_prefix = "nh_",
                          verbose = TRUE,
                          residuals = FALSE) {
   stopifnot(is.data.table(DT))
   stopifnot(is.character(y_name) && length(y_name) == 1L)
   stopifnot(is.character(x_names) && length(x_names) >= 1L)

   if (is.null(margins)) margins <- character() else margins <- as.character(margins)

   grp_cols <- c("sample", margins)
   if (!all(grp_cols %chin% names(DT))) {
     stop("Missing grouping columns: ", paste(setdiff(grp_cols, names(DT)), collapse = ", "))
   }
   if (!(fold_col %chin% names(DT))) stop("Missing fold column: ", fold_col)

   hat_col <- paste0(y_name, hat_suffix)
   diag_usedK <- paste0(diag_prefix, y_name, "_K_used")
   diag_nok   <- paste0(diag_prefix, y_name, "_n_ok")
   diag_note  <- paste0(diag_prefix, y_name, "_note")

   DT[, c(hat_col, diag_usedK, diag_nok, diag_note) := {
     # outcome (force numeric; map logical -> 0/1)
     y_raw <- get(y_name)
     y <- as.numeric(y_raw)
     if (is.logical(y_raw)) y <- y - 1

     fold <- get(fold_col)

     w  <- if (is.null(weight_name)) NULL else get(weight_name)
     cl <- if (is.null(cluster_name)) NULL else get(cluster_name)

     n <- .N
     preds <- rep(NA_real_, n)
     note <- ""
     n_ok <- 0L

     # K used (constant within group if provided; otherwise inferred)
     K_used <- if (!is.null(k_col) && (k_col %chin% names(DT))) as.integer(get(k_col)[1L]) else NA_integer_

     # IMPORTANT: build each ok-component explicitly as length-.N logical
     ok_y    <- is.finite(y)
     ok_fold <- is.finite(as.numeric(fold))          # robust even if fold is integer/factor-like
     ok_x    <- stats::complete.cases(.SD)           # .SD is DT subset, returns length .N
     ok_w    <- if (is.null(w))  rep(TRUE, n) else is.finite(w)
     ok_cl   <- if (is.null(cl)) rep(TRUE, n) else !is.na(cl)

     ok <- ok_y & ok_fold & ok_x & ok_w & ok_cl
     n_ok <- sum(ok)

     if (is.na(K_used)) K_used <- as.integer(uniqueN(fold[ok]))

     if (n_ok < 2L) {
       mu <- if (is.null(w)) mean(y, na.rm = TRUE) else weighted.mean(y, w = w, na.rm = TRUE)
       preds[ok_y] <- mu
       note <- "fallback_mean_too_few_ok"
     } else {
       y_ok <- y[ok]

       if (length(unique(y_ok)) < 2L) {
         preds[ok] <- y_ok[1L]
         note <- "fallback_constant_y"
       } else if (!is.finite(K_used) || K_used < 2L) {
         w_ok <- if (is.null(w)) NULL else w[ok]
         mu <- if (is.null(w_ok)) mean(y_ok, na.rm = TRUE) else weighted.mean(y_ok, w = w_ok, na.rm = TRUE)
         preds[ok] <- mu
         note <- "fallback_mean_K<2"
       } else {
         X_ok   <- as.matrix(.SD[ok])                 # now safe: 2D matrix
         f_ok   <- as.integer(fold[ok])
         ok_idx <- which(ok)

         w_ok  <- if (is.null(w)) NULL else w[ok]
         cl_ok <- if (is.null(cl)) NULL else cl[ok]

         for (k in sort(unique(f_ok))) {
           te <- f_ok == k
           tr <- !te
           if (!any(te) || !any(tr)) next

           rf <- do.call(
             regression_forest,
             c(
               list(
                 X = X_ok[tr, , drop = FALSE],
                 Y = y_ok[tr],
                 sample.weights = if (is.null(w_ok)) NULL else w_ok[tr],
                 clusters       = if (is.null(cl_ok)) NULL else cl_ok[tr]
               ),
               forest_opts
             )
           )

           preds[ok_idx[te]] <- predict(rf, X_ok[te, , drop = FALSE])$predictions
         }

         note <- "ok_crossfit"
       }
     }

     list(preds, K_used, n_ok, note)
   }, by = c("sample", margins), .SDcols = x_names]

   if (verbose) {
     bad <- DT[get(diag_note) != "ok_crossfit",
               unique(.SD),
               .SDcols = c("sample", margins, diag_note, diag_usedK, diag_nok)]
     if (nrow(bad) > 0L) {
       warning(sprintf("crossfit_hat(%s): %d group(s) used fallback.", y_name, nrow(bad)))
     }
   }

   invisible(DT)
 }
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
     args$Y <- as.numeric(.sd[[y_name]])
     if (model %in% c("cf","iv")) {yhat <- paste0(y_name, ".hat"); if (yhat %in% names(.sd)) args$Y.hat <- as.numeric(.sd[[yhat]])}
   }
   if (!is.null(w_name)) {
     if (!w_name %in% names(.sd)) stop("Missing W: ", w_name)
     args$W <- as.numeric(.sd[[w_name]])
     if (model %in% c("cf","iv")) {what <- paste0(w_name, ".hat"); if (what %in% names(.sd)) args$W.hat <- as.numeric(.sd[[what]])}
   }
   if (!is.null(z_name)) {
     if (!z_name %in% names(.sd)) stop("Missing Z: ", z_name)
     args$Z <- as.numeric(.sd[[z_name]])
     if (model=="iv")  {zhat <- paste0(z_name, ".hat"); if (zhat %in% names(.sd)) args$Z.hat <- as.numeric(.sd[[zhat]])}
   }
   if (!is.null(weight_col))  { if (!weight_col  %in% names(.sd)) stop("Missing weight: ",  weight_col);  args$sample.weights <- as.numeric(.sd[[weight_col]]) }
   if (!is.null(cluster_col)) { if (!cluster_col %in% names(.sd)) stop("Missing cluster: ", cluster_col); args$clusters       <- as.numeric(.sd[[cluster_col]]) }
   c(args, forest_opts)
 }


 fit_models <- function(
    data,
    condition,                 # string: column name in data OR a fixed tag ("rf","cf","iv")
    get_spec,
    materialize_args,
    wcol = NULL, ccol = NULL,
    Qcol = "Q", Zcol = "Z", Dcol = "D",
    rows = NULL,               # integer row indices to process; if NULL, do all rows
    cross_within_subset = TRUE,# if TRUE, write pred_o only to rows in `rows`
    variance = TRUE,  # <- NEW: also compute pred_var / pred_o_var via GRF
    cleanup = TRUE
 ) {
   stopifnot(data.table::is.data.table(data))

   if (is.null(rows)) rows <- seq_len(nrow(data))
   if (!length(rows)) return(invisible(data))

   # stable ids for write-back
   if (!(".__rowid__" %chin% names(data))) {
     data[, .__rowid__ := .I]
     on.exit(try(data[, .__rowid__ := NULL], silent = TRUE), add = TRUE)
   }

   # work view on the subset only
   DT <- data[rows]
   DT[, idx__ := .__rowid__]

   # ensure grouping columns exist
   if (!("grp_id" %chin% names(DT))) DT[, grp_id := 1L]
   if (!("sample" %chin% names(DT))) DT[, sample := 1L]

   cond_arg <- condition

   DT[, {
     # decide key for this group
     key <- if (is.character(..cond_arg) && length(..cond_arg) == 1L) {
       if ((..cond_arg) %chin% names(.SD)) .SD[[..cond_arg]][1L] else ..cond_arg
     } else stop("`condition` must be a single string (column name or fixed tag).")

     sp  <- get_spec(key)

     idx <- idx__
     gid <- grp_id[1L]
     s0  <- sample[1L]

     # ---- fit according to spec ----
     if (sp$model == "rf") {
       args <- materialize_args(
         .sd         = .SD,
         model       = "rf",
         y_name      = Qcol,
         x_names     = sp$xvars,
         forest_opts = c(get0("forest_optsC", inherits=TRUE, ifnotfound=list()),
                         get0("Cparameters",  inherits=TRUE, ifnotfound=list())),
         weight_col  = wcol,
         cluster_col = ccol
       )
       fit <- do.call(grf::regression_forest, args)
       pin <- predict(fit, estimate.variance = variance)
       pred_in <- as.numeric(if (is.list(pin)) pin$predictions else pin)
       pred_in_var <- if (variance && is.list(pin) && !is.null(pin$variance.estimates))
         as.numeric(pin$variance.estimates) else rep(NA_real_, length(idx))
       scores_in <- .SD[[Qcol]]  # RF: scores = outcome

     } else if (sp$model == "cf") {
       args <- materialize_args(
         .sd         = .SD,
         model       = "cf",
         y_name      = Qcol,
         w_name      = Zcol,
         x_names     = sp$xvars,
         forest_opts = c(get0("forest_optsC", inherits=TRUE, ifnotfound=list()),
                         get0("Cparameters",  inherits=TRUE, ifnotfound=list())),
         weight_col  = wcol,
         cluster_col = ccol
       )
       fit <- do.call(grf::causal_forest, args)
       pin <- predict(fit, estimate.variance = variance)
       pred_in <- as.numeric(if (is.list(pin)) pin$predictions else pin)
       pred_in_var <- if (variance && is.list(pin) && !is.null(pin$variance.estimates))
         as.numeric(pin$variance.estimates) else rep(NA_real_, length(idx))
       scores_in <- as.numeric(get_scores(fit))

     } else if (sp$model == "iv") {
       args <- materialize_args(
         .sd         = .SD,
         model       = "iv",
         y_name      = Qcol,
         w_name      = Dcol,
         z_name      = Zcol,
         x_names     = sp$xvars,
         forest_opts = c(get0("forest_optsC", inherits=TRUE, ifnotfound=list()),
                         get0("Cparameters",  inherits=TRUE, ifnotfound=list())),
         weight_col  = wcol,
         cluster_col = ccol
       )
       fit <- do.call(grf::instrumental_forest, args)
       pin <- predict(fit, estimate.variance = variance)
       pred_in <- as.numeric(if (is.list(pin)) pin$predictions else pin)
       pred_in_var <- if (variance && is.list(pin) && !is.null(pin$variance.estimates))
         as.numeric(pin$variance.estimates) else rep(NA_real_, length(idx))
       scores_in <- as.numeric(get_scores(fit))

     } else stop(sprintf("Unknown model '%s'", sp$model))

     # ---- write back (only this subset) ----
     data.table::set(data, i = idx, j = "pred",     value = pred_in)
     data.table::set(data, i = idx, j = "scores",   value = scores_in)
     if (variance) {
       data.table::set(data, i = idx, j = "pred_var", value = pred_in_var)
     }

     # ---- cross-sample predictions (pred_o, pred_o_var) ----
     # decide which "other" rows to fill
     if (cross_within_subset) {
       other_idx <- data[.__rowid__ %in% rows & grp_id == gid & sample != s0, which = TRUE]
     } else {
       other_idx <- data[grp_id == gid & sample != s0, which = TRUE]
     }

     if (length(other_idx)) {
       X_other <- as.matrix(data[other_idx, sp$xvars, with = FALSE])
       pout <- predict(fit, X_other, estimate.variance = variance)
       pred_o <- as.numeric(if (is.list(pout)) pout$predictions else pout)
       data.table::set(data, i = other_idx, j = "pred_o", value = pred_o)
       if (variance && is.list(pout) && !is.null(pout$variance.estimates)) {
         pred_o_var <- as.numeric(pout$variance.estimates)
         data.table::set(data, i = other_idx, j = "pred_o_var", value = pred_o_var)
       } else if (variance) {
         data.table::set(data, i = other_idx, j = "pred_o_var", value = NA_real_)
       }
     }

     NULL
   }, by = .(grp_id, sample)]

   if (cleanup) data[, c("idx__","grp_id") := NULL]
   invisible(data)
 }




 forest_test <- function(
    data,
    cluster = NULL,         # string: name of the cluster column (accept NULL?)
    weight = NULL,           # NULL, scalar, or string (weight column)
    sample = "sample",
    pred   = "pred",
    pred_o = "pred_o",
    scores = "scores",
    minsize = 50L,
    onetest=FALSE
 ) {

   # Order by training score
   setorderv(data, cols = c(sample,pred))

   # N per sample (running index in current order)
   data[, N := seq_len(.N), by = sample]

   # Within-cluster running totals (per current order)
   if (is.null(weight)==FALSE) {
     data[, `:=`(a = weight * scores, b = weight),env=list(weight=weight,scores=scores)]
     data[, `:=`(WgY = cumsum(a), Wg = cumsum(b)), by = c(sample, cluster)]
   }
   else {
     data[, `:=`(a = scores, b = N),env=list(scores=scores)]
     data[, `:=`(WgY = cumsum(scores), Wg = N), by = c(sample, cluster)]
   }
   # Global running totals and running mean by sample
   data[, `:=`(SW = cumsum(b), SWY = cumsum(a)), by = sample]
   data[, m := SWY / SW, by = sample]

   # Per-row deltas (avoid shifts)
   data[, `:=`(
     dTA2 =  WgY^2 - (WgY - a)^2,
     dTB2 =  Wg^2  - (Wg  - b)^2,
     dTAB =  WgY*Wg - (WgY - a)*(Wg - b)
   )]

   # Cumulative aggregates across rows (by sample)
   data[, `:=`(TA2 = cumsum(dTA2), TB2 = cumsum(dTB2), TAB = cumsum(dTAB)), by = sample]

   # Number of unique clusters seen so far (by sample, current order)
   data[, G := cumsum(!duplicated(cluster)), by = sample,env=list(cluster=cluster)]

   # CRV1 (small-sample adj.) SE of weighted mean
   data[, sumS2 := (TA2 - 2*m*TAB + (m^2)*TB2) / (SW^2)]
   data[, se := sqrt((G / pmax.int(G - 1L, 1L)) * sumS2)]
   data[G < 2, se := NA_real_]

   # t-stat for H0: mean = 0
   data[, t := m / se]

   # --- Dual-constraint tau ensuring >= minsize clusters in BOTH samples ---
   # Train-side cutoff by sample
   tau_tr <- data[G >= minsize, .(tau_tr = min(pred)), by = sample,env=list(pred=pred)]

   # Estimation-side cutoff by sample (order by pred_o, count clusters)
   setorderv(data, cols = c(sample, pred_o))
   tau_est <- data[cumsum(!duplicated(cluster)) >= minsize, .(tau_est = min(pred_o)), by = sample,env=list(cluster=cluster,pred_o=pred_o)]

   # Combine into cross-sample constraint: for two samples s1,s2
   tau=c(max(as.numeric(c(tau_tr[1,2],tau_est[2,2]))),max(as.numeric(c(tau_tr[2,2],tau_est[1,2]))))
   data[, tau := fifelse(sample==1,tau[1],tau[2]),env=list(sample=sample)]

   # --- Train-sample optimal cutoff: choose pred with minimizes t among pred > tau ---
   res <- data[pred >= tau&G>=minsize, .SD[which.min(t)], by = sample,
               .SDcols = c("G","N","m","se","t", pred),env=list(pred=pred)]
   res[,train:=TRUE]
   setcolorder(res,c("train","sample"))

   # Replace dual-constraint tau by the optimal cutoff from the *other* sample
   data[, tau := fifelse(sample==1,res[2,pred],res[1,pred]),env=list(sample=sample,pred=pred)]

   # ----- Final CR1 test on estimation sample: pred_o <= tau -----
   if (is.null(cluster)==FALSE) cl=as.formula(paste0("~",cluster)) else cl=NULL
   if (is.null(weight)==FALSE) wg=as.formula(paste0("~",weight)) else wg=NULL

   if (onetest==FALSE) {
     fit <- feols(as.formula(paste0(scores, " ~ i(", sample, ") - 1")),
        data = data[pred_o<=tau,env=list(pred_o=pred_o)], vcov = cl, weights = wg)
     byv=sample
   } else {
     fit <- feols(scores~1,
                  data = data[pred_o<=tau,env=list(pred_o=pred_o)], vcov = cl, weights = wg)
    byv=NULL
     }

   # Counts on estimation sample
   if (is.null(cluster)==FALSE) GN <- data[pred_o<=tau, .(G.est = uniqueN(cluster), N.est = .N), by = byv,env=list(pred_o=pred_o,cluster=cluster)]
   else GN <- data[pred_o<=tau, .(G.est = .N, N.est = .N), by = byv,env=list(pred_o=pred_o)]

   # Assemble output
   if (onetest==TRUE) test=c(FALSE,NA,GN,fit$coeftable[1:3],NA) else {
     test=as.data.table(cbind(GN,fit$coeftable[,1:3]))
     test[,train:=FALSE]
     test[,tau_cutoff:=NA]
     setcolorder(test,c("train","sample"))
   }
   res=rbind(res,test,use.names=FALSE)
   colnames(res)=c("train","sample","G","N","coef","stderr","t","tau cutoff")

   data[, c("N","a","b","WgY","Wg","dTA2","dTB2","dTAB","TA2","TB2","TAB","sumS2","se","G","t","m","SW","SWY") := NULL]
   return(res)

 }
