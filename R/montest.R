  #' @param infile Path to the input file
  #' @return A matrix of the infile
  #' @export

  montest=function(data,D,Z,X=NULL,Y=NULL,W=NULL,treefraction=0.5,test=NULL,inner_folds=5,inner_nuisance=FALSE,pool=TRUE,
                   normalize.Z=TRUE,stack=NULL,treetype="forest",aipw_clip=1e-3,gridpoints=NULL,
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

    ##OUTER SPLIT
    make_group_folds(data,K = 2,cluster_name = cluster, fold_col = "sample",verbose = TRUE)

    ##OPTIONAL INNER SPLIT
    if (is.null(inner_folds)==FALSE) {
      make_group_folds(data,K = inner_folds,cluster_name = cluster,fold_col = "cf_fold",verbose = TRUE,by_col="sample")
      foldname="cf_fold"
      } else {
      foldname=NULL
    }

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
            forest_opts = Zparameters
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
        data[,outcome:=Y,by=c("id",margins)]
        data[,maxlevsY:=maxlevs,by=c("id",margins)]
        margins=c(margins,"outcome")
        Y="Ystack"
      } else if (sum(test %in% c("BP","K"))>0) maxlevsY=maxlevs

      ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
      if ((sum(test %in% c("MW","AHS"))>0)&Y.res==TRUE)  {
        if (length(test)>1) i=which(data$condition %in% c("MW","AHS")) else i=NULL
        crossfit_hat(
          data,
          i=i,
          y_name = Y,
          x_names = X,
          margins = margins,
          weight_name = weight,
          forest_opts = Yparameters
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
        if (length(test)>1) i=which(data$condition=="K") else i=NULL
        crossfit_hat(
          data,
          i=i,
          y_name = D,
          x_names = X,
          margins = margins,
          weight_name = weight,
          forest_opts = Zparameters
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
          if (length(test)>1) i=which(data$condition %in%  c("simple","BP","K","BPK")) else i=NULL
          crossfit_hat(
            data,
            i=i,
            y_name = "Q",
            x_names = X,
            margins = margins,
            weight_name = weight,
            forest_opts = Qparameters
          )
        }

        if ("AHS" %in% test) {
          if (length(test)>1) i=which(data$condition=="AHS") else i=NULL
          crossfit_hat(
            data,
            i=i,
            y_name = "Q",
            x_names = c(X,paste0(Y,".res")),
            margins = margins,
            weight_name = weight,
            forest_opts = Qparameters
          )
          }

        ##Split BP and K conditions if doing both
        if ((sum(test %in% c("BP","K"))>0)) {
            if (sum(test %in% c("BP","K"))==2) data=data[rep(seq(.N),1+1*(condition=="BPK"))]
            data[condition=="BPK",condition:=test[test %in% c("BP","K")],by=c("id",margins)]
        }

    ########## ESTIMATE ALL CAUSAL/REGRESSION/IV FORESTS AND  predict in/out of sample ##########

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
              forest_opts = Cparameters,
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
                   forest_opts = Cparameters,
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
                   forest_opts = Cparameters,
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
                   forest_opts = Cparameters,
                   aipw_clip=aipw_clip)
      }


        ######################################## FIND OPTIMAL SUBSET TO TEST AND TEST IN OPPOSITE SAMPLE #####################
        res=forest_test(data,cluster=cluster,weight=weight,minsize=minsize,pool=pool,gridpoints=gridpoints)

        ##Global as comp
        if (is.null(cluster)==FALSE) cl=as.formula(paste0("~",cluster)) else cl=NULL
        if (is.null(weight)==FALSE) wg=as.formula(paste0("~",weight)) else wg=NULL
        global=data[,.(.N,uniqueN(id_)),by=sample]
        global=cbind(global,feols(scores~i(sample)-1,data=data,cluster=cl,weight=wg)$coeftable[,-4])
        colnames(global)=c("sample","N","G","coef","stderr","t")

        w_vec <- if (is.null(weight)) rep(1, nrow(data)) else data[[weight]]
        cols <- get0("XW", inherits = TRUE)

        if (pool==TRUE) byv=NULL else byv="sample"
        data[sample==1,tau:=res[2,"tau cutoff"]]
        data[sample==2,tau:=res[1,"tau cutoff"]]
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
      if (pool==TRUE) {
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
      if (pool==FALSE) {
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
                              fold_col = "fold",
                              by_col = NULL,        # NULL = global; else within groups
                              diag_prefix = "cf_",
                              verbose = TRUE) {
   stopifnot(data.table::is.data.table(DT))
   stopifnot(is.numeric(K) && K >= 2)
   stopifnot(is.character(fold_col), length(fold_col) == 1L)

   if (!is.null(by_col)) {
     stopifnot(is.character(by_col), length(by_col) == 1L, by_col %chin% names(DT))
   }
   if (!is.null(cluster_name)) {
     stopifnot(is.character(cluster_name), length(cluster_name) == 1L, cluster_name %chin% names(DT))
   }

   has_diag <- !is.null(diag_prefix)
   if (has_diag) {
     k_col <- paste0(diag_prefix, "K")
     g_col <- paste0(diag_prefix, "G")
     n_col <- paste0(diag_prefix, "n_ok")
   }

   # helper: create folds for one group
   make_one <- function(n, cl = NULL) {
     ok <- rep(TRUE, n)
     n_ok <- n

     fold <- rep(NA_integer_, n)

     if (!is.null(cl)) {
       ok <- ok & !is.na(cl)
       n_ok <- sum(ok)

       if (n_ok == 0L) {
         if (has_diag) return(list(fold, 0L, 0L, 0L))
         return(fold)
       }

       cl_ok <- cl[ok]
       ucl <- unique(cl_ok)
       G <- length(ucl)
       Kg <- as.integer(min(K, G))

       if (Kg < 2L) {
         fold[ok] <- 1L
         if (has_diag) return(list(fold, Kg, G, n_ok))
         return(fold)
       }

       fold_cl <- sample(rep_len(seq_len(Kg), G))
       fold[ok] <- as.integer(fold_cl[match(cl_ok, ucl)])

       if (has_diag) return(list(fold, Kg, G, n_ok))
       return(fold)
     }

     # no clustering
     if (n_ok == 0L) {
       if (has_diag) return(list(fold, 0L, NA_integer_, 0L))
       return(fold)
     }

     Kg <- as.integer(min(K, n_ok))
     if (Kg < 2L) {
       fold[ok] <- 1L
       if (has_diag) return(list(fold, Kg, NA_integer_, n_ok))
       return(fold)
     }

     fold[ok] <- as.integer(sample(rep_len(seq_len(Kg), n_ok)))

     if (has_diag) return(list(fold, Kg, NA_integer_, n_ok))
     fold
   }

   if (is.null(by_col)) {
     # ---------- GLOBAL SPLIT ----------
     n <- nrow(DT)
     cl <- if (!is.null(cluster_name)) DT[[cluster_name]] else NULL

     if (has_diag) {
       DT[, c(fold_col, k_col, g_col, n_col) := make_one(n, cl)]
     } else {
       DT[, (fold_col) := make_one(n, cl)]
     }

     if (verbose) {
       message(sprintf(
         "Cross-fitting folds: created '%s' with up to %d folds (global).",
         fold_col, K
       ))
     }

   } else {
     # ---------- WITHIN-GROUP SPLIT ----------
     if (has_diag) {
       DT[, c(fold_col, k_col, g_col, n_col) := {
         cl <- if (!is.null(cluster_name)) get(cluster_name) else NULL
         make_one(.N, cl)
       }, by = by_col]
     } else {
       DT[, (fold_col) := {
         cl <- if (!is.null(cluster_name)) get(cluster_name) else NULL
         make_one(.N, cl)
       }, by = by_col]
     }

     if (verbose) {
       message(sprintf(
         "Cross-fitting folds: created '%s' with up to %d folds within '%s'.",
         fold_col, K, by_col
       ))
     }
   }

   invisible(DT)
 }





 # ========= Helper: Predict nuissance using regression-forest - cross-fit across sample parts =========
 crossfit_hat <- function(DT,
                               i = NULL,
                               y_name,
                               x_names,
                               margins = NULL,
                               weight_name = NULL,
                               forest_opts = list(),
                               hat_suffix = ".hat",
                               verbose = FALSE) {

   stopifnot(data.table::is.data.table(DT))
   stopifnot(is.character(y_name), length(y_name) == 1L)
   stopifnot(is.character(x_names), length(x_names) >= 1L)
   stopifnot("sample" %chin% names(DT))
   stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

   if (is.null(i)) i <- DT[, .I]
   stopifnot(is.integer(i) || is.numeric(i))
   i <- as.integer(i)
   n_i <- length(i)
   if (n_i == 0L) return(invisible(DT))

   if (is.null(margins)) margins <- character()
   grp_cols <- as.character(margins)

   missing_cols <- setdiff(
     c("sample", y_name, x_names, grp_cols, if (!is.null(weight_name)) weight_name),
     names(DT)
   )
   if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

   hat_col <- paste0(y_name, hat_suffix)
   if (!(hat_col %chin% names(DT))) DT[, (hat_col) := NA_real_]

   # ---- Precompute X for rows i ONCE ----
   X_all <- as.matrix(DT[i, ..x_names])

   # rowid within i subset for fast X indexing
   rid_col <- ".__rid__"
   DT[i, (rid_col) := seq_len(n_i)]
   getX <- function(idx) X_all[DT[[rid_col]][idx], , drop = FALSE]

   # pull vectors once
   sample_all <- DT[["sample"]]
   y_all <- DT[[y_name]]
   w_all <- if (!is.null(weight_name)) DT[[weight_name]] else NULL

   # output buffer aligned with i
   preds_buf <- rep(NA_real_, n_i)

   # map DT-row indices -> position in buffer (safer than integer(max(i)))
   pos_map <- function(idx) match(idx, i)

   # Fit on idx_tr, predict on idx_te, using precomputed X
   fit_predict_idx <- function(idx_tr, idx_te) {
     n_te <- length(idx_te)
     preds <- rep(NA_real_, n_te)

     y_tr <- as.numeric(y_all[idx_tr])
     w_tr <- if (is.null(w_all)) NULL else as.numeric(w_all[idx_tr])

     # fallbacks: too few training rows
     if (length(y_tr) < 2L) {
       mu <- if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
       preds[] <- mu
       return(preds)
     }

     # fallback: constant outcome (treat NA as not helping)
     uy <- unique(y_tr[is.finite(y_tr)])
     if (length(uy) < 2L) {
       preds[] <- y_tr[which(is.finite(y_tr))[1L]]
       return(preds)
     }

     rf <- do.call(
       grf::regression_forest,
       c(list(
         X = getX(idx_tr),
         Y = y_tr,
         sample.weights = w_tr
       ), forest_opts)
     )

     preds[] <- as.numeric(predict(rf, getX(idx_te))$predictions)
     preds
   }

   # Build group index lists once
   idx_all <- i
   if (length(grp_cols) == 0L) {
     group_list <- list(idx_all)
   } else {
     group_dt <- DT[i, ..grp_cols]
     gid <- data.table::groupid(group_dt)
     group_list <- split(idx_all, gid)
   }

   # Main loop over groups
   for (g in seq_along(group_list)) {
     idx_g <- group_list[[g]]
     if (length(idx_g) == 0L) next

     idx1 <- idx_g[sample_all[idx_g] == 1L]
     idx2 <- idx_g[sample_all[idx_g] == 2L]

     # 1 -> 2
     if (length(idx2) > 0L) {
       p_12 <- fit_predict_idx(idx1, idx2)
       preds_buf[pos_map(idx2)] <- p_12
     }

     # 2 -> 1
     if (length(idx1) > 0L) {
       p_21 <- fit_predict_idx(idx2, idx1)
       preds_buf[pos_map(idx1)] <- p_21
     }

     if (verbose && (g %% 25L == 0L)) {
       message(sprintf("crossfit_hat_fast: processed %d/%d groups", g, length(group_list)))
     }
   }

   # single writeback and cleanup
   DT[i, (hat_col) := preds_buf]
   DT[, (rid_col) := NULL]

   invisible(DT)
 }

 crossfit_hat_slow <- function(DT,
                          i = NULL,
                          y_name,
                          x_names,
                          margins = NULL,
                          weight_name = NULL,
                          forest_opts = list(),
                          hat_suffix = ".hat") {
   stopifnot(data.table::is.data.table(DT))
   stopifnot(is.character(y_name), length(y_name) == 1L)
   stopifnot(is.character(x_names), length(x_names) >= 1L)
   stopifnot("sample" %chin% names(DT))
   stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

   if (is.null(i)) i <- DT[, .I]
   stopifnot(is.integer(i) || is.numeric(i))
   i <- as.integer(i)

   if (is.null(margins)) margins <- character()
   grp_cols <- as.character(margins)

   missing_cols <- setdiff(
     c("sample", y_name, x_names, grp_cols, if (!is.null(weight_name)) weight_name),
     names(DT)
   )
   if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

   hat_col <- paste0(y_name, hat_suffix)
   if (!(hat_col %chin% names(DT))) DT[, (hat_col) := NA_real_]

   # work on subset, but write back to DT using i
   subDT <- DT[i]

   fit_predict <- function(train_dt, test_dt) {
     y_tr <- as.numeric(train_dt[[y_name]])
     X_tr <- as.matrix(train_dt[, ..x_names])
     X_te <- as.matrix(test_dt[, ..x_names])
     w_tr <- if (is.null(weight_name)) NULL else as.numeric(train_dt[[weight_name]])

     preds <- rep(NA_real_, nrow(test_dt))

     # fallback: too few training rows
     if (length(y_tr) < 2L) {
       mu <- if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
       preds[] <- mu
       return(preds)
     }

     # fallback: constant outcome
     if (length(unique(y_tr)) < 2L) {
       preds[] <- y_tr[1L]
       return(preds)
     }

     rf <- do.call(
       grf::regression_forest,
       c(list(
         X = X_tr,
         Y = y_tr,
         sample.weights = w_tr
       ), forest_opts)
     )

     preds[] <- predict(rf, X_te)$predictions
     preds
   }

   # sdcols without NULLs (important inside data.table)
   sdcols <- c("sample", y_name, x_names, grp_cols)
   if (!is.null(weight_name)) sdcols <- c(sdcols, weight_name)
   sdcols <- unique(sdcols)

   if (length(grp_cols) == 0L) {
     dt1 <- subDT[sample == 1L]
     dt2 <- subDT[sample == 2L]

     # predictions in subset order
     preds_sub <- rep(NA_real_, nrow(subDT))
     preds_sub[subDT[["sample"]] == 2L] <- fit_predict(dt1, dt2)  # 1 -> 2
     preds_sub[subDT[["sample"]] == 1L] <- fit_predict(dt2, dt1)  # 2 -> 1

     DT[i, (hat_col) := preds_sub]
   } else {
     # compute predictions within each margins cell of the subset,
     # returning a vector aligned with subDT's row order
     preds_sub <- subDT[, {
       dt1 <- .SD[sample == 1L]
       dt2 <- .SD[sample == 2L]

       preds <- rep(NA_real_, .N)
       if (nrow(dt2) > 0L) preds[.SD[["sample"]] == 2L] <- fit_predict(dt1, dt2)
       if (nrow(dt1) > 0L) preds[.SD[["sample"]] == 1L] <- fit_predict(dt2, dt1)
       list(preds = preds)
     }, by = grp_cols, .SDcols = sdcols][["preds"]]

     DT[i, (hat_col) := preds_sub]
   }

   invisible(DT)
 }



 # ========= Helper: Estimate causal/regression/instrumental forests =========
 fit_models <- function(DT,
                             i = NULL,
                             forest_type = c("causal", "regression", "instrumental"),
                             y_name,
                             x_names,
                             margins = NULL,
                             w_name = NULL,
                             z_name = NULL,
                             folds = NULL,              # name of fold-id column (within-sample crossfit for pred only)
                             weight_name = NULL,
                             cluster_name = NULL,
                             forest_opts = NULL,
                             aipw_clip = 1e-3,
                             verbose = FALSE) {

   stopifnot(data.table::is.data.table(DT))
   forest_type <- match.arg(forest_type)

   stopifnot("sample" %chin% names(DT))
   stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

   if (!is.null(folds)) stopifnot(is.character(folds), length(folds) == 1L, folds %chin% names(DT))
   if (!is.null(aipw_clip)) {
     stopifnot(is.numeric(aipw_clip), length(aipw_clip) == 1L, is.finite(aipw_clip), aipw_clip > 0, aipw_clip < 1)
   }

   if (is.null(i)) i <- DT[, .I]
   i <- as.integer(i)
   n_i <- length(i)
   if (n_i == 0L) return(invisible(DT))

   if (is.null(margins)) margins <- character()
   grp_cols <- as.character(margins)

   if (is.null(forest_opts)) forest_opts <- list()

   # Output cols
   if (!("pred"   %chin% names(DT))) DT[, pred   := NA_real_]
   if (!("pred_o" %chin% names(DT))) DT[, pred_o := NA_real_]
   if (!("scores" %chin% names(DT))) DT[, scores := NA_real_]

   # Nuisance names (must exist for causal/IV)
   y_hat <- paste0(y_name, ".hat")
   w_hat <- if (!is.null(w_name)) paste0(w_name, ".hat") else NULL
   z_hat <- if (!is.null(z_name)) paste0(z_name, ".hat") else NULL

   if (forest_type %in% c("causal", "instrumental")) {
     stopifnot(!is.null(w_name), w_name %chin% names(DT))
     stopifnot(y_hat %chin% names(DT), w_hat %chin% names(DT))
   }
   if (forest_type == "instrumental") {
     stopifnot(!is.null(z_name), z_name %chin% names(DT))
     stopifnot(z_hat %chin% names(DT))
   }

   # ---- Precompute X for rows i ONCE ----
   X_all <- as.matrix(DT[i, ..x_names])

   # rowid within the i-subset
   rid_col <- ".__rid__"
   DT[i, (rid_col) := seq_len(n_i)]

   # helper: extract X rows for DT-row indices
   getX <- function(idx) X_all[DT[[rid_col]][idx], , drop = FALSE]

   # Pull commonly used vectors once
   sample_all <- DT[["sample"]]
   folds_all  <- if (!is.null(folds)) DT[[folds]] else NULL
   wgt_all    <- if (!is.null(weight_name)) DT[[weight_name]] else NULL
   cl_all     <- if (!is.null(cluster_name)) DT[[cluster_name]] else NULL

   y_all <- DT[[y_name]]
   if (forest_type %in% c("causal", "instrumental")) {
     w_all    <- DT[[w_name]]
     yhat_all <- DT[[y_hat]]
     what_all <- DT[[w_hat]]
   }
   if (forest_type == "instrumental") {
     z_all    <- DT[[z_name]]
     zhat_all <- DT[[z_hat]]
   }

   # ---- Forest builder taking DT-row indices ----
   build_forest_idx <- function(type, idx) {
     X <- getX(idx)
     Y <- as.numeric(y_all[idx])

     sw <- if (is.null(wgt_all)) NULL else as.numeric(wgt_all[idx])
     cl <- if (is.null(cl_all)) NULL else cl_all[idx]

     if (type == "regression") {
       do.call(grf::regression_forest,
               c(list(X = X, Y = Y, sample.weights = sw, clusters = cl), forest_opts))
     } else if (type == "causal") {
       W <- as.numeric(w_all[idx])
       do.call(grf::causal_forest,
               c(list(X = X, Y = Y, W = W,
                      Y.hat = as.numeric(yhat_all[idx]),
                      W.hat = as.numeric(what_all[idx]),
                      sample.weights = sw,
                      clusters = cl), forest_opts))
     } else { # instrumental
       W <- as.numeric(w_all[idx])
       Z <- as.numeric(z_all[idx])
       do.call(grf::instrumental_forest,
               c(list(X = X, Y = Y, W = W, Z = Z,
                      Y.hat = as.numeric(yhat_all[idx]),
                      W.hat = as.numeric(what_all[idx]),
                      Z.hat = as.numeric(zhat_all[idx]),
                      sample.weights = sw,
                      clusters = cl), forest_opts))
     }
   }

   # ---- Within-sample predictions (pred) ----
   within_pred_idx <- function(idx_s) {
     Xs <- getX(idx_s)

     if (is.null(folds_all)) {
       fit <- build_forest_idx(forest_type, idx_s)
       return(as.numeric(predict(fit, Xs)$predictions))
     }

     f <- folds_all[idx_s]
     u <- unique(f)
     p <- rep(NA_real_, length(idx_s))

     for (k in u) {
       te_pos <- which(f == k)
       tr_pos <- which(f != k)
       fit_k <- build_forest_idx(forest_type, idx_s[tr_pos])
       p[te_pos] <- as.numeric(predict(fit_k, Xs[te_pos, , drop = FALSE])$predictions)
     }
     p
   }

   # ---- AIPW scores using OUT-OF-SAMPLE tau (pred_o) ----
   compute_aipw_vec <- function(idx_s, tau_vec) {
     n <- length(idx_s)
     if (forest_type == "regression") return(as.numeric(y_all[idx_s]))
     if (forest_type != "causal")     return(rep(NA_real_, n))

     Y <- as.numeric(y_all[idx_s])
     m <- as.numeric(yhat_all[idx_s])
     e <- as.numeric(what_all[idx_s])
     W <- as.numeric(w_all[idx_s])

     W <- as.numeric(W > 0.5)

     # if e not in [0,1], assume logit index and convert
     if (any(e < 0 | e > 1, na.rm = TRUE)) e <- plogis(e)

     if (!is.null(aipw_clip)) e <- pmin(pmax(e, aipw_clip), 1 - aipw_clip)

     tau <- as.numeric(tau_vec)
     stopifnot(length(tau) == n)

     m1 <- m + (1 - e) * tau
     m0 <- m - e * tau

     tau + (W / e) * (Y - m1) - ((1 - W) / (1 - e)) * (Y - m0)
   }

   # ---- Build group index lists once ----
   idx_all <- i
   if (length(grp_cols) == 0L) {
     group_list <- list(idx_all)
   } else {
     group_dt <- DT[i, ..grp_cols]
     gid <- data.table::groupid(group_dt)
     group_list <- split(idx_all, gid)
   }

   # ---- Output buffers ----
   pred_buf   <- rep(NA_real_, n_i)
   pred_o_buf <- rep(NA_real_, n_i)
   score_buf  <- rep(NA_real_, n_i)

   # map DT-row indices -> positions in buffers
   # (safer than integer(max(i)) and fast enough)
   pos_map <- function(idx) match(idx, i)

   # ---- Main loop ----
   for (g in seq_along(group_list)) {
     idx_g <- group_list[[g]]
     if (length(idx_g) == 0L) next

     idx1 <- idx_g[sample_all[idx_g] == 1L]
     idx2 <- idx_g[sample_all[idx_g] == 2L]

     p1 <- if (length(idx1)) within_pred_idx(idx1) else numeric()
     p2 <- if (length(idx2)) within_pred_idx(idx2) else numeric()

     fit1 <- if (length(idx1)) build_forest_idx(forest_type, idx1) else NULL
     fit2 <- if (length(idx2)) build_forest_idx(forest_type, idx2) else NULL

     p1_o <- if (!is.null(fit2) && length(idx1)) as.numeric(predict(fit2, getX(idx1))$predictions) else numeric()
     p2_o <- if (!is.null(fit1) && length(idx2)) as.numeric(predict(fit1, getX(idx2))$predictions) else numeric()

     sc1 <- if (length(idx1)) compute_aipw_vec(idx1, p1_o) else numeric()
     sc2 <- if (length(idx2)) compute_aipw_vec(idx2, p2_o) else numeric()

     if (length(idx1)) {
       pos <- pos_map(idx1)
       pred_buf[pos]   <- p1
       pred_o_buf[pos] <- p1_o
       score_buf[pos]  <- sc1
     }
     if (length(idx2)) {
       pos <- pos_map(idx2)
       pred_buf[pos]   <- p2
       pred_o_buf[pos] <- p2_o
       score_buf[pos]  <- sc2
     }

     if (verbose && (g %% 25L == 0L)) {
       message(sprintf("fit_models_fast: processed %d/%d groups", g, length(group_list)))
     }
   }

   # ---- Single writeback ----
   DT[i, `:=`(pred = pred_buf, pred_o = pred_o_buf, scores = score_buf)]

   # cleanup (IMPORTANT: no i when deleting columns)
   DT[, (rid_col) := NULL]

   invisible(DT)
 }

 fit_models_slow <- function(DT,
                        i = NULL,
                        forest_type = c("causal", "regression", "instrumental"),
                        y_name,
                        x_names,
                        margins = NULL,
                        w_name = NULL,
                        z_name = NULL,
                        folds = NULL,              # name of fold-id column (within-sample crossfit for pred only)
                        weight_name = NULL,
                        cluster_name = NULL,
                        forest_opts = NULL,
                        aipw_clip = 1e-3) {

   stopifnot(data.table::is.data.table(DT))
   forest_type <- match.arg(forest_type)

   stopifnot("sample" %chin% names(DT))
   stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

   if (!is.null(folds)) {
     stopifnot(is.character(folds), length(folds) == 1L, folds %chin% names(DT))
   }

   if (is.null(i)) i <- DT[, .I]
   i <- as.integer(i)

   if (is.null(margins)) margins <- character()
   grp_cols <- as.character(margins)

   if (is.null(forest_opts)) forest_opts <- list()

   # Output columns (only filled on rows i)
   if (!("pred"   %chin% names(DT))) DT[, pred   := NA_real_]
   if (!("pred_o" %chin% names(DT))) DT[, pred_o := NA_real_]
   if (!("scores" %chin% names(DT))) DT[, scores := NA_real_]

   # Nuisances (must already exist for causal/IV)
   y_hat <- paste0(y_name, ".hat")
   w_hat <- if (!is.null(w_name)) paste0(w_name, ".hat") else NULL
   z_hat <- if (!is.null(z_name)) paste0(z_name, ".hat") else NULL

   if (forest_type %in% c("causal", "instrumental")) {
     stopifnot(!is.null(w_name), w_name %chin% names(DT))
     stopifnot(y_hat %chin% names(DT), w_hat %chin% names(DT))
   }
   if (forest_type == "instrumental") {
     stopifnot(!is.null(z_name), z_name %chin% names(DT))
     stopifnot(z_hat %chin% names(DT))
   }

   build_forest <- function(type, train_dt) {
     X <- as.matrix(train_dt[, ..x_names])
     Y <- as.numeric(train_dt[[y_name]])

     sw <- if (is.null(weight_name)) NULL else as.numeric(train_dt[[weight_name]])
     cl <- if (is.null(cluster_name)) NULL else train_dt[[cluster_name]]

     if (type == "regression") {
       do.call(
         grf::regression_forest,
         c(list(X = X, Y = Y, sample.weights = sw, clusters = cl), forest_opts)
       )
     } else if (type == "causal") {
       W <- as.numeric(train_dt[[w_name]])
       do.call(
         grf::causal_forest,
         c(list(
           X = X, Y = Y, W = W,
           Y.hat = as.numeric(train_dt[[y_hat]]),
           W.hat = as.numeric(train_dt[[w_hat]]),
           sample.weights = sw,
           clusters = cl
         ), forest_opts)
       )
     } else { # instrumental
       W <- as.numeric(train_dt[[w_name]])
       Z <- as.numeric(train_dt[[z_name]])
       do.call(
         grf::instrumental_forest,
         c(list(
           X = X, Y = Y, W = W, Z = Z,
           Y.hat = as.numeric(train_dt[[y_hat]]),
           W.hat = as.numeric(train_dt[[w_hat]]),
           Z.hat = as.numeric(train_dt[[z_hat]]),
           sample.weights = sw,
           clusters = cl
         ), forest_opts)
       )
     }
   }

   # Work on subset
   subDT <- DT[i]

   # .SDcols without NULLs
   sdcols <- c("sample", y_name, x_names, grp_cols)
   if (!is.null(folds)) sdcols <- c(sdcols, folds)
   if (forest_type %in% c("causal", "instrumental")) sdcols <- c(sdcols, w_name, y_hat, w_hat)
   if (forest_type == "instrumental") sdcols <- c(sdcols, z_name, z_hat)
   if (!is.null(weight_name)) sdcols <- c(sdcols, weight_name)
   if (!is.null(cluster_name)) sdcols <- c(sdcols, cluster_name)
   sdcols <- unique(sdcols)

   # Within-sample tau prediction: OOB-ish if folds is NULL (via predict on training X),
   # or explicit K-fold within sample if folds is provided.
   within_pred <- function(dt_s) {
     X <- as.matrix(dt_s[, ..x_names])

     if (is.null(folds)) {
       fit <- build_forest(forest_type, dt_s)
       return(as.numeric(predict(fit, X)$predictions))
     }

     f <- dt_s[[folds]]
     p <- rep(NA_real_, nrow(dt_s))
     for (k in unique(f)) {
       te <- (f == k)
       tr <- !te
       fit_k <- build_forest(forest_type, dt_s[tr])
       p[te] <- as.numeric(predict(fit_k, X[te, , drop = FALSE])$predictions)
     }
     p
   }

   # Canonical AIPW score for binary treatment:
   # psi = tau + W/e*(Y-m1) - (1-W)/(1-e)*(Y-m0)
   # where m1 = m + (1-e)*tau, m0 = m - e*tau, and m = E[Y|X], e = P(W=1|X).
   compute_aipw <- function(dt_s, tau_vec) {
     if (forest_type == "regression") {
       return(as.numeric(dt_s[[y_name]]))
     }
     if (forest_type != "causal") {
       return(rep(NA_real_, nrow(dt_s)))
     }

     Y <- as.numeric(dt_s[[y_name]])
     m <- as.numeric(dt_s[[y_hat]])
     e <- as.numeric(dt_s[[w_hat]])
     W <- as.numeric(dt_s[[w_name]])

     # clip propensity
     #if (is.null(aipw_clip)==FALSE) e <- pmin(pmax(e, aipw_clip), 1 - aipw_clip)

     tau <- as.numeric(tau_vec)
     m1 <- m + (1 - e) * tau
     m0 <- m - e * tau

     tau + (W / e) * (Y - m1) - ((1 - W) / (1 - e)) * (Y - m0)
   }

   compute_group <- function(dtg) {
     dt1 <- dtg[sample == 1L]
     dt2 <- dtg[sample == 2L]

     # pred: within-sample (fold crossfit if folds provided)
     p1 <- if (nrow(dt1)) within_pred(dt1) else numeric()
     p2 <- if (nrow(dt2)) within_pred(dt2) else numeric()

     # pred_o: across-sample (fit on full sample 1/2 within this margins cell)
     fit1_full <- if (nrow(dt1)) build_forest(forest_type, dt1) else NULL
     fit2_full <- if (nrow(dt2)) build_forest(forest_type, dt2) else NULL

     X1 <- if (nrow(dt1)) as.matrix(dt1[, ..x_names]) else NULL
     X2 <- if (nrow(dt2)) as.matrix(dt2[, ..x_names]) else NULL

     p1_o <- if (!is.null(fit2_full)) as.numeric(predict(fit2_full, X1)$predictions) else numeric() # 2 -> 1
     p2_o <- if (!is.null(fit1_full)) as.numeric(predict(fit1_full, X2)$predictions) else numeric() # 1 -> 2

     # scores: AIPW, computed using OUT-OF-SAMPLE tau (pred_o)
     sc1 <- if (nrow(dt1)) compute_aipw(dt1, p1_o) else numeric()
     sc2 <- if (nrow(dt2)) compute_aipw(dt2, p2_o) else numeric()

     # Align to dtg row order
     out_pred   <- rep(NA_real_, nrow(dtg))
     out_pred_o <- rep(NA_real_, nrow(dtg))
     out_scores <- rep(NA_real_, nrow(dtg))

     s <- dtg[["sample"]]
     if (nrow(dt1)) {
       out_pred[s == 1L]   <- p1
       out_pred_o[s == 1L] <- p1_o
       out_scores[s == 1L] <- sc1
     }
     if (nrow(dt2)) {
       out_pred[s == 2L]   <- p2
       out_pred_o[s == 2L] <- p2_o
       out_scores[s == 2L] <- sc2
     }

     list(pred = out_pred, pred_o = out_pred_o, scores = out_scores)
   }

   if (length(grp_cols) == 0L) {
     res <- compute_group(subDT)
     DT[i, `:=`(pred = res$pred, pred_o = res$pred_o, scores = res$scores)]
   } else {
     res <- subDT[, compute_group(.SD), by = grp_cols, .SDcols = sdcols]
     DT[i, `:=`(
       pred   = res[["pred"]],
       pred_o = res[["pred_o"]],
       scores = res[["scores"]]
     )]
   }

   invisible(DT)
 }



 ####### FIND OPTIMAL CUTOFF AND TEST #############
 forest_test <- function(
    data,
    cluster = NULL,          # string or NULL
    weight  = NULL,          # NULL or string (weight column)
    sample  = "sample",
    pred    = "pred",
    pred_o  = "pred_o",
    scores  = "scores",      # AIPW scores used for BOTH training + testing
    minsize = 50L,
    pool = TRUE,
    gridpoints = NULL        # NULL => all cutoffs; else integer, e.g. 100 (weighted percentiles of pred)
 ) {
   stopifnot(data.table::is.data.table(data))

   # column name strings
   sample_col <- sample
   pred_col   <- pred
   pred_o_col <- pred_o
   scores_col <- scores
   weight_col <- weight

   stopifnot(sample_col %chin% names(data))
   stopifnot(pred_col   %chin% names(data))
   stopifnot(pred_o_col %chin% names(data))
   stopifnot(scores_col %chin% names(data))
   stopifnot(is.numeric(minsize), length(minsize) == 1L, minsize >= 1L)

   if (!is.null(weight_col)) stopifnot(weight_col %chin% names(data))
   if (!is.null(cluster))    stopifnot(cluster %chin% names(data))

   if (!is.null(gridpoints)) {
     stopifnot(is.numeric(gridpoints), length(gridpoints) == 1L, is.finite(gridpoints), gridpoints >= 2)
     gridpoints <- as.integer(gridpoints)
   }
   grid_on <- !is.null(gridpoints)

   # ---- Work on a small DT with only needed columns (major speed + avoids modifying 'data') ----
   # Use integer cluster id for faster duplicated/uniqueN.
   dt <- data.table::data.table(
     sample = as.integer(data[[sample_col]]),
     pred   = as.numeric(data[[pred_col]]),
     pred_o = as.numeric(data[[pred_o_col]]),
     score  = as.numeric(data[[scores_col]])
   )

   if (!is.null(weight_col)) {
     dt[, w := as.numeric(data[[weight_col]])]
   } else {
     dt[, w := 1.0]
   }

   if (!is.null(cluster)) {
     # integerize cluster ids (fast)
     dt[, cl := as.integer(factor(data[[cluster]], exclude = NULL))]
   } else {
     # iid fallback: each row is its own cluster
     dt[, cl := .I]
   }

   stopifnot(all(dt$sample %in% c(1L, 2L)))

   # Helpers -------------------------------------------------------------

   # CRV1 SE for weighted mean using cluster sums:
   # theta = sum_g U_g / sum_g W_g, where U_g = sum_i w_i * score_i, W_g = sum_i w_i
   # se = sqrt(G/(G-1) * sum_g (U_g - theta W_g)^2) / |sum_g W_g|
   crv1_mean <- function(dsub) {
     # cluster sums
     gb <- dsub[, .(U = sum(w * score), W = sum(w)), by = cl]
     G <- nrow(gb)
     U <- sum(gb$U)
     W <- sum(gb$W)
     theta <- U / W

     ug <- gb$U - theta * gb$W
     se <- sqrt((G / pmax.int(G - 1L, 1L)) * sum(ug^2)) / abs(W)
     if (G < 2L || !is.finite(se)) se <- NA_real_

     list(coef = theta, se = se, t = theta / se, G = G, N = nrow(dsub))
   }

   # Find (weighted) percentile candidate indices in current order within each sample.
   # Returns a logical vector cand (same length as dsub) indicating candidate rows.
   make_grid_candidates <- function(dsub, gridpoints) {
     n <- nrow(dsub)
     if (n <= gridpoints) return(rep(TRUE, n))

     cw <- cumsum(dsub$w)
     totw <- cw[n]
     if (!is.finite(totw) || totw <= 0) return(rep(TRUE, n))

     probs <- (seq_len(gridpoints)) / (gridpoints + 1)  # avoid 0/1 endpoints
     targets <- probs * totw
     idx <- vapply(targets, function(tt) which(cw >= tt)[1L], integer(1))
     idx <- unique(pmin(pmax(idx, 1L), n))

     cand <- rep(FALSE, n)
     cand[idx] <- TRUE
     cand
   }

   # --------------------------------------------------------------------
   # TRAINING objective computed as in your current function:
   # order by (sample, pred), compute running mean + CRV1 SE via cumulative cluster contributions.
   # --------------------------------------------------------------------

   # Sort once by sample,pred
   data.table::setorder(dt, sample, pred)

   # Running counts (for reporting)
   dt[, N := seq_len(.N), by = sample]

   # a = w*score, b = w
   dt[, `:=`(a = w * score, b = w)]

   # Within-cluster running totals, per sample (in current pred-order)
   dt[, `:=`(WgY = cumsum(a), Wg = cumsum(b)), by = .(sample, cl)]

   # Global running totals per sample
   dt[, `:=`(SW = cumsum(b), SWY = cumsum(a)), by = sample]
   dt[, m := SWY / SW, by = sample]

   # Per-row deltas for cluster-robust variance recursion
   dt[, `:=`(
     dTA2 =  WgY^2 - (WgY - a)^2,
     dTB2 =  Wg^2  - (Wg  - b)^2,
     dTAB =  WgY * Wg - (WgY - a) * (Wg - b)
   )]

   # Cumulative aggregates per sample
   dt[, `:=`(TA2 = cumsum(dTA2), TB2 = cumsum(dTB2), TAB = cumsum(dTAB)), by = sample]

   # Number of clusters seen so far (in this pred-order)
   dt[, G := cumsum(!duplicated(cl)), by = sample]

   # CRV1 SE of weighted mean
   dt[, sumS2 := (TA2 - 2*m*TAB + (m^2)*TB2) / (SW^2)]
   dt[, se := sqrt((G / pmax.int(G - 1L, 1L)) * sumS2)]
   dt[G < 2L, se := NA_real_]
   dt[, t := m / se]

   # --------------------------------------------------------------------
   # Grid candidates (optional)
   # --------------------------------------------------------------------
   if (grid_on) {
     dt[, cand := make_grid_candidates(.SD, gridpoints), by = sample]
   } else {
     dt[, cand := TRUE]
   }

   # --------------------------------------------------------------------
   # Dual constraint tau_vec:
   #   tau_tr: smallest pred where running clusters >= minsize (train-side)
   #   tau_est: smallest pred_o where running clusters >= minsize (est-side)
   # --------------------------------------------------------------------

   # train-side: dt already ordered by pred
   tau_tr <- dt[G >= minsize, .(tau_tr = min(pred, na.rm = TRUE)), by = sample]

   # est-side: order by pred_o within sample and count clusters
   dt_est <- dt[, .(sample, pred_o, cl, w)]  # small view
   data.table::setorder(dt_est, sample, pred_o)
   dt_est[, G_est := cumsum(!duplicated(cl)), by = sample]
   tau_est <- dt_est[G_est >= minsize, .(tau_est = min(pred_o, na.rm = TRUE)), by = sample]

   # If constraints are missing in a sample, fail early (better than silent NA bugs)
   if (nrow(tau_tr) < 2L || nrow(tau_est) < 2L) {
     stop("minsize is too large: cannot reach minsize clusters in train or est order for both samples.")
   }

   tau_vec <- c(
     max(as.numeric(c(tau_tr[sample == 1L, tau_tr], tau_est[sample == 2L, tau_est]))),
     max(as.numeric(c(tau_tr[sample == 2L, tau_tr], tau_est[sample == 1L, tau_est])))
   )

   # --------------------------------------------------------------------
   # Choose best cutoff in each sample among eligible rows:
   # eligible: pred >= tau_constraint, G >= minsize, cand==TRUE
   # --------------------------------------------------------------------

   # eligibility constraint tau (per sample) for searching within that sample
   dt[, tau_constraint := data.table::fifelse(sample == 1L, tau_vec[1], tau_vec[2])]

   elig <- (dt$pred >= dt$tau_constraint) & (dt$G >= minsize) & (dt$cand)

   res <- dt[elig, .SD[which.min(t)], by = sample, .SDcols = c("G", "N", "m", "se", "t", "pred")]
   res[, train := TRUE]
   data.table::setcolorder(res, c("train", "sample"))

   # robust extraction with fallback to constraint if one side is missing
   cut1 <- res[sample == 1L, pred]
   cut2 <- res[sample == 2L, pred]
   if (length(cut1) == 0L || !is.finite(cut1)) cut1 <- tau_vec[1]
   if (length(cut2) == 0L || !is.finite(cut2)) cut2 <- tau_vec[2]

   # swapped testing cutoffs: sample1 uses cut2, sample2 uses cut1
   dt[, tau_test := data.table::fifelse(sample == 1L, cut2, cut1)]

   # --------------------------------------------------------------------
   # TESTING set and test statistic (FAST; no fixest)
   # --------------------------------------------------------------------
   dtest <- dt[pred_o <= tau_test]

   if (!pool) {
     # per sample
     out1 <- crv1_mean(dtest[sample == 1L])
     out2 <- crv1_mean(dtest[sample == 2L])

     test_out <- data.table::data.table(
       train = FALSE,
       sample = c(1L, 2L),
       G = c(out1$G, out2$G),
       N = c(out1$N, out2$N),
       coef = c(out1$coef, out2$coef),
       stderr = c(out1$se, out2$se),
       t = c(out1$t, out2$t),
       `tau cutoff` = NA_real_
     )
   } else {
     # pooled
     out <- crv1_mean(dtest)

     test_out <- data.table::data.table(
       train = FALSE,
       sample = NA_integer_,
       G = out$G,
       N = out$N,
       coef = out$coef,
       stderr = out$se,
       t = out$t,
       `tau cutoff` = NA_real_
     )
   }

   # --------------------------------------------------------------------
   # OUTPUT ASSEMBLY (match your schema)
   # --------------------------------------------------------------------
   train_out <- data.table::data.table(
     train = TRUE,
     sample = res$sample,
     G = res$G,
     N = res$N,
     coef = res$m,
     stderr = res$se,
     t = res$t,
     `tau cutoff` = res$pred
   )

   res_out <- data.table::rbindlist(list(train_out, test_out), use.names = TRUE, fill = TRUE)

   # optional p-values (one-sided for negative effects): H0: >=0 vs H1:<0
   res_out[, p.raw := stats::pnorm(t)]

   res_out
 }

 forest_test_slow <- function(
    data,
    cluster = NULL,          # string or NULL
    weight  = NULL,          # NULL or string (weight column)
    sample  = "sample",
    pred    = "pred",
    pred_o  = "pred_o",
    scores  = "scores",      # AIPW scores used for BOTH training + testing
    minsize = 50L,
    pool = TRUE,
    gridpoints = NULL        # NULL => all cutoffs; else integer, e.g. 100 (weighted percentiles of pred)
 ) {

   stopifnot(data.table::is.data.table(data))

   # --- freeze column-name arguments to avoid NSE collisions ---
   sample_col <- sample
   pred_col   <- pred
   pred_o_col <- pred_o
   scores_col <- scores
   weight_col <- weight

   stopifnot(is.character(sample_col), length(sample_col) == 1L, sample_col %chin% names(data))
   stopifnot(is.character(pred_col),   length(pred_col)   == 1L, pred_col   %chin% names(data))
   stopifnot(is.character(pred_o_col), length(pred_o_col) == 1L, pred_o_col %chin% names(data))
   stopifnot(is.character(scores_col), length(scores_col) == 1L, scores_col %chin% names(data))
   stopifnot(is.numeric(minsize) && minsize >= 1)

   if (!is.null(weight_col)) {
     stopifnot(is.character(weight_col), length(weight_col) == 1L, weight_col %chin% names(data))
   }

   if (!is.null(gridpoints)) {
     stopifnot(is.numeric(gridpoints), length(gridpoints) == 1L, is.finite(gridpoints), gridpoints >= 2)
     gridpoints <- as.integer(gridpoints)
   }
   grid_on <- !is.null(gridpoints)

   # cluster handling: if NULL, create per-row clusters (iid-like fallback)
   tmp_cluster <- FALSE
   if (is.null(cluster)) {
     cluster_col <- ".__cl__"
     tmp_cluster <- TRUE
     data[, (cluster_col) := .I]
   } else {
     stopifnot(is.character(cluster), length(cluster) == 1L, cluster %chin% names(data))
     cluster_col <- cluster
   }

   # Unique internal columns
   tau_col   <- ".__tau__"
   cand_col  <- ".__cand__"
   wtau_col  <- ".__wtau__"
   cwtau_col <- ".__cwtau__"

   # ---------------------------
   # TRAINING: choose cutoff to minimize t-stat of mean(AIPW score)
   # ---------------------------

   # Order by training predictor within sample
   data.table::setorderv(data, cols = c(sample_col, pred_col))

   # Running index per sample (for reporting only)
   data[, N := seq_len(.N), by = sample_col]

   # Weighted contributions for mean(scores)
   if (!is.null(weight_col)) {
     data[, `:=`(
       a = get(weight_col) * get(scores_col),
       b = get(weight_col)
     )]
   } else {
     data[, `:=`(
       a = get(scores_col),
       b = 1.0
     )]
   }

   # Within-cluster running totals (per current order)
   data[, `:=`(WgY = cumsum(a), Wg = cumsum(b)), by = c(sample_col, cluster_col)]

   # Global running totals and running mean by sample
   data[, `:=`(SW = cumsum(b), SWY = cumsum(a)), by = sample_col]
   data[, m := SWY / SW, by = sample_col]

   # Per-row deltas (avoid shifts)
   data[, `:=`(
     dTA2 =  WgY^2 - (WgY - a)^2,
     dTB2 =  Wg^2  - (Wg  - b)^2,
     dTAB =  WgY*Wg - (WgY - a)*(Wg - b)
   )]

   # Cumulative aggregates across rows (by sample)
   data[, `:=`(TA2 = cumsum(dTA2), TB2 = cumsum(dTB2), TAB = cumsum(dTAB)), by = sample_col]

   # Number of unique clusters seen so far (by sample, current order)
   data[, G := cumsum(!duplicated(get(cluster_col))), by = sample_col]

   # CRV1 (small-sample adj.) SE of weighted mean
   data[, sumS2 := (TA2 - 2*m*TAB + (m^2)*TB2) / (SW^2)]
   data[, se := sqrt((G / pmax.int(G - 1L, 1L)) * sumS2)]
   data[G < 2, se := NA_real_]

   # t-stat for H0: mean = 0
   data[, t := m / se]

   # ---------------------------
   # Optional grid restriction: only consider cutoffs at (weighted) percentiles of pred
   # ---------------------------
   if (grid_on) {
     if (!is.null(weight_col)) {
       data[, (wtau_col) := get(weight_col)]
     } else {
       data[, (wtau_col) := 1.0]
     }

     probs <- (seq_len(gridpoints)) / (gridpoints + 1)  # avoid 0/1 endpoints

     # NOTE: dynamic assignment must use c(colnames) := list(...)
     data[, c(cwtau_col, cand_col) := {
       cw <- cumsum(get(wtau_col))
       n  <- .N
       cand <- rep(FALSE, n)

       if (n <= gridpoints) {
         cand[] <- TRUE
       } else {
         totw <- cw[n]
         if (!is.finite(totw) || totw <= 0) {
           cand[] <- TRUE
         } else {
           targets <- probs * totw
           idx <- vapply(targets, function(tt) which(cw >= tt)[1L], integer(1))
           idx <- unique(pmin(pmax(idx, 1L), n))
           cand[idx] <- TRUE
         }
       }
       list(cw, cand)
     }, by = sample_col]
   }

   # ---------------------------
   # Dual-constraint tau ensuring >= minsize clusters in BOTH samples
   # ---------------------------

   # Train-side cutoff by sample (ordered by pred already)
   tau_tr <- data[G >= minsize, .(tau_tr = min(get(pred_col), na.rm = TRUE)), by = sample_col]

   # Estimation-side cutoff by sample (order by pred_o, count clusters)
   data.table::setorderv(data, cols = c(sample_col, pred_o_col))
   data[, G_est_order := cumsum(!duplicated(get(cluster_col))), by = sample_col]
   tau_est <- data[G_est_order >= minsize, .(tau_est = min(get(pred_o_col), na.rm = TRUE)), by = sample_col]

   # Combine into cross-sample constraint (assumes exactly two samples coded 1/2)
   tau_vec <- c(
     max(as.numeric(c(tau_tr[1, tau_tr], tau_est[2, tau_est]))),
     max(as.numeric(c(tau_tr[2, tau_tr], tau_est[1, tau_est])))
   )

   # Create tau column NOW (before any get(tau_col) is called)
   data[, (tau_col) := fifelse(get(sample_col) == 1L, tau_vec[1], tau_vec[2])]

   # Back to training order for picking best cutoff
   data.table::setorderv(data, cols = c(sample_col, pred_col))

   # Eligible rows for optimization (and optional grid restriction)
   elig <- (data[[pred_col]] >= data[[tau_col]]) & (data[["G"]] >= minsize)
   if (grid_on) elig <- elig & data[[cand_col]]

   # Train-sample optimal cutoff: choose pred that minimizes t among eligible rows
   res <- data[elig,
               .SD[which.min(t)],
               by = sample_col,
               .SDcols = c("G", "N", "m", "se", "t", pred_col)]
   res[, train := TRUE]
   data.table::setcolorder(res, c("train", sample_col))

   # Replace tau by the optimal cutoff from the *other* sample (assumes two samples 1/2)
   cut1 <- res[get(sample_col) == 1L, get(pred_col)]
   cut2 <- res[get(sample_col) == 2L, get(pred_col)]
   data[, (tau_col) := data.table::fifelse(get(sample_col) == 1L, cut2, cut1)]

   # ---------------------------
   # TESTING: pred_o <= tau (using SAME AIPW scores column)
   # ---------------------------

   dtest <- data[get(pred_o_col) <= get(tau_col)]

   # Counts on testing set
   byv <- if (pool) NULL else sample_col
   GN <- dtest[, .(G.est = data.table::uniqueN(get(cluster_col)), N.est = .N), by = byv]

   cl <- if (!is.null(cluster)) as.formula(paste0("~", cluster)) else NULL
   wg <- if (!is.null(weight_col)) as.formula(paste0("~", weight_col)) else NULL

   if (!pool) {
     fml <- stats::as.formula(paste0(scores_col, " ~ i(", sample_col, ") - 1"))
     fit <- fixest::feols(fml, data = dtest, vcov = cl, weights = wg)

     test <- data.table::as.data.table(cbind(GN, fit$coeftable[, 1:3]))
     test[, train := FALSE]
     test[, tau_cutoff := NA_real_]
     data.table::setcolorder(test, c("train", sample_col))
   } else {
     fml <- stats::as.formula(paste0(scores_col, " ~ 1"))
     fit <- fixest::feols(fml, data = dtest, vcov = cl, weights = wg)

     test <- data.table::data.table(
       train = FALSE,
       G = GN$G.est,
       N = GN$N.est,
       coef = fit$coeftable[1, 1],
       stderr = fit$coeftable[1, 2],
       t = fit$coeftable[1, 3],
       `tau cutoff` = NA_real_
     )
     test[, (sample_col) := NA_integer_]
     data.table::setcolorder(test, c("train", sample_col))
   }

   # ---------------------------
   # OUTPUT ASSEMBLY
   # ---------------------------

   # Training rows
   train_out <- res[, .(
     train = TRUE,
     G = G,
     N = N,
     coef = m,
     stderr = se,
     t = t,
     `tau cutoff` = get(pred_col)
   )]
   train_out[, (sample_col) := res[[sample_col]] ]
   data.table::setcolorder(train_out, c("train", sample_col))

   # Testing rows into same schema
   if (!pool) {
     test_out <- test[, .(
       train = FALSE,
       G = G.est,
       N = N.est,
       coef = Estimate,
       stderr = `Std. Error`,
       t = `t value`,
       `tau cutoff` = NA_real_
     )]
     test_out[, (sample_col) := test[[sample_col]] ]
     data.table::setcolorder(test_out, c("train", sample_col))
   } else {
     test_out <- test[, .(
       train = FALSE,
       G = G,
       N = N,
       coef = coef,
       stderr = stderr,
       t = t,
       `tau cutoff` = `tau cutoff`
     )]
     test_out[, (sample_col) := test[[sample_col]] ]
     data.table::setcolorder(test_out, c("train", sample_col))
   }

   res_out <- data.table::rbindlist(list(train_out, test_out), use.names = TRUE, fill = TRUE)

   # ---------------------------
   # CLEANUP
   # ---------------------------

   drop_cols <- c(
     "N","a","b","WgY","Wg","dTA2","dTB2","dTAB","TA2","TB2","TAB",
     "sumS2","se","G","t","m","SW","SWY","G_est_order", tau_col
   )
   if (grid_on) drop_cols <- c(drop_cols, cand_col, wtau_col, cwtau_col)

   drop_cols <- drop_cols[drop_cols %chin% names(data)]
   if (length(drop_cols)) data[, (drop_cols) := NULL]
   if (tmp_cluster) data[, (cluster_col) := NULL]

   res_out
 }

