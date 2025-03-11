  #' @param infile Path to the input file
  #' @return A matrix of the infile
  #' @export

  montest=function(data,D,Z,X=NULL,Y=NULL,treefraction=0.5,control="forest",test=NULL,
                   normalize=TRUE,stack="all",treetype="forest",se=FALSE,saveforest=F,
                   gridtypeY="equidistant",gridtypeD="equidistant",gridtypeZ="equidistant",
                   Ysubsets = 4, Dsubsets = 4,Zsubsets=4,Y.res=TRUE,
                   weight=NULL,cluster=NULL,num.trees=2000,seed=10101,minsize=50,
                   maxrankcp=5,prune=TRUE,cp=0,alpha=0.05,preselect="nonpositive", ##CART options
                   Zparameters=NULL,Yparameters=NULL,Qparameters=NULL,Dparameters=NULL,Cparameters=NULL,
                   tune.Qparameters="none",tune.Zparameters="none",tune.Cparameters="none",tune.Yparameters="none",tune.Dparameters="none",
                   tune.num.trees=200,tune.num.reps=50,tune.num.draws=1000,tunetype="one" ##tuning options
                   ){

    time=proc.time()
    set.seed(seed)

    ##check input
    if (is.null(test)==TRUE) {
      if (is.null(Y)==TRUE) {
        test="simple"
      } else test="all"
    }

    test=match.arg(test,c("simple","BP","K","MW","AHS","all"),several.ok=TRUE)
    if ("all" %in% test) {
      test=c("simple","BP","MW","AHS")
    }

    if (is.null(Y)==TRUE) {
      if (sum(!test %in% "simple")>0) {
      stop("Other tests than simple may not be used when Y is not specified. Specify the Y argument or use test=simple")
      }
    }

    if ((Ysubsets<=1)|(Dsubsets<=1)|(Zsubsets<=1)) stop("Ysubsets, Dsubsets and Zsubsets must be integers larger than 1")

    if (is.null(stack)==FALSE) stack=match.arg(stack,c("all","Z","D","Y","Ybin","above","condition","none"),several.ok=TRUE)
    if (("all" %in% stack)&("none" %in% stack)) stop("Do not specify both all and none in stack")
    if ("all" %in% stack) stack=c("Z","D","Y","Ybin","above","condition")
    if ("none" %in% stack) stack=NULL
    if (("Ybin" %in% stack)&(!"Y" %in% stack)) stop("Cannot stack across bins of Y (stack includes Ybin) if not also stacking across outcomes (stack includes Y)")

    if (("condition" %in% stack)&(sum(c("MW","AHS") %in% test)>0)&(sum(c("simple","BP","K") %in% test)>0)&treetype=="CART") {
      stop("Cannot stack conditions when testing MW or AHS together with simple, BP or K conditions using the CART approach - right hand side variables differ!")
    }

    if (sum(sum(grepl("Z.hat",colnames(data))))) stop("Variable name beginning with Z.hat discovered, reserved for internal use. Please rename.")
    if (sum(sum(grepl("D.hat",colnames(data))))) stop("Variable name beginning with D.hat discovered, reserved for internal use.. Please rename.")
    if (sum(sum(grepl("Q.hat",colnames(data))))) stop("Variable name beginning with Q.hat discovered, reserved for internal use.. Please rename.")
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

    if (test=="simple"&is.null(Y)==FALSE) Y=NULL

    ##Prepare data
    W=colnames(data)[!colnames(data) %in% c(Y,X,Z,D,weight,cluster)]
    XW=c(X,W)

    if (sum(c("MW","AHS") %in% test)>0) {XWY=c(X,W,Y)} else {XWY=XW}

    data=data.table(data)
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

    ##Discretize Z, D and Y into subsets

    ##NB: SHOULD APPLY WEIGHTS WHEN CALCULATING QUANTILES

    if (is.null(Dsubsets)==FALSE) { ##bin treatment
      if (nrow(unique(data[,D,with=F]))>Dsubsets){
        if (gridtypeD=="equidistant") {
          data[,D:=as.numeric(cut(D, breaks = seq(from = min(D) - 0.001, to = max(D) + 0.001, length.out = Dsubsets + 1)))-1,env=list(D=D)]
        }
        else {
          data[,D:=as.numeric(cut(D, breaks = c(-Inf, quantile(D, seq(1/Dsubsets, 1 - 1/Dsubsets, 1/Dsubsets)), Inf)))-1,env=list(D=D)]
        }
      }
    }

    if (is.null(Zsubsets)==FALSE) { ##bin instrument
      if (nrow(unique(data[,Z,with=F]))>Zsubsets){
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
      if (nrow(unique(data[,Yno,with=F]))>Ysubsets) {
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
    if (("Z" %in% stack)&(!"Y" %in% stack)&K>1) stop("Cannot stack across outcomes (stack includes Y) if not also stacking margins of Z (stack includes Z) and Z is multivalued")

    results=c();tunable.Cparams=c()

    if (J>1&(sum(test %in% "simple")!=length(test))) stop("Multivalued treatments and testable conditions other than simple not supported")

    if (J==1&K==1&is.null(X)==TRUE&is.null(W)==TRUE&is.null(Y)==TRUE) {
      stop("Nothing to test with a binary treatment, a binary instrument and no other variables in data.")
    }

    ######################## DETERMINE MARGINS OF D,Z,Y,Ybin,condition,above for stacking or looping ######################

    stackmargins=nostackmargins=c();list=list();index=data.frame()
    if (K>1) { ##margins of Z
      if ("Z" %in% stack) stackmargins=c(stackmargins,"zmargin") else nostackmargins=c(nostackmargins,"Z")
      list=append(list,list(zmargin=1:K))
    }

    if (J>1) { ##margins of D
      if ("D" %in% stack) stackmargins=c(stackmargins,"dmargin") else nostackmargins=c(nostackmargins,"D")
      list=append(list,list(dmargin=1:J))
    }

    if (L>1) { ##Outcomes
      if ("Y" %in% stack) stackmargins=c(stackmargins,"outcome") else nostackmargins=c(nostackmargins,"outcome")
      list=append(list,list(outcome=Y))
    }

    if (length(test)>1) { ##testable conditions
      if ("condition" %in% stack) stackmargins=c(stackmargins,"condition") else nostackmargins=c(nostackmargins,"condition")
      list=append(list,list(condition=test))
    }

    if (length(list)>0) index=do.call(CJ,list)

    if (sum(test %in% c("K","BP"))>0) { ##bins of Y
      if (length(list)==0) {
        index=data.table(ybin=0:(maxlevs-1))
      } else {
        index=index[rep(seq(.N),ifelse(condition %in% c("K","BP"),rep(maxlevs,each=length(Y)),1))]
        index[condition %in% c("K","BP"),ybin:=sequence(maxlevs),by=c(stackmargins[!stackmargins %in% "outcome"],nostackmargins[!nostackmargins %in% "outcome"])]
        #index[!condition %in% c("K","BP"),ybin:=0]
      }
      if ("Ybin" %in% stack) stackmargins=c(stackmargins,"ybin") else nostackmargins=c(nostackmargins,"ybin")
    }

     if (sum(test %in% c("MW","K","BP"))>0) { ##above/below for two-eq. conditions
         if (nrow(index)==0) {
           index=data.table(above=0:1)
        } else {
            if (length(test)>1) {
              index=index[rep(seq(.N),1+(condition %in% c("MW","BP","K")))]
              index[condition %in% c("MW","BP","K"),above:=seq(.N)-1*(condition %in% c("MW","BP","K")),by=c(stackmargins,nostackmargins)]
              #index[!condition %in% c("MW","BP","K"),above:=1,by=c(stackmargins,nostackmargins)]
            } else {
              index=index[rep(seq(.N),2)]
              index[,above:=seq(.N)-1,by=c(stackmargins,nostackmargins)]
            }
        }
          if ("above" %in% stack) stackmargins=c(stackmargins,"above") else nostackmargins=c(nostackmargins,"above")
        }

      if (is.null(c(stackmargins,nostackmargins))==FALSE) {
        setorderv(index,cols=c(nostackmargins,stackmargins))
        setcolorder(index,c(nostackmargins,stackmargins))
        index[,dotest:=seq_len(.N)==.N,by=eval(nostackmargins)]
        loopmax=max(nrow(index),1)
        } else loopmax=1


    ######################## Estimate Z.hat, D.hat and residualize Y (if not stacking) #########

       if (("zmargin" %in% nostackmargins)|length(stackmargins)==0) {  ##estimate Z.hat for all margins
        for (z in 1:K) {
          data[Z==z|Z==z-1,paste0("Z.hat",z):=do.call(regression_forest,append(list(X=.SD[,..X],Y=as.numeric(.SD[,Z]==z),
             num.trees=max(50,num.trees/4),tune.parameters=tune.Zparameters,tune.num.trees=tune.num.trees,
             tune.num.reps=tune.num.reps),Zparameters))$prediction,by=sample,env=list(Z=Z)]
          if (normalize==TRUE) { #Normalize propensity scores
            data[Z==z|Z==z-1,paste0("Z.hat",z):=.N*get(paste0("Z.hat",z))/sum((Z==z)/get(paste0("Z.hat",z))),by=sample,env=list(Z=Z)]
          }
        }
       }

          ##RESIDUALIZE Y separately in each sample of Z=z|Z=z-1 if testing MW or AHS and using Y.res=TRUE
          if ((sum(test %in% c("MW","AHS"))>0)&("outcome" %in% nostackmargins))  {
            for (z in 1:K) {
            for (Yno in Y) {
              data[Z==z|Z==z-1,paste0(eval(Yno),".res",z):=Yvar-do.call(regression_forest,append(list(X=.SD[,..XW],Y=as.numeric(.SD[,Yvar]),
                  num.trees=max(50,num.trees/4),tune.parameters=tune.Yparameters,tune.num.trees=tune.num.trees,
                  tune.num.reps=tune.num.reps),Yparameters))$prediction,by=sample,env=list(Yvar=Yno,Z=Z)]
            }
         }
      }

    ######################## STACK DATA ############################
    if (length(stackmargins)>0) {
    tmpmargins=c()

    ##HANDLE Z LOOPS OR STACKS
      if (! "zmargin" %in% nostackmargins) {
        if ("zmargin" %in% stackmargins) { ##Expand to all margins of Z
          data=data[rep(seq(.N),1+1*(Z>0&Z<K)),env=list(Z=Z)]
          data[,zmargin:=seq(.N)+Z-1*(Z>0),by=id,env=list(Z=Z)]
          data[,Z:=Z>=zmargin,env=list(Z=Z)]
          tmpmargins=c(tmpmargins,"zmargin")
        }
      ##estimate Z.hat for each margin in stacked data
        if (is.null(X)==FALSE) {
        XZ=c(X,Z)
        data[,Z.hat:=do.call(regression_forest,append(list(X=.SD[,..X],Y=as.numeric(.SD[,Z]),
            num.trees=max(50,num.trees/4),tune.parameters=tune.Zparameters,tune.num.trees=tune.num.trees,
            tune.num.reps=tune.num.reps),Zparameters))$prediction,by=c("sample",tmpmargins),.SDcols=XZ,env=list(Z=Z)]
        if (normalize==TRUE) { #Normalize propensity scores
          data[,Z.hat:=.N*Z.hat/sum(Z/Z.hat),by=c("sample",tmpmargins),env=list(Z=Z)]
        }
      } else {
        data[,Z.hat:=(mean(Z)*.N-Z)/(.N-1),by=c("sample",tmpmargins),env=list(Z=Z)] ##leave one out mean
      }
    }

    ##STACK MULTIPLE OUTCOMES
    if (("outcome" %in% stackmargins)) {
      data=melt(data, measure = list(paste0(Y,rep(".bin",length(Y))),Y),value.name = c("Ystack.bin", "Ystack"),variable.name="outcome")
      data[,outcome:=Y,by=c("id",tmpmargins)]
      data[,maxlevsY:=maxlevs,by=c("id",tmpmargins)]
      tmpmargins=c(tmpmargins,"outcome")
      Y="Ystack"
    }

    ##RESIDUALIZE Y in stacked data if testing MW or AHS and using Y.res=TRUE
    if ((sum(test %in% c("MW","AHS"))>0)&Y.res==TRUE)  {
      XY=c(X,Y)
      data[,Y.res:=Y-do.call(regression_forest,append(list(X=.SD[,..XW],as.numeric(.SD[,Y]),
           num.trees=max(50,num.trees/4),tune.parameters=tune.Yparameters,tune.num.trees=tune.num.trees,
           tune.num.reps=tune.num.reps),Yparameters))$prediction,by=c("sample",tmpmargins),env=list(Y=Y)]
    }

    ##STACK MARGINS OF D
    if ("dmargin" %in% stackmargins) { ##Stack data for all margins of D
          data=cbind(CJ(dmargin=(1:J),id=data$id),data[,!"id"])
          data[,D:=D>=dmargin,env=list(D=D)]
          tmpmargins=c(tmpmargins,"dmargin")
        }

    ##Estimate D.hat if test includes "K"
    if ("K" %in% test) {
      XWD=c(X,W,D)
      if (J==1|"D" %in% stackmargins) {
      data[,D.hat:=do.call(regression_forest,append(list(X=.SD[,..XW],
          Y=.SD[,D],num.trees=max(num.trees/4,50),
          tune.parameters=tune.Dparameters,tune.num.trees=tune.num.trees,
          tune.num.reps=tune.num.reps),Dparameters))$prediction,
          by=c("sample",tmpmargins),.SDcols=XWD,env=list(D=D)]
      }
    }

      ##Expand multiple conditions for testing (except K + BP, which has same def of Q - expand later)
        if (length(test)>1) {
          if ("condition" %in% stackmargins) {
            if (sum(test %in% c("BP","K"))==2) testexpand=c(test[!test %in% c("BP","K")],"BPK")
            else testexpand=test
            if (length(testexpand)>1) {
              data=data[rep(seq(.N),length(testexpand))]
              data[,condition:=testexpand,by=c("id",tmpmargins)]
            } else {
              data[,condition:=ifelse(length(test)==2,"BPK",test)]
            }
          tmpmargins=c(tmpmargins,"condition")
            }
        } else {condition=test}

      ##Expand to a=0,1 for BP, K and MW conditions
        if (sum(test %in% c("MW","BP","K"))>0) {
          if (!"above" %in% nostackmargins) {
            if (!"condition" %in% nostackmargins) {
              data=data[rep(seq(.N),1+condition %in% c("BPK","MW","BP","K"))]
              data[condition %in% c("BPK","MW","BP","K"),above:=seq(.N)-1*(condition %in% c("BPK","MW","BP","K")),by=c("id",tmpmargins)]
            } else {
              data=data[rep(seq(.N),2)]
              data[,above:=0:1,by=c("id",tmpmargins)]
            }
            tmpmargins=c(tmpmargins,"above")
          }
        }


      ##Expand to all groups of Ybin for BP, K conditions
         if (sum(test %in% c("BP","K"))>0) {
           if (!"outcome" %in% nostackmargins) {
            if (!"condition" %in% nostackmargins) {
              data=data[rep(seq(.N), 1+(maxlevsY-1)*(condition %in% c("BPK","K","BP")))]
              data[condition %in% c("BP","K","BPK"),ybin:=(1:maxlevsY),by=c("id",tmpmargins)]
            } else {
              data=data[rep(seq(.N), length(levsY))]
              data[,ybin:=levsY,by=c("id",tmpmargins)]
            }
            tmpmargins=c(tmpmargins,"ybin")
           }
         }

      ##Create outcome variable Q in stacked data

      if (is.null(nostackmargins)==TRUE) {
      if (sum(test %in% c("simple","AHS"))>0) data[condition %in% c("simple","AHS"),Q:=D,env=list(D=D)]
      if (sum(test %in% c("BPK","K","BP"))>0) data[condition %in% c("BPK","BP","K"),Q:=above*(D*(name==ybin))-(1-above)*(1-D)*(name==ybin),env=list(D=D,name=paste0(Y,".bin"))]
      if ("MW" %in% test) data[condition=="MW",Q:=above*((1-Z.hat)*D*Z-Z.hat*D*(1-Z))+(1-above)*(Z.hat*(1-D)*(1-Z)-(1-Z.hat)*(1-D)*Z),env=list(Z=Z,D=D)]
      if ("K" %in% test) {
          data[condition %in% c("BPK","K"),D.hat:=D.hat*above+(1-D.hat)*(1-above)]
          data[condition %in% c("BPK","K"),D:=D*above+(1-D)*(1-above),env=list(D=D)]
          ##reverse treatment indicator for above==0 and condition=="K"
      }

      ##drop groups/margins/conditions/outcomes with constant Q
      data[,rQ:=max(Q)-min(Q),by=c("sample",tmpmargins)]
      data=data[,rQ:=min(rQ),by=tmpmargins]
      data[,rQ:=NULL]

      ##Estimate Q.hat in stacked data
      if (sum(test %in% c("simple","BP","K"))>0) {
        SDcols=c(XW,"Q","id")
        data[condition %in% c("simple","BPK","BP","K"),Q.hat:=do.call(regression_forest,append(list(X=.SD[,..XW],Y=.SD[,Q],
            num.trees=max(num.trees/4,50),tune.parameters=tune.Qparameters,
            tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps),Qparameters))$prediction,
            by=c("sample",tmpmargins),.SDcols=SDcols]
      }
        if ("AHS" %in% test) {
        SDcols=c(X,W,Y,"Q","id")
        XWY=c(X,W,Y)
        data[condition=="AHS",Q.hat:=do.call(regression_forest,append(list(X=.SD[,..XWY],Y=.SD[,Q],
           num.trees=max(num.trees/4,50),tune.parameters=tune.Qparameters,
           tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps),Qparameters))$prediction,
      by=c("sample",tmpmargins),.SDcols=SDcols]
        }

      ##Split BP and K conditions if doing both
      if ((sum(test %in% c("BP","K"))>0)&("condition" %in% stack)) {
          if (sum(test %in% c("BP","K"))==2) data=data[rep(seq(.N),1+1*(condition=="BPK"))]
          data[condition=="BPK",condition:=test[test %in% c("BP","K")],by=c("id",tmpmargins)]
      }
      }

    }


    ########## LOOP OVER INDEX OF MARGINS TO ESTIMATE #####################
    if (saveforest==T) forests=list()
    for (l in 1:loopmax) {
      if (K>1) z=index[l,zmargin] else z=1
      if (J>1) d=index[l,dmargin] else d=1
      if (sum(c("BP","K","MW") %in% test)>0) a=index[l,above] else a=NA
      if (sum(c("BP","K") %in% test)>0) y=index[l,ybin] else y=NA
      if (length(test)>1) cond=index[l,condition] else cond=test
      if (L>1) outcome=index[l,outcome] else outcome=Y

      if (cond %in% c("AHS","MW")) vars=XWY
      else vars=XW

      data[,use:=0]
      if (is.null(stackmargins)==FALSE) data[index[l],use:=1,on=stackmargins] else data[,use:=1]
      if ("zmargin" %in% nostackmargins) data[!(Z==z|Z==z-1),use:=0,env=list(Z=Z)]


      if (("Z" %in% nostackmargins)|length(stackmargins)==0) { ##Fix Z.hat
        data[use==1,Z.hat:=get(paste0("Z.hat",z))]
        }

      if ((cond %in% c("MW","AHS"))&(("outcome" %in% nostackmargins))) { #fix Y.res
        data[use==1,Y.res:=get(paste0(outcome,".res",z))]
      }

      if (is.null(nostackmargins)==FALSE|(is.null(c(nostackmargins,stackmargins))==TRUE)) { # define Q and estimate Q.hat, estimate D.hat if relevant
        ####NEEED TO ESTIMATE D.hat if using KITAGAWA + no stacking!!!

        ##define Q if not stacking
        if (cond=="MW") {
          data[use==1,Q:=a*((1-Z.hat)*(D>=d)*(Z==z)-Z.hat*D*(1-(Z==z)))+(1-a)*(Z.hat*(1-(D>=d))*(1-(Z==z))-(1-Z.hat)*(1-(D>=d))*(Z==z)),env=list(Z=Z,D=D)]
           } else if (cond %in% c("BP","K")) {
          data[use==1,Q:=a*(D>=d)*(name==y)-(1-a)*(1-(D>=d))*(name==y),env=list(D=D,name=paste0(outcome,".bin"))]
        } else {
          data[use==1,Q:=1*(D>=d),env=list(D=D)]
        }

        if (cond=="K") {
          data[use==1,D.hat:=do.call(regression_forest,append(list(X=.SD[,..XW],
            Y=.SD[,D],num.trees=max(num.trees/4,50),
            tune.parameters=tune.Dparameters,tune.num.trees=tune.num.trees,
            tune.num.reps=tune.num.reps),Dparameters))$prediction,
             by="sample",.SDcols=XWD,env=list(D=D)]
        }

        ##Determine if no variation in Q in either sample
        if (sd(data[sample==1&use==1,Q])==0|sd(data[sample==2&use==1,Q])==0) novar=TRUE else novar=FALSE

        ##estimate Q.hat
        if (novar==FALSE) {
          data[use==1,Q.hat:=do.call(regression_forest,append(list(X=.SD[,..vars],Y=.SD[,Q],
          num.trees=max(num.trees/4,50),tune.parameters=tune.Qparameters,
          tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps),Qparameters))$prediction,
          by=c("sample")]
        }

      } else {
        novar=FALSE
      }

      ##estimate causal/instrumental/regression forests

      if ("zmargin" %in% stackmargins) z=1; if ("dmargin" %in% stackmargins) d=1
      if (novar==FALSE) {
        for (i in 1:2) {
            if (cond=="K") {
            forest=do.call(instrumental_forest,append(list(X=data[use==1&sample==i,vars,with=F],
                Y=data[use==1&sample==i,Q],Y.hat=data[use==1&sample==i,Q.hat],
                Z=data[use==1&sample==i,Z==z,env=list(Z=Z)],Z.hat=data[use==1&sample==i,Z.hat],
                W=data[use==1&sample==i,D>=d,env=list(D=D)],W.hat=data[use==1&sample==i,D.hat],
                compute.oob.predictions=compute.oob.predictions,num.trees=max(num.trees,50),
                tune.parameters=tune.Cparameters,tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps),Cparameters))
          } else if (cond=="MW") {
            forest=do.call(regression_forest,append(list(X=data[use==1&sample==i,vars,with=F],
                Y=data[use==1&sample==i,Q],
                compute.oob.predictions=compute.oob.predictions,num.trees=max(num.trees,50),
                tune.parameters=tune.Cparameters,tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps),Cparameters))
          } else  {
            forest=do.call(causal_forest,append(list(X=data[use==1&sample==i,vars,with=F],
                Y=data[use==1&sample==i,Q],Y.hat=data[use==1&sample==i,Q.hat],
                W=data[use==1&sample==i,Z==z,env=list(Z=Z)],W.hat=data[use==1&sample==i,Z.hat],
                compute.oob.predictions=compute.oob.predictions,num.trees=max(num.trees,50),
                tune.parameters=tune.Cparameters,tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps),Cparameters))
          }
          if (saveforest==T) forests=append(forests,list(forest))

          if (cond!="MW") data[use==1&sample==i,scores:=get_scores(forest)]
          else data[use==1&sample==i,scores:=Q]
          if (treetype=="forest") {
            if (se==FALSE) {
              data[use==1&sample==i,tau:=forest$predictions]
              data[use==1&sample!=i,taupred := predict(forest,newdata=data[use==1&sample!=i,vars,with=F])$prediction]
            } else {
              data[use==1&sample==i,c("tau","var_tau"):=predict(forest,estimate.variance=T)[,1:2]]
              data[use==1&sample!=i,c("taupred","var_taupred") := predict(forest,newdata=data[use==1&sample!=i,vars,with=F],estimate.variance=T)[1:2]]
            }
          }
          if (tune.Cparameters!="none"&((tunetype=="all")|tunetype=="one"&is.null(tunable.Cparams)==T)) {
            tunable.Cparams=rbind(tunable.Cparams,c(dmargin=d,ybin=y,above=a,zmargin=z,condition=cond,sample=i,unlist(forest$tunable.params)))
            if (tunetype=="one") {
              Cparameters=forest$tunable.params
            }
            } else Cparameters=forest$tunable.params
         } ##end i loop
      }

      if (nrow(index)>0) dotest=index[l,dotest] else dotest=TRUE
      if (dotest==TRUE) { ##TEST HERE!
        if (treetype=="forest") { ##test using forest approach
          if (se==TRUE) {
            data[,tau:=tau/sqrt(var_tau)]
            data[,taupred:=taupred/sqrt(var_taupred)]
          }
          setorderv(data,cols=c("sample","tau"))

          data[(Z==z|Z==z-1),N:=seq_len(.N),by=sample,env=list(Z=Z)]
          data[(Z==z|Z==z-1),m:=cumsum(scores)/N,by=sample,env=list(Z=Z)]
            if (is.null(cluster)==TRUE) {
              data[(Z==z|Z==z-1),m2:=cumsum(scores^2)/seq_len(.N),by=sample,env=list(Z=Z)]
              data[(Z==z|Z==z-1),G:=N,env=list(Z=Z)]
              data[(Z==z|Z==z-1),se:=sqrt((m2-m^2)/(N-1)),by=sample,env=list(Z=Z)]
            } else {
              data[(Z==z|Z==z-1),dyiyj:=2*cumsum(scores)*scores-scores^2,by=c(cluster,"sample"),env=list(Z=Z)]
              data[(Z==z|Z==z-1),Ng:=seq_len(.N),by=c(cluster,"sample"),env=list(Z=Z)]
              data[(Z==z|Z==z-1),G:=cumsum(!duplicated(get(cluster))),by=sample,env=list(Z=Z)]
              data[(Z==z|Z==z-1),dmy:=cumsum(scores)+scores*(Ng-1),by=c(cluster,"sample"),env=list(Z=Z)]
              data[(Z==z|Z==z-1),se:=sqrt((G/(G-1))*(cumsum(dyiyj)-2*m*cumsum(dmy)+m^2*cumsum(2*Ng-1)))/seq_len(.N),by=sample,env=list(Z=Z)]
            }

            data[(Z==z|Z==z-1),t:=m/se,env=list(Z=Z)]

            res=data[G>=minsize&(Z==z|Z==z-1), .SD[which.min(t)], by = sample,.SDcols=c("G","N","m","se","t","tau"),env=list(Z=Z)]

            if (is.null(cluster)==FALSE) clust=as.formula(paste0("~",cluster)) else clust=NULL

            browser()
            if (nrow(data[(Z==z|Z==z-1)&sample==2&taupred<=res[1,tau],env=list(Z=Z)])>1) {
              if (sd(data[(Z==z|Z==z-1)&sample==2&taupred<=res[1,tau],scores,env=list(Z=Z)])>0) {
                fe=feols(scores~1,data=data[(Z==z|Z==z-1)&sample==2&taupred<=res[sample==1,tau],env=list(Z=Z)],cluster=clust)
                test1=c(length(unique(data[(Z==z|Z==z-1)&sample==2&taupred<=res[sample==1,tau],G,env=list(Z=Z)])),fe$nobs,coeftable(fe)[1:3])
                } else test1=rep(NA,5)
              } else test1 = rep(NA,5)
            if (nrow(data[(Z==z|Z==z-1)&sample==1&taupred<=res[2,tau],env=list(Z=Z)])>1) {
              if (sd(data[(Z==z|Z==z-1)&sample==1&taupred<=res[2,tau],scores,env=list(Z=Z)])>0) {
                fe=feols(scores~1,data=data[(Z==z|Z==z-1)&sample==1&taupred<=res[sample==2,tau],env=list(Z=Z)],cluster=clust)
                test2=c(length(unique(data[(Z==z|Z==z-1)&sample==1&taupred<=res[sample==2,tau],G,env=list(Z=Z)])),fe$nobs,coeftable(fe)[1:3])
              } else test2 = rep(NA,5)
            } else test2 = rep(NA,5)

            res=cbind(res,rbind(test1,test2))
            if (length(nostackmargins)>0) results=rbind(results,cbind(index[l,..nostackmargins],res)) else results=rbind(results,res)

        } else { ##test using CART approach
          for (i in 1:2) {
            treevec=c(XWY,stackmargins,"sample","id","scores")
            esttree=esttree(data=data[(Z==z|Z==z-1),..treevec,env=list(Z=Z)],testsample=i,cp=cp,maxrankcp=maxrankcp,alpha=alpha,prune=prune,minsize=minsize,preselect=preselect,cluster=cluster)

            if (length(nostackmargins)>0) {
              results=rbind(results,cbind(index[l,..nostackmargins],sample=rep(i,nrow(esttree$res)),esttree$res))
            } else {
              results=rbind(results,cbind(sample=rep(i,nrow(esttree$res)),esttree$res))
            }
            ##trees[[i]]=esttree$tree

        } ##end i loop
      } ##end CART loop
    } ##end test loop
    } ##end margins loop


      ##NAMING and postprocessing of results

      if (treetype=="forest") {
                taunames="tau cutoff";leafnames=NULL
                } else {
                  taunames="relevant";leafnames="leaf"
                }
      colnames(results)=c(nostackmargins,"sample",leafnames,"G.train","N.train","coef.train","stderr.train","t.train",taunames,"G.est","N.est","coef.est","stderr.est","t.est")

      results[is.na(t.est)==FALSE,p.raw:=pnorm(t.est)]
      for (m in c("holm","hochberg","BH","BY")) {
        results[is.na(t.est)==FALSE,paste0("p.",m):=p.adjust(p.raw,method=m)]
      }

      minp=apply(results[is.na(t.est)==FALSE,c("p.raw","p.holm","p.hochberg","p.BH","p.BY")],2,min)
      return=list(results=results,minp=minp)
      #if (treetype=="CART"){return=append(return,list(trees=trees))}

      time=proc.time()-time
      if (nrow(index)>0) return=append(return,list(margins=index[,1:(ncol(index)-1)]))
      if (tune.Cparameters=="none") {
        Cout=unlist(Cparameters)
      } else Cout=tunable.Cparams
      return=append(return,list(time=time,Cparameters=Cout))
      if (saveforest==T) return=append(return,list(forests=forests))
      return(return)
  }

##HELPER FUNCTIONS
 esttree=function(data,testsample,cp,maxrankcp,alpha,prune,minsize,preselect,cluster=NULL){
    tree=rpart(scores~.,data=as.data.frame(data[sample==testsample,!c("id","sample")]),method="anova",cp=cp,minbucket=minsize) # run tree with transformed outcome

    maxrankcp=min(maxrankcp, length(tree$cp[, 4]))
    maxcp=tree$cp[maxrankcp, 1]
    if (is.null(cluster)==FALSE) clust=as.formula(paste0("~",cluster)) else clust=NULL
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

    res=data[is.na(leaf)==FALSE,data.table(cbind(G=max(G),N=.N,feols(scores~1,data=.SD,vcov=clust)$coeftable)),by=.(leaf,sample)]
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
