

############
BAYESIAN SHRINKAGE STUFF

##Empirical Bayes shrinkage of predictions
if (shrink==TRUE) { ##NEED TO ACCEPTS WEIGHTS
  calc_stats <- function(DT) {
    mu  <- mean(DT[["scores"]],   na.rm = TRUE)
    v1  <- var(DT[["pred"]],      na.rm = TRUE)
    e1  <- mean(DT[["pred_var"]], na.rm = TRUE)
    n1  <- sum(!is.na(DT[["pred"]]))
    v2  <- var(DT[["pred_o"]],      na.rm = TRUE)
    e2  <- mean(DT[["pred_o_var"]], na.rm = TRUE)
    n2  <- sum(!is.na(DT[["pred_o"]]))
    sig2_tau <- pmax(weighted.mean(c(v1 - e1, v2 - e2), w = c(n1, n2), na.rm = TRUE),0)
    list(mu = mu, sig2_tau = sig2_tau, e1_bar = e1, e2_bar = e2)
  }

  if (is.null(margins)) {
    # no grouping: one-row stats
    tmp <- calc_stats(data)
    stats <- data.table(mu = tmp$mu, sig2_tau = tmp$sig2_tau,
                        e1_bar = tmp$e1_bar, e2_bar = tmp$e2_bar)
    # attach constants
    data[, `:=`(mu = stats$mu[1L],
                sig2_tau = stats$sig2_tau[1L],
                e1_bar = stats$e1_bar[1L],
                e2_bar = stats$e2_bar[1L])]
  } else {
    # grouped stats by margins
    stats <- data[, {
      tmp <- calc_stats(.SD)
      .(mu = tmp$mu, sig2_tau = tmp$sig2_tau, e1_bar = tmp$e1_bar, e2_bar = tmp$e2_bar)
    }, by = margins]
    # join back on those keys
    data <- stats[data, on = margins]
  }

  # 2) Join and create per-obs noise vars
  data[, `:=`(
    se2_pred   = fifelse(is.na(pred_var),   e1_bar, pred_var),
    se2_pred_o = fifelse(is.na(pred_o_var), e2_bar, pred_o_var)
  )]

  # 3) Weights
  data[, `:=`(
    w_pred   = fifelse(sig2_tau + pmax(se2_pred,   0) > 0,
                       sig2_tau / (sig2_tau + pmax(se2_pred,   0)), 0),
    w_pred_o = fifelse(sig2_tau + pmax(se2_pred_o, 0) > 0,
                       sig2_tau / (sig2_tau + pmax(se2_pred_o, 0)), 0)
  )]

  # 4) Shrunken estimates
  data[, `:=`(
    pred   = mu + ((1-shrink.alpha)*w_pred+shrink.alpha)   * (pred   - mu),
    pred_o = mu + ((1-shrink.alpha)*w_pred_o+shrink.alpha)  * (pred_o - mu)
  )]

  # (optional) tidy up helpers
  data[, c("e1_bar","e2_bar","se2_pred","se2_pred_o","pred_var","pred_o_var","sig2_tau","w_pred","w_pred_o") := NULL]
}


########3
##########   If not stacking: Loop #####################
if (stack==FALSE) {

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

  if (sum(test %in% c("MW","K","BP"))>0) { ##above/below for two-eq. conditions
    margins=c(margins,"above")
    if (nrow(index)==0) index=data.table(above=0:1) else {
      if (length(test)>1) {
        index=index[rep(seq(.N),1+(condition %in% c("MW","BP","K")))]
        index[condition %in% c("MW","BP","K"),above:=seq(.N)-1*(condition %in% c("MW","BP","K")),by=margins]
      } else {
        index=index[rep(seq(.N),2)]
        index[,above:=seq(.N)-1,by=margins]
      }
    }
  }

  if (is.null(margins)==FALSE) {
    setorderv(index,cols=margins)
    setcolorder(index,margins)
    index[,dotest:=seq_len(.N)==.N,by=eval(margins)]
    loopmax=max(nrow(index),1)
  } else loopmax=1

  ################# PRE-ESTIMATE Q.hat, D.hat, Y.res etc so not to do it over and over ##############


  ################# RUN THE LOOP, TEST AND ADD TO RESULTS ###########################################
  for (l in 1:loopmax) {
    if (K>1) z=index[l,zmargin] else z=1
    if (J>1) d=index[l,dmargin] else d=1
    if (sum(c("BP","K","MW") %in% test)>0) a=index[l,above] else a=NA
    if (sum(c("BP","K") %in% test)>0) y=index[l,ybin] else y=NA
    if (length(test)>1) cond=index[l,condition] else cond=test
    if (L>1) outcome=index[l,outcome] else outcome=Y


    ##if not stacking: set Z.hat, define Q, estimate Q.hat
    if (stack==FALSE) {
      if (K>1) data[Z==z|Z==z-1,Z.hat:=get(paste0("Z.hat",z)),env=list(Z=Z)]

      if (cond %in% c("simple","AHS")) data[Z==z|Z==z-1,Q:=D,env=list(D=D)]
      else if (cond %in% c("BPK","K","BP")) data[Z==z|Z==z-1,Q:=a*((D>=d)*(name==y))-(1-a)*(1-(D>=d))*(name==y),env=list(D=D,name=paste0(Y,".bin"))]
      else if ("MW" %in% test) data[Z==z|Z==z-1,Q:=a*((1-Z.hat)*(D>=d)*(Z>=z)-Z.hat*(D>=d)*(1-(Z>=z)+(1-a)*(Z.hat*(1-(D>=d))*(1-(Z>=z))-(1-Z.hat)*(1-(D>=d))*(Z>=z)),env=list(Z=Z,D=D)]
      if ("K" %in% test) {
        data[Z==z|Z==z-1,D.hat:=D.hat*a+(1-D.hat)*(1-a)]
        data[Z==z|Z==z-1,D:=D*a+(1-D)*(1-a),env=list(D=D)]
        ##reverse treatment indicator for above==0 and condition=="K"
      }
      if (cond %in% c("AHS","MW")) {
        data[Z==z|Z==z-1,Q.hat:=do.call(regression_forest,append(list(X=.SD[,..XWY],Y=.SD[,Q],
                                                                      tune.parameters=tune.Qparameters,sample.weights=.SD[,weight],clusters=.SD[,cluster]),Qparameters,treeopts))$prediction,env=list(weight=weight,cluster=cluster),by=sample]
      } else {
        data[Z==z|Z==z-1,Q.hat:=do.call(regression_forest,append(list(X=.SD[,..XW],Y=.SD[,Q],
                                                                      tune.parameters=tune.Qparameters,sample.weights=.SD[,weight],clusters=.SD[,cluster]),Qparameters,treeopts))$prediction,env=list(weight=weight,cluster=cluster),by=sample]
      }

    }

    ##IF Q depends on loop, define Q and estimate Q.hat





    if ("condition" %in% nostackmargins)|("above" %in% nostackmargins) {
      if (index{l,condition]=="MW") {
        data[use==1,Q:=a*((1-Z.hat)*(D>=d)*(Z==z)-Z.hat*D*(1-(Z==z)))+(1-a)*(Z.hat*(1-(D>=d))*(1-(Z==z))-(1-Z.hat)*(1-(D>=d))*(Z==z)),env=list(Z=Z,D=D)]
      } else if (index[l,condition] %in% c("BP","K")) {
        data[use==1,Q:=a*(D>=d)*(name==y)-(1-a)*(1-(D>=d))*(name==y),env=list(D=D,name=paste0(outcome,".bin"))]
      } else {
        data[use==1,Q:=1*(D>=d),env=list(D=D)]
      }
      }
      ##IF Z not in stackmargn & K>1, set correct Z.hat. and define Z, set correct Y.res, D.hat if relevant

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

      if (is.null(nostackmargins)==FALSE|(is.null(c(nostackmargins,stackmargins))==TRUE)) {
        # define Q and estimate Q.hat, estimate D.hat if relevant
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
                                                                   Y=.SD[,D],tune.parameters=tune.Dparameters,sample.weights=.SD[,weight],clusters=.SD[,cluster]),Dparameters,treopts))$prediction,
               by="sample",env=list(D=D,cluster=cluster,weight=weight)]
        }

        ##Determine if no variation in Q in either sample
        if (sd(data[sample==1&use==1,Q])==0|sd(data[sample==2&use==1,Q])==0) novar=TRUE else novar=FALSE

        ##estimate Q.hat
        if (novar==FALSE) {
          data[use==1,Q.hat:=do.call(regression_forest,append(list(X=.SD[,..vars],Y=.SD[,Q],
                                                                   tune.parameters=tune.Qparameters,sample.weights=.SD[,weight],clusters=.SD[,cluster]),Qparameters,treeopts))$prediction,
               by=c("sample"),env=list(cluster=cluster,weight=weight)]
        }

      } else {
        novar=FALSE
      }

      ######## STOPPED HERE - NEED TO PERFORM BY MARGIN INSTEAD OF LOOPING!!
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
                                                           tune.parameters=tune.Cparameters,tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps, sample.weights=data[use==1&sample==i,weight],clusters=data[use==1&sample==i,cluster]),Cparameters))
          } else if (cond=="MW") {
            forest=do.call(regression_forest,append(list(X=data[use==1&sample==i,vars,with=F],
                                                         Y=data[use==1&sample==i,Q],
                                                         compute.oob.predictions=compute.oob.predictions,num.trees=max(num.trees,50),
                                                         tune.parameters=tune.Cparameters,tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps,
                                                         sample.weights=data[use==1&sample==i,weight],clusters=data[use==1&sample==i,cluster]),Cparameters))
          } else  {
            forest=do.call(causal_forest,append(list(X=data[use==1&sample==i,vars,with=F],
                                                     Y=data[use==1&sample==i,Q],Y.hat=data[use==1&sample==i,Q.hat],
                                                     W=data[use==1&sample==i,Z==z,env=list(Z=Z)],W.hat=data[use==1&sample==i,Z.hat],
                                                     compute.oob.predictions=compute.oob.predictions,num.trees=max(num.trees,50),
                                                     tune.parameters=tune.Cparameters,tune.num.trees=tune.num.trees,tune.num.reps=tune.num.reps,
                                                     sample.weights=data[use==1&sample==i,weight],clusters=data[use==1&sample==i,cluster]),Cparameters))
          }

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


          if (nrow(data[(Z==z|Z==z-1)&sample==2&taupred<=res[1,tau],env=list(Z=Z)])>1) {
            if (sd(data[(Z==z|Z==z-1)&sample==2&taupred<=res[1,tau],scores,env=list(Z=Z)])>0) {
              fe=feols(scores~1,data=data[(Z==z|Z==z-1)&sample==2&taupred<=res[sample==1,tau],env=list(Z=Z)],cluster=clust,weights=weight)
              test1=c(length(unique(data[(Z==z|Z==z-1)&sample==2&taupred<=res[sample==1,tau],G,env=list(Z=Z)])),fe$nobs,coeftable(fe)[1:3])
            } else test1=rep(NA,5)
          } else test1 = rep(NA,5)
          if (nrow(data[(Z==z|Z==z-1)&sample==1&taupred<=res[2,tau],env=list(Z=Z)])>1) {
            if (sd(data[(Z==z|Z==z-1)&sample==1&taupred<=res[2,tau],scores,env=list(Z=Z)])>0) {
              fe=feols(scores~1,data=data[(Z==z|Z==z-1)&sample==1&taupred<=res[sample==2,tau],env=list(Z=Z)],cluster=clust,weights=weight)
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




    res=cbind(res,data[pred_o<=tau, { ##PERFORM TEST
      # drop NAs in the target within this group
      n <- length(scores)
      ybar <- mean(scores)

      # cluster summaries
      clsum <- data.table(cluster = cluster, scores = scores)[, .(ng = .N, ybar_g = mean(scores)), by = cluster]
      G <- nrow(clsum)

      # Liangb
