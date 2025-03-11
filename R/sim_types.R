    library(montest)
    library(foreach)
    library(doParallel)

    n=2000
    R=1000

    cols=c("dgp","treetype","stack")
    speclist= data.frame(matrix(nrow = 0, ncol = length(cols)))
    colnames(speclist)=cols

    for (dgp in 10:20) {
      for (treetype in c("CART","forest")) {
        for (stack in c("none","all")) {
          if (dgp<=13&stack=="all") next #nothing to stack
            speclist[nrow(speclist)+1,cols]=c(dgp,treetype,stack)
          }
      }
    }

    if (restart==F) {
      start=max(readRDS("sim_types.Rdata")[,"rep"])+1
    } else start=1

    cl <- makeCluster(min(nrow(speclist),detectCores()/2),outfile="outfile.out")
    registerDoParallel(cl)
    pb <- txtProgressBar(max = R,style=3,initial=start-1)
    for (r in (1:R)) {
      sim_types = foreach (i = 1:nrow(speclist),.combine=rbind,.packages=c("montest","fixest","grf","rpart","treeClust","mvtnorm","data.table","rpart","LATEtest")) %dopar% {
        set.seed(r)
        dgp=speclist[i,"dgp"]
        if (dgp==10) data=simtypes_fs(J=1,K=1,n=n,complier.probs=matrix(0,nrow=1,ncol=1))
        if (dgp==11) data=simtypes_fs(J=1,K=1,n=n,complier.probs=matrix(0.15,nrow=1,ncol=1))
        if (dgp==12) data=simtypes_fs(J=1,K=1,n=n,complier.probs=matrix(0.15,nrow=1,ncol=1),defier.probs=matrix(0.25,nrow=1,ncol=1),condition="Xvar1<qnorm(0.15)")
        if (dgp==13) data=simtypes_fs(J=1,K=1,n=n,complier.probs=matrix(0.15,nrow=1,ncol=1),defier.probs=matrix(0.25,nrow=1,ncol=1),condition="abs(Xvar1)>qnorm(0.925)")
        if (dgp==14) data=simtypes_fs(J=1,K=2,n=n,complier.probs=matrix(0.115,nrow=1,ncol=2),defier.probs=matrix(0.3,nrow=1,ncol=2),condition="Xvar1<qnorm(0.15)")
        if (dgp==15) data=simtypes_fs(J=2,K=1,n=n,complier.probs=matrix(0.115,nrow=2,ncol=1),defier.probs=matrix(0.2,nrow=2,ncol=1),condition="Xvar1<qnorm(0.15)")
        if (dgp==16) data=simtypes_fs(J=2,K=2,n=n,complier.probs=matrix(0.07,nrow=2,ncol=2),defier.probs=matrix(0.125,nrow=2,ncol=2),condition="Xvar1<qnorm(0.15)")
        if (dgp==17) data=simtypes_fs(J=1,K=2,n=n,complier.probs=matrix(c(0.225,0),nrow=1,ncol=2),defier.probs=matrix(c(0,0.125),nrow=1,ncol=2))
        if (dgp==18) data=simtypes_fs(J=2,K=1,n=n,complier.probs=matrix(c(0.225,0),nrow=2,ncol=1),defier.probs=matrix(c(0,0.09),nrow=2,ncol=1))
        if (dgp==19) data=simtypes_fs(J=2,K=2,n=n,complier.probs=matrix(c(rep(0.1,3),0),nrow=2,ncol=2),defier.probs=matrix(c(rep(0,3),0.125),nrow=2,ncol=2))
        if (dgp==20) data=simtypes_fs(J=2,K=2,n=n,complier.probs=matrix(c(rep(0.06,4)),nrow=2,ncol=2),defier.probs=matrix(c(rep(0,3),0.3),nrow=2,ncol=2),condition="Xvar1<qnorm(0.15)")

        est=try(montest(data=data,D="D",Z="Z",X=c("Xvar1","Xvar2","Xvar3"),stack=speclist[i,"stack"],treetype=speclist[i,"treetype"],test="simple"))
        if (!inherits(est, "try-error")) {
          out = data.frame(rep=r,speclist[i,],t(est$minp),t(est$time[1:3]),timestamp=Sys.time())
          out
        } else {
          out=data.frame(r,speclist[i,],t(rep(NA,8)),Sys.time())
          colnames(out)=c("rep",colnames(speclist),"p.raw","p.holm","p.hochberg","p.BH","p.BY","user.self","sys.self","elapsed","timestamp")
          out
        }
      }

      if (restart==F|r>1) sim_types=rbind(data.frame(readRDS("sim_types.Rdata")),sim_types)
      saveRDS(sim_types,file="sim_types.Rdata")
      setTxtProgressBar(pb,r)
    }

  stopCluster(cl)
