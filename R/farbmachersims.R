ppool=c()
p=c()
pcf=c()
pnocf=c()
n=2000
pb <- txtProgressBar(min = 0, max = 300, style = 3)
for (r in 1:300) {
  set.seed(r)
  data=fct_datasim(setup="A", dgp=2)
  pnocf=c(pnocf,montest(data=data,D="D",Z="Z",X=c("Xvar1","Xvar2","Xvar3"),test="simple",pool=TRUE,inner.folds=5,crossfit.forest=FALSE)$minp)
  pcf=c(pcf,montest(data=data,D="D",Z="Z",X=c("Xvar1","Xvar2","Xvar3"),test="simple",pool=TRUE,inner.folds=5)$minp)
  ppool=c(ppool,montest(data=data,D="D",Z="Z",X=c("Xvar1","Xvar2","Xvar3"),test="simple",pool=TRUE,inner.folds=5)$minp)
  p=c(p,montest(data=data,D="D",Z="Z",X=c("Xvar1","Xvar2","Xvar3"),test="simple",inner.folds=5,pool=FALSE)$minp[4])
  setTxtProgressBar(pb, r)
}
close(pb)


 #setwd("montest")
  #setwd("../")
  #library(devtools)
  #devtools::install()
  library(montest)
  source("dgp_setup_farbmacher.R")
  library(foreach)
  library(doParallel)

  subsets <- 4
  siglevel <- 0.05
  n=3000
  R=1000

  cols=c("dgp","dtype","montest","farbspec","test","stack","treetype")
  speclist= data.frame(matrix(nrow = 0, ncol = length(cols)))
  colnames(speclist)=cols

  i=0
  for (dgp in 0:5) {
    for (dtype in c("A","B")) {
      speclist[nrow(speclist)+1,cols]=c(dgp,dtype,0,1,"BP","none","CART")
      speclist[nrow(speclist)+1,cols]=c(dgp,dtype,1,1,"BP","none","CART")
      for (test in c("simple","BP","MW","AHS","all")) {
        for (stack in c("none","all")) {
          if (test %in% c("simple","AHS")&stack=="all") next ##nothing to stack
          for (treetype in c("CART","forest")) {
            if (test=="all"&stack=="all"&treetype=="CART") next ##Cannot stack with conditions that has Y on both sides + CART approach  for (test in c("simple","BP","MW","AHS","all")) {
            speclist[nrow(speclist)+1,cols]=c(dgp,dtype,1,0,test,stack,treetype)
            i=i+1
          }
        }
      }
    }
  }




if (restart==F) {
  start=max(readRDS("farbmachersims.Rdata")[,"rep"])+1
  } else start=1

cl <- makeCluster(min(nrow(speclist),detectCores()/2),outfile="outfile.out")
registerDoParallel(cl)

pb <- txtProgressBar(max = R,style=3,initial=start-1)

for (r in start:R) {
  farbmachersims = foreach (i = 1:nrow(speclist),.combine=rbind,.packages=c("montest","fixest","LATEtest","grf","rpart","treeClust","mvtnorm","data.table","rpart")) %dopar% {
    set.seed(r)
    data <- fct_datasim(setup=speclist[i,"dtype"], dgp=speclist[i,"dgp"])
      if (speclist[i,"montest"]==0) { ##LATEtest
        time=proc.time()
        farb=LATEtest(data=data, covars=c("Xvar1","Xvar2","Xvar3"), subsets=subsets, alpha=siglevel)
        time=proc.time()-time
        out=data.frame(r,speclist[i,],1-pnorm(farb$results$Tmax[2]),farb$results$pvalue_Holm[,2],NA,farb$results$pvalue_BenjHochberg[,2],farb$results$pvalue_BenjYekutieli[,2],t(time[1:3]),timestamp=Sys.time())
        colnames(out)=c("rep",colnames(speclist),"p.raw","p.holm","p.hochberg","p.BH","p.BY","user.self","sys.self","elapsed","timestamp")
        out
      } else if (speclist[i,"farbspec"]==1) {    ##COMPARABLE MONTEST
        est=try(montest(data=data,D="D",Z="Z",Y="Y",X=c("Xvar1","Xvar2","Xvar3"),stack="none",treetype="CART",test="BP",Cparameters=list(min.node.size=100),minsize=200))
      if (!inherits(est,"try-error")) {
        data.frame(rep=r,speclist[i,],t(est$minp),t(est$time[1:3]),timestamp=Sys.time())
         } else {
           out=data.frame(rep=r,speclist[i,],t(rep(NA,8)),timestamp=Sys.time())
           colnames(out)=c("rep",colnames(speclist),"p.raw","p.holm","p.hochberg","p.BH","p.BY","user.self","sys.self","elapsed","timestamp")
           out
         }
      } else { ##MONTEST
        test=speclist[i,"test"]
        stack=speclist[i,"stack"]
        treetype=speclist[i,"treetype"]
        est=try(montest(data=data,D="D",Z="Z",Y="Y",X=c("Xvar1","Xvar2","Xvar3"),stack=stack,treetype=treetype,test=test))
        if (!inherits(est,"try-error")) {
          data.frame(rep=r,speclist[i,],t(est$minp),t(est$time[1:3]),timestamp=Sys.time())
      } else {
        out=data.frame(rep=r,speclist[i,],t(rep(NA,8)),timestamp=Sys.time())
        colnames(out)=c("rep",colnames(speclist),"p.raw","p.holm","p.hochberg","p.BH","p.BY","user.self","sys.self","elapsed","timestamp")
        out
      }
     }
  }

if (restart==F|r>1) farbmachersims=rbind(data.frame(readRDS("farbmachersims.Rdata")),farbmachersims)
saveRDS(farbmachersims,file="farbmachersims.Rdata")
setTxtProgressBar(pb,r)
}
stopCluster(cl)

