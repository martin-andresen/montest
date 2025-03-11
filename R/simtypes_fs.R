simtypes_fs=function(n=3000,J=1,K=1,L=3,condition=NULL,defier.probs=matrix(0,nrow=J,ncol=K),complier.probs=matrix(1/(J*K+J+1),nrow=J,ncol=K),drop=TRUE) {

    if (nrow(defier.probs)!=J|ncol(defier.probs)!=K) stop("argument defier.probs must take a matrix of size J x K")
    if (sum(defier.probs>1)>0|sum(defier.probs<0)>0|sum(defier.probs)>1) stop("matrix defier.probs must contain complier shares between 0 and 1 at each margin, and must sum to less than 1")

    if (nrow(complier.probs)!=J|ncol(complier.probs)!=K) stop("argument complier.probs must take a matrix of size J x K")
    if (sum(complier.probs>1)>0|sum(complier.probs<0)>0|sum(complier.probs)>1) stop("matrix complier.probs must contain complier shares between 0 and 1 at each margin")

  probs=matrix(NA,ncol=K+2,nrow=0)
  colnames(probs)=c(paste0("D",(0:K)),"p")
  for (j in 1:J) {
      for (k in 1:K) {
        if (complier.probs[j,k]>0) probs=rbind(probs,c(rep(j-1,k),rep(j,K-k+1),complier.probs[j,k]))
      }
    }
  if (sum(defier.probs)>0&is.null(condition)==T) {
    for (j in 1:J) {
      for (k in 1:K) {
        if (defier.probs[j,k]>0) probs=rbind(probs,c(rep(j,k),rep(j-1,K-k+1),defier.probs[j,k]))
      }
    }
  }

  if (L>0) {
    data=data.table(matrix(rnorm(n * L), n, L))
    colnames(data)[1:L]=paste0("Xvar",(1:L))
  }

  p=sum(probs[,"p"])
  if (p<1) {
    for (j in 0:J) {
      probs=rbind(probs,c(rep(j,K+1),(1-p)/(J+1)))
    }
  }

  if (is.null(condition)==F) {
    data[,v:=eval(parse(text = condition))]
    data=rbind(data[v==FALSE],data[v==TRUE])
    types=probs[sample(nrow(probs),size=sum(data[,v==FALSE]),replace=TRUE,prob=probs[,"p"]),paste0("D",0:K)]

    if (sum(defier.probs)>0) {
      probs=matrix(NA,ncol=K+2,nrow=0)
      colnames(probs)=c(paste0("D",(0:K)),"p")
      for (j in 1:J) {
        for (k in 1:K) {
          if (defier.probs[j,k]>0) probs=rbind(probs,c(rep(j,k),rep(j-1,K-k+1),defier.probs[j,k]))
        }
      }
    }

    p=sum(probs[,"p"])
    if (p<1) {
      for (j in 0:J) {
        probs=rbind(probs,c(rep(j,K+1),(1-p)/(J+1)))
      }
    }

    types=rbind(types,types=probs[sample(nrow(probs),size=sum(data[,v==TRUE]),replace=TRUE,prob=probs[,"p"]),paste0("D",0:K)])
    data=cbind(data,types)

    } else data=cbind(data,probs[sample(nrow(rbind(probs)),size=n,replace=TRUE,prob=probs[,"p"]),paste0("D",0:K)])


    data$Z=sample(0:K,replace=T,size=n)

    for (k in 0:K) {
      data[Z==k,D:=get(paste0("D",k))]
    }
  if (drop==TRUE) {
    data[,c(paste0("D",0:K))]=NULL
    if (is.null(condition)==FALSE) data=data[,v:=NULL]
  }
  return(as.data.frame(data))
}
