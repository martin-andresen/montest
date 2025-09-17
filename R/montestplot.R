montestplot<-function(object,sample=NULL,margins=NULL) {
  if (is.null(sample)==TRUE) s=object$minsample
  if (is.null(margins)==FALSE) { ##Plot margins
    marginnames=colnames(object$margins)
    marginnames=marginnames[!marginnames %in% c("sample","share","share_all")]
    match.arg(margins,c(marginnames,"all"),several.ok=TRUE)
    if (margins=="all") margins=marginnames
    numplots=length(marginnames)+1
    op <- par(mfrow = c(ceiling(numplots/2), 2))
    on.exit(par(op), add = TRUE)
  }

  ## Plot X means
    vals=(object$Xmeans-object$Xmeans_all)/object$XSD
    names=colnames(object$Xmeans)[-1]
    vals=as.numeric(vals[s,-1])
    barplot(vals, names.arg = names, las = 2,
            main = "Standardized differences in means",
            ylab = "(test sample - full sample) / SD full sample")


  if (is.null(margins)==FALSE) { ##Plot margins
        data=object$margins
        data=data[sample==s]



        for (m in margins) {
          vals = data[,lapply(.SD,sum), by=m,.SDcols=c("share","share_all"),env=list(m=m)]
          mat <- t(as.matrix(vals[,c("share","share_all")]))
          rownames(mat) <- c("Test", "Full")

          cols <- c("grey70", "steelblue")   # two series colors
          lbl=vals[[m]]
          lbl[is.na(lbl)] <- "NA"
          barplot(
            mat,
            beside   = TRUE,
            names.arg= lbl,
            las      = 2,
            ylab     = "Share of sample",
            col      = cols,
            main   = paste("Margin:",m)
          )
          legend("topright",
                 legend = rownames(mat),     # c("v1","v2")
                 fill   = cols,
                 title = "Sample",
                 bty    = "n")
        }
        par(op)
  }


}
