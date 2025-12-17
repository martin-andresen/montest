montestplot<-function(object,sample=NULL,margins=NULL,numX=10) {
  if (is.null(sample)==TRUE) s=object$minsample
  if (is.null(margins)==FALSE) { ##Determine margins
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
    vals=vals[s,-1]
    setcolorder(vals,names(vals)[order(abs(unlist(vals[1])), decreasing = TRUE)])
    names=names(vals)
    vals=as.numeric(vals)
    if (length(vals)>numX)vals=vals[1:numX]
    barplot(vals, names.arg = names,
            main = "Standardized differences in means",
            ylab = "(test mean - full mean) / test SD")

  ##Plot margins
  if (is.null(margins)==FALSE) {
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
