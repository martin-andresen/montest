#' Plot montest results
#'
#' \code{montestplot()} visualizes the covariate imbalance of the testing subset selected
#' by \code{montest()} and, optionally, the shares of the selected subset across pooled
#' margin dimensions.
#'
#' The first panel is a bar plot of standardized differences in covariate means between
#' the testing subset and the full sample, using information stored in
#' \code{object$Xmeans}, \code{object$Xmeans_all}, and \code{object$XSD}. If margin-share
#' information is available and requested, additional panels compare the test-sample share
#' with the full-sample share across the specified margin variables.
#'
#' @param object An object returned by \code{montest()} for a single search method,
#'   expected to contain elements such as \code{results}, \code{Xmeans},
#'   \code{Xmeans_all}, \code{XSD}, and optionally a table of margin shares.
#' @param sample Optional sample index to plot. If \code{NULL}, the function attempts to
#'   infer the relevant sample from the object.
#' @param margins Optional character vector of margin names to visualize. Use
#'   \code{"all"} to plot all available stored margins.
#' @param numX Maximum number of covariates to display in the standardized-difference plot.
#'
#' @details
#' Covariates are ordered by the absolute magnitude of their standardized mean difference,
#' and only the top \code{numX} are plotted. When \code{margins} is specified, the
#' function produces additional bar plots comparing the share of the testing subset with
#' the share of the full sample within each requested margin category.
#'
#' This function expects a fairly specific object structure. In the current source, the
#' plotting code appears to assume components named \code{Xmeans}, \code{Xmeans_all},
#' \code{XSD}, and a table of margin shares. If the returned \code{montest()} object uses
#' a different nesting convention, a wrapper method or small revision may be needed.
#'
#' @return
#' Invisibly returns the bar midpoints from the standardized-difference plot. The main
#' purpose of the function is side-effect plotting.
#'
#' @examples
#' \dontrun{
#' out <- montest(
#'   data = mydata,
#'   D = "D",
#'   Z = "Z",
#'   Y = "Y",
#'   X = c("x1", "x2", "x3"),
#'   test = "simple"
#' )
#'
#' # Plot default summary
#' montestplot(out$forest)
#'
#' # Plot standardized differences plus all stored margins
#' montestplot(out$forest, margins = "all", numX = 8)
#' }
#'
#' @seealso montest
#' @export

montestplot<-function(object,sample=NULL,margins=NULL,numX=10) {
  if (dim(object$results)[1]>2) {
    sno=1
    } else {
    if (is.null(sample)==TRUE) sno=object$minsample else sno=sample
  }
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
    vals=vals[sno,-1]
    setcolorder(vals,names(vals)[order(abs(unlist(vals[1])), decreasing = TRUE)])
    names=names(vals)
    vals=as.numeric(vals)
    if (length(vals)>numX) {
      vals=vals[1:numX]
      names=names[1:numX]
    }

    ##find right margin
    theta <- 45 * pi/180
    w <- strwidth(names, units = "inches", cex = 1)   # unrotated width
    t <- strheight(names, units = "inches", cex = 1)  # unrotated height

    # Vertical projection of a rotated rectangle (good approximation for text)
    needed_in <- max(w * sin(theta) + t * cos(theta))

    # Convert inches to "lines" used by par(mar=...)
    line_in <- par("csi")  # character size (inches) per line
    bottom_lines <- ceiling(needed_in / line_in) + 1  # +1 line padding

    op <- par(mar = c(bottom_lines, 4, 2, 1))

    ##plot
    barplot = barplot(vals, names.arg = names,
            main = "Standardized differences in means",
            ylab = "(test mean - full mean) / test SD",
            xaxt= "n")
    text(
      x      = barplot,
      y      = par("usr")[3] - 0.05,
      labels = names,
      srt    = 45,
      adj    = 1,
      xpd    = NA
    )
    barplot

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
