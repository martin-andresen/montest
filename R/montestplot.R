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
#' with the full-sample share across the specified margin variables, using \code{object$shares}.
#'
#' @param object An object returned by \code{montest()} for a single search method,
#'   expected to contain elements such as \code{results}, \code{Xmeans},
#'   \code{Xmeans_all}, \code{XSD}, and optionally \code{shares}.
#' @param sample Optional sample index to plot. If \code{NULL}, the function attempts to
#'   infer the relevant sample from the object.
#' @param margins Optional character vector of margin names to visualize. Use
#'   \code{"all"} to plot all available stored margins.
#' @param numX Maximum number of covariates to display in the standardized-difference plot.
#'
#' @details
#' \code{Xmeans}, \code{Xmeans_all}, \code{XSD}, and \code{shares} are all keyed by
#' whichever of \code{montest()}'s pool/select dimensions (\code{zmargin}, \code{dval},
#' \code{yval}, \code{condition}, \code{equation}, \code{sample}) were \emph{not} pooled
#' when the object was produced -- i.e. their row structure depends on the \code{pool}/
#' \code{select} choices made in the original \code{montest()} call, and can range from a
#' single pooled row to one row per sample/margin combination. \code{montestplot()}
#' recovers this key structure from the column names actually present, locates the tested
#' subset with the strongest violation (smallest test-sample p-value in
#' \code{object$results}, optionally restricted to \code{sample}), and plots the matching
#' row(s) rather than assuming a fixed layout.
#'
#' Covariates are ordered by the absolute magnitude of their standardized mean difference,
#' and only the top \code{numX} are plotted. When \code{margins} is specified, the
#' function produces additional bar plots comparing the share of the testing subset with
#' the share of the full sample within each requested, pooled margin category (from
#' \code{object$shares}).
#'
#' \code{testtype = "CART"} results from \code{montest()} do not currently record
#' \code{Xmeans}/\code{Xmeans_all}/\code{XSD}/\code{shares}, so the corresponding panels
#' are unavailable for such objects.
#'
#' @return
#' Invisibly returns the bar midpoints from the standardized-difference plot, or
#' \code{NULL} if that panel could not be produced. The main purpose of the function is
#' side-effect plotting.
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

montestplot <- function(object, sample = NULL, margins = NULL, numX = 10) {

  ## montest()'s fixed vocabulary of pool/select dimensions. Xmeans/
  ## Xmeans_all/XSD/shares/results are all keyed by whichever of these were
  ## *not* pooled -- this lets us recover that key structure purely from the
  ## column names actually present, instead of assuming a fixed layout.
  DIMS <- c("zmargin", "dval", "yval", "condition", "equation", "sample")

  has_X <- !is.null(object$Xmeans) && !is.null(object$Xmeans_all) && !is.null(object$XSD)
  has_shares <- !is.null(object$shares)

  if (!has_X && !has_shares) {
    stop(
      "`object` has neither covariate-balance information (Xmeans/Xmeans_all/XSD) ",
      "nor margin-share information (shares) to plot. Note that testtype = \"CART\" ",
      "results currently record neither."
    )
  }

  key_cols <- intersect(names(if (has_X) object$Xmeans else object$shares), DIMS)

  ## ---- Locate the tested subset with the strongest violation ----
  target <- NULL
  if (length(key_cols) > 0L) {
    res <- object$results
    if (is.null(res) || !all(c("train", "p.raw") %in% names(res))) {
      stop("`object$results` with `train`/`p.raw` columns is required to locate the tested subset.")
    }
    tst <- res[res$train == FALSE & is.finite(res$p.raw)]
    if (!is.null(sample) && "sample" %in% names(tst)) {
      ## Copied to a differently-named variable because data.table's
      ## `[.data.table` NSE would otherwise resolve a bare `sample` symbol
      ## to the data.table's own `sample` column, not this function's
      ## `sample` argument, making the filter a silent no-op.
      sample_arg <- sample
      tst <- tst[tst$sample %in% sample_arg]
    }
    if (nrow(tst) == 0L) {
      stop("No finite test-sample p-values found for the requested `sample`.")
    }
    winner <- tst[which.min(tst$p.raw)]
    target <- winner[, intersect(key_cols, names(winner)), with = FALSE]
  }

  ## Row indices of `dt` whose key columns match `target` (matching only on
  ## columns present in both; NA matches NA).
  match_rows <- function(dt, target) {
    if (is.null(dt) || is.null(target) || ncol(target) == 0L) return(seq_len(NROW(dt)))
    cols <- intersect(names(target), names(dt))
    if (length(cols) == 0L) return(seq_len(NROW(dt)))
    ok <- rep(TRUE, NROW(dt))
    for (cl in cols) {
      tv <- target[[cl]][1L]
      dv <- dt[[cl]]
      ok <- ok & (if (is.na(tv)) is.na(dv) else (!is.na(dv) & dv == tv))
    }
    which(ok)
  }

  ## ---- Panel 1: standardized differences in covariate means ----
  bar_mid <- NULL

  if (has_X) {
    x_cols <- setdiff(names(object$Xmeans), key_cols)

    rows <- if (length(key_cols) == 0L) 1L else match_rows(object$Xmeans, target)

    if (length(rows) == 0L) {
      ## The winning cell's test statistic was pooled across sample halves
      ## (montest()'s adaptive pool = "sample" design), but Xmeans/XSD are
      ## always kept by-sample, so no single row matches exactly. Fall back
      ## to averaging the per-sample rows sharing the remaining keys.
      target_nosample <- target[, setdiff(names(target), "sample"), with = FALSE]
      rows <- match_rows(object$Xmeans, target_nosample)
      if (length(rows) == 0L) {
        stop("Could not locate the tested subset within Xmeans/Xmeans_all/XSD.")
      }
      message(
        "The selected subset's test statistic was pooled across sample halves; ",
        "showing covariate balance averaged across the matching sample-specific rows."
      )
    }

    Xm  <- colMeans(as.matrix(object$Xmeans[rows, x_cols, with = FALSE]))
    Xa  <- colMeans(as.matrix(object$Xmeans_all[rows, x_cols, with = FALSE]))
    Xsd <- colMeans(as.matrix(object$XSD[rows, x_cols, with = FALSE]))

    vals <- (Xm - Xa) / Xsd
    vals <- vals[order(abs(vals), decreasing = TRUE)]
    names_x <- names(vals)
    vals <- as.numeric(vals)

    if (length(vals) > numX) {
      vals <- vals[1:numX]
      names_x <- names_x[1:numX]
    }
  }

  ## ---- Resolve margin panels (test-sample vs. full-sample share) ----
  if (!is.null(margins)) {
    if (!has_shares) {
      message("No margin-share information (`object$shares`) is available; skipping margin panels.")
      margins <- NULL
    } else {
      marginnames <- setdiff(names(object$shares), c(key_cols, "share", "share_all"))
      if (length(marginnames) == 0L) {
        message("No pooled margins were recorded on `object`; skipping margin panels.")
        margins <- NULL
      } else {
        margins <- match.arg(margins, c(marginnames, "all"), several.ok = TRUE)
        if (identical(margins, "all")) margins <- marginnames
      }
    }
  }

  numplots <- (if (has_X) 1L else 0L) + length(margins)
  if (numplots > 1L) {
    op_outer <- par(mfrow = c(ceiling(numplots / 2), 2))
    on.exit(par(op_outer), add = TRUE)
  }

  if (has_X) {
    ## find right margin
    theta <- 45 * pi / 180
    w <- strwidth(names_x, units = "inches", cex = 1)   # unrotated width
    t <- strheight(names_x, units = "inches", cex = 1)  # unrotated height

    # Vertical projection of a rotated rectangle (good approximation for text)
    needed_in <- max(w * sin(theta) + t * cos(theta))

    # Convert inches to "lines" used by par(mar=...)
    line_in <- par("csi")  # character size (inches) per line
    bottom_lines <- ceiling(needed_in / line_in) + 1  # +1 line padding

    op_mar <- par(mar = c(bottom_lines, 4, 2, 1))
    on.exit(par(op_mar), add = TRUE)

    bar_mid <- barplot(vals, names.arg = names_x,
            main = "Standardized differences in means",
            ylab = "(test mean - full mean) / test SD",
            xaxt = "n")
    text(
      x      = bar_mid,
      y      = par("usr")[3] - 0.05,
      labels = names_x,
      srt    = 45,
      adj    = 1,
      xpd    = NA
    )
  }

  ## ---- Margin panels ----
  if (length(margins) > 0L) {
    shares_row <- object$shares
    match_cols <- intersect(key_cols, names(shares_row))
    if (length(match_cols) > 0L && !is.null(target)) {
      idx <- match_rows(shares_row, target[, intersect(names(target), match_cols), with = FALSE])
      if (length(idx) == 0L) {
        t2 <- target[, setdiff(intersect(names(target), match_cols), "sample"), with = FALSE]
        idx <- match_rows(shares_row, t2)
      }
      if (length(idx) > 0L) shares_row <- shares_row[idx]
    }

    for (m in margins) {
      vals_m <- shares_row[, lapply(.SD, sum), by = m, .SDcols = c("share", "share_all"), env = list(m = m)]
      mat <- t(as.matrix(vals_m[, c("share", "share_all")]))
      rownames(mat) <- c("Test", "Full")

      cols <- c("grey70", "steelblue")   # two series colors
      lbl <- vals_m[[m]]
      lbl[is.na(lbl)] <- "NA"
      barplot(
        mat,
        beside    = TRUE,
        names.arg = lbl,
        las       = 2,
        ylab      = "Share of sample",
        col       = cols,
        main      = paste("Margin:", m)
      )
      legend("topright",
             legend = rownames(mat),
             fill   = cols,
             title  = "Sample",
             bty    = "n")
    }
  }

  invisible(bar_mid)
}
