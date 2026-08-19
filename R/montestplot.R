#' Plot montest results
#'
#' \code{montestplot()} visualizes the covariate imbalance of the testing subset selected
#' by \code{montest()}, the cutoff-selection curve that produced it, and, optionally, the
#' shares of the selected subset across pooled margin dimensions.
#'
#' The first panel is a bar plot of standardized differences in covariate means between
#' the testing subset and the full sample, using information stored in
#' \code{object$Xmeans}, \code{object$Xmeans_all}, and \code{object$XSD}. The second panel
#' (when \code{grid = TRUE}) plots the cutoff-search curve(s) (test statistic \code{t}
#' against candidate cutoff \code{tau}) from \code{object$grid}, with the actually-selected
#' cutoff(s) from \code{object$results}' train rows marked as vertical lines. If
#' margin-share information is available and requested, additional panels compare the
#' test-sample share with the full-sample share across the specified margin variables,
#' using \code{object$shares}.
#'
#' @param object An object returned by \code{montest()} for a single search method,
#'   expected to contain elements such as \code{results}, \code{Xmeans},
#'   \code{Xmeans_all}, \code{XSD}, \code{grid}, and optionally \code{shares}.
#' @param sample Optional sample index to plot. If \code{NULL}, the function attempts to
#'   infer the relevant sample from the object.
#' @param margins Optional character vector of margin names to visualize. Use
#'   \code{"all"} to plot all available stored margins.
#' @param grid Logical, default \code{TRUE}. If \code{TRUE}, adds a panel plotting the
#'   cutoff-selection curve(s) from \code{object$grid}. Requires \code{testtype = "forest"}
#'   results (see Details); set to \code{FALSE} to omit this panel.
#' @param numX Maximum number of covariates to display in the standardized-difference plot.
#'
#' @details
#' \code{Xmeans}, \code{Xmeans_all}, \code{XSD}, \code{shares}, \code{grid}, and
#' \code{results} are all keyed by whichever of \code{montest()}'s pool/select dimensions
#' (\code{zmargin}, \code{dval}, \code{yval}, \code{condition}, \code{equation},
#' \code{sample}) were \emph{not} pooled when the object was produced -- i.e. their row
#' structure depends on the \code{pool}/\code{select} choices made in the original
#' \code{montest()} call, and can range from a single pooled row to one row per
#' sample/margin combination. \code{montestplot()} recovers this key structure from
#' \code{object$results}' column names, locates the tested subset with the strongest
#' violation (smallest test-sample p-value, optionally restricted to \code{sample}), and
#' plots the matching row(s) rather than assuming a fixed layout.
#'
#' Covariates are ordered by the absolute magnitude of their standardized mean difference,
#' and only the top \code{numX} are plotted.
#'
#' The cutoff-selection panel (\code{grid = TRUE}) shows the train-sample search that
#' produced the selected subset. \code{object$grid}'s \code{sample} column identifies the
#' \emph{training} half that searched a given curve (the opposite of the half it was later
#' tested in). If the winning subset's test statistic was pooled across sample halves,
#' both training samples' curves and cutoffs are drawn together in one panel; otherwise
#' only the one training half whose cutoff was actually tested against the winning
#' (holdout) sample is drawn -- i.e. \code{sample} here, like elsewhere in this function,
#' refers to the holdout/test half, with the matching training half inferred as the
#' other one. \code{testtype = "CART"} results do not record \code{$grid} (or
#' \code{Xmeans}/\code{Xmeans_all}/\code{XSD}/\code{shares}); requesting \code{grid = TRUE}
#' (the default) on such an object raises an error rather than silently skipping the panel.
#'
#' When \code{margins} is specified, the function produces additional bar plots comparing
#' the share of the testing subset with the share of the full sample within each
#' requested, pooled margin category (from \code{object$shares}).
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
#' # Plot default summary: standardized differences + cutoff-selection curve
#' montestplot(out$forest)
#'
#' # Add all stored margin-share panels
#' montestplot(out$forest, margins = "all", numX = 8)
#'
#' # Skip the cutoff-selection panel
#' montestplot(out$forest, grid = FALSE)
#' }
#'
#' @seealso montest
#' @export

montestplot <- function(object, sample = NULL, margins = NULL, grid = TRUE, numX = 10) {

  ## montest()'s fixed vocabulary of pool/select dimensions. Xmeans/
  ## Xmeans_all/XSD/shares/grid/results are all keyed by whichever of these
  ## were *not* pooled -- this lets us recover that key structure purely
  ## from the column names actually present, instead of assuming a fixed
  ## layout.
  DIMS <- c("zmargin", "dval", "yval", "condition", "equation", "sample")

  ## `$grid` (and Xmeans/Xmeans_all/XSD/shares) are only ever recorded for
  ## testtype = "forest" results -- CART_test() never builds them. Since
  ## grid = TRUE is the default, a bare montestplot(cart_object) call would
  ## otherwise fail with a generic "nothing to plot" message; naming CART
  ## specifically here is far more actionable. `object$options$testtype`
  ## (added alongside `object$call` in montest()'s return value) makes this
  ## precise instead of guessing from `is.null(object$grid)` alone.
  if (isTRUE(grid)) {
    testtype_used <- object$options$testtype
    if (identical(testtype_used, "CART")) {
      stop(
        "montestplot()'s cutoff-selection panel (grid = TRUE, the default) requires ",
        "testtype = \"forest\" results; `object` was estimated with testtype = \"CART\", ",
        "which does not record $grid. Re-fit with testtype = \"forest\", or call ",
        "montestplot(object, grid = FALSE) to skip this panel.",
        call. = FALSE
      )
    }
    if (is.null(object$grid)) {
      stop(
        "montestplot()'s cutoff-selection panel (grid = TRUE, the default) requires ",
        "`object$grid`, which is not present on this object. Call ",
        "montestplot(object, grid = FALSE) to skip this panel.",
        call. = FALSE
      )
    }
  }

  has_grid  <- isTRUE(grid) && !is.null(object$grid)
  has_X     <- !is.null(object$Xmeans) && !is.null(object$Xmeans_all) && !is.null(object$XSD)
  has_shares <- !is.null(object$shares)

  if (!has_X && !has_shares && !has_grid) {
    stop(
      "`object` has no covariate-balance information (Xmeans/Xmeans_all/XSD), no ",
      "cutoff-selection grid, and no margin-share information (shares) to plot. Note ",
      "that testtype = \"CART\" results currently record none of these."
    )
  }

  res <- object$results
  if (is.null(res) || !all(c("train", "p.raw") %in% names(res))) {
    stop("`object$results` with `train`/`p.raw` columns is required to locate the tested subset.")
  }
  key_cols <- intersect(names(res), DIMS)

  ## ---- Locate the tested subset with the strongest violation ----
  target <- NULL
  if (length(key_cols) > 0L) {
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

  ## ---- Panel 1 (values): standardized differences in covariate means ----
  vals <- names_x <- NULL

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

  ## ---- Panel 2 (values): cutoff-selection curve(s) ----
  grid_panel <- NULL

  if (has_grid) {
    g  <- object$grid
    tr <- res[res$train == TRUE & res$relevant == 1L]

    cell_cols <- setdiff(key_cols, "sample")
    target_cell <- if (!is.null(target)) {
      target[, intersect(names(target), cell_cols), with = FALSE]
    } else {
      NULL
    }

    g_sub  <- g[match_rows(g, target_cell)]
    tr_sub <- tr[match_rows(tr, target_cell)]

    ## `sample` on `target` (when present) is the winning row's *holdout*
    ## sample -- NA means the winning cell's test statistic was pooled
    ## across sample halves. `g`/`tr`'s own `sample` column instead means
    ## the *training* half (each cutoff is always searched in one half and
    ## tested in the other), so a non-pooled holdout sample `s` maps to
    ## training sample `3 - s`.
    winner_sample <- if (!is.null(target) && "sample" %in% names(target)) target$sample[1L] else NA_integer_
    pooled_sample <- is.na(winner_sample)

    samples_to_plot <- if (pooled_sample) {
      sort(unique(tr_sub$sample))
    } else {
      3L - as.integer(winner_sample)
    }
    samples_to_plot <- samples_to_plot[!is.na(samples_to_plot)]

    gp <- g_sub[g_sub$sample %in% samples_to_plot]

    if (length(samples_to_plot) == 0L || nrow(gp) == 0L) {
      stop(
        "Could not locate matching cutoff-search data in `object$grid`/`object$results` ",
        "for the selected subset."
      )
    }

    main_title <- "Cutoff selection"
    if (!is.null(target_cell) && ncol(target_cell) > 0L) {
      main_title <- paste0(
        main_title, ": ",
        paste(paste(names(target_cell), as.character(unlist(target_cell)), sep = "="), collapse = ", ")
      )
    }

    grid_panel <- list(
      gp = gp, tr_sub = tr_sub, samples_to_plot = samples_to_plot, main_title = main_title
    )
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

  numplots <- (if (has_X) 1L else 0L) + (if (!is.null(grid_panel)) 1L else 0L) + length(margins)
  if (numplots > 1L) {
    op_outer <- par(mfrow = c(ceiling(numplots / 2), 2))
    on.exit(par(op_outer), add = TRUE)
  }

  ## ---- Panel 1 (draw): standardized differences in covariate means ----
  bar_mid <- NULL

  if (has_X) {
    ## find right margin
    theta <- 45 * pi / 180
    w <- strwidth(names_x, units = "inches", cex = 1)   # unrotated width
    t <- strheight(names_x, units = "inches", cex = 1)  # unrotated height

    ## Small physical gap between the axis and the rotated label's anchor,
    ## in inches (see below for why this must be a physical, not a data-
    ## unit, distance).
    gap_in <- 0.05

    # Vertical projection of a rotated rectangle (good approximation for text),
    # plus the anchor gap -- both physical distances below the axis.
    needed_in <- max(w * sin(theta) + t * cos(theta)) + gap_in

    # Convert inches to "lines" used by par(mar=...)
    line_in <- par("csi")  # character size (inches) per line
    bottom_lines <- ceiling(needed_in / line_in) + 1  # +1 line padding

    op_mar <- par(mar = c(bottom_lines, 4, 2, 1))
    on.exit(par(op_mar), add = TRUE)

    bar_mid <- barplot(vals, names.arg = names_x,
            main = "Standardized differences in means",
            ylab = "(test mean - full mean) / test SD",
            xaxt = "n")

    ## The label anchor needs to sit a small, fixed physical distance below
    ## the axis (matching the `gap_in` budgeted into the margin above). A
    ## plain `par("usr")[3] - 0.05` offset is in *data* units, so its
    ## physical size swings with the y-axis's scale -- for montest()'s
    ## standardized differences (typically well under 1 in magnitude), that
    ## alone can eat most or all of the reserved margin before any text is
    ## drawn, pushing (with adj = 1, srt = 45 extending the label down and
    ## to the left of its anchor) most of each label off the bottom of the
    ## device. Converting through inches keeps the gap scale-invariant.
    y_axis_in <- grconvertY(par("usr")[3], from = "user", to = "inches")
    y_anchor  <- grconvertY(y_axis_in - gap_in, from = "inches", to = "user")

    text(
      x      = bar_mid,
      y      = y_anchor,
      labels = names_x,
      srt    = 45,
      adj    = 1,
      xpd    = NA
    )
  }

  ## ---- Panel 2 (draw): cutoff-selection curve(s) ----
  if (!is.null(grid_panel)) {
    gp <- grid_panel$gp
    tr_sub <- grid_panel$tr_sub
    samples_to_plot <- grid_panel$samples_to_plot

    cols_pal <- c("black", "steelblue")

    plot(
      gp$tau, gp$t, type = "n",
      xlab = expression(tau), ylab = "t",
      main = grid_panel$main_title
    )

    leg_labels <- character(length(samples_to_plot))
    for (i in seq_along(samples_to_plot)) {
      s <- samples_to_plot[i]
      gs <- gp[gp$sample == s]
      data.table::setorder(gs, tau)
      lines(gs$tau, gs$t, col = cols_pal[((i - 1L) %% length(cols_pal)) + 1L])

      cut_s <- tr_sub[tr_sub$sample == s, tau_cutoff]
      if (length(cut_s) > 0L && is.finite(cut_s[1L])) {
        abline(v = cut_s[1L], col = cols_pal[((i - 1L) %% length(cols_pal)) + 1L], lty = 2)
      } else {
        message("No train-row cutoff found for sample ", s, "; curve shown without a cutoff line.")
      }
      leg_labels[i] <- paste("Sample", s)
    }

    if (length(samples_to_plot) > 1L) {
      legend("topright", legend = leg_labels,
             col = cols_pal[((seq_along(samples_to_plot) - 1L) %% length(cols_pal)) + 1L],
             lty = 1, bty = "n")
    }
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
