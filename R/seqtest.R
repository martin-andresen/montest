#' Tests of sequential-treatment IV conditions with two treatments
#'
#' \code{seqtest()} is a thin wrapper around \code{\link{montest}} for an instrument Z and two
#' (sequential) treatments D1 and D2, written \code{Y ~ X | FE | D1 + D2 ~ Z}. It parses the formula,
#' translates each requested condition into the appropriate \code{montest()} call and passes
#' all other arguments straight through. The conditions are
#' \describe{
#'   \item{\code{"KRD"}}{The Kwan-Roth conditions with D1 as the treatment and D2 as the only outcome,
#'     i.e. \code{montest(D2 ~ X | FE | D1 ~ Z, condition = "KR")}.}
#'   \item{\code{"KRDY"}}{The Kwan-Roth conditions with D1 as the treatment and (D2, Y) jointly as the
#'     outcome, i.e. \code{montest(D2 + Y ~ X | FE | D1 ~ Z, condition = "KR")}. The sets A range over the
#'     joint support of the binned outcomes; the \code{yval} labels in \code{$results} list the
#'     (D2, Y) tuples in each set and \code{$Ylookup} of the fit decodes them. Requires Y.}
#'   \item{\code{"FSD"}}{The first stage difference condition
#'     \eqn{E[D1|Z=1]-E[D1|Z=0] \ge E[D2|Z=1]-E[D2|Z=0]}, tested as the simple first stage condition for the
#'     constructed treatment D1 - D2 (shifted by 1 and scored linearly), i.e.
#'     \code{montest(~ X | FE | D1 - D2 + 1 ~ Z, condition = "simple", linearD = TRUE)}. D2 need not be
#'     nested in D1. Rejected unless D1 and D2 are both binary.}
#'   \item{\code{"MWD"}}{The Mourifie-Wan conditions with D1 as the treatment and D2 as the outcome,
#'     i.e. \code{montest(D2 ~ X | D1 ~ Z, condition = "MW")}. As for \code{montest}'s \code{"MW"}, Z must be
#'     binary and \code{fml} may not contain fixed effects; D1 must be binary (an error otherwise).}
#'   \item{\code{"MWDY"}}{As \code{"MWD"} with (D2, Y) as outcomes, i.e.
#'     \code{montest(D2 + Y ~ X | D1 ~ Z, condition = "MW")}; each outcome is residualized separately (see
#'     \code{Y.res}) and all enter the forest. Requires Y.}
#' }
#'
#' @param fml A formula \code{Y ~ X | FE | D1 + D2 ~ Z}. Y may be omitted (\code{~ X | D1 + D2 ~ Z}) unless
#'   \code{"KRDY"} is requested, and may contain several variables joined by \code{+}.
#' @param data A \code{data.frame} or \code{data.table}.
#' @param condition Character vector, any of \code{"KRD"}, \code{"KRDY"}, \code{"FSD"}, \code{"MWD"},
#'   \code{"MWDY"} (or \code{"all"}, which fails if any of them is infeasible for the data).
#'   Defaults to \code{"KRD"}, plus \code{"KRDY"} if Y is given and \code{"FSD"} if D1 and D2 are both binary;
#'   the MW conditions are never run by default.
#' @param ... Further arguments passed unchanged to \code{\link{montest}} (e.g. \code{stack = FALSE} to run
#'   the margins sequentially).
#'
#' @return An object of class \code{"seqtest"}. With one condition, the full \code{montest} output for
#'   that condition plus \code{condition} and \code{call}. With several conditions, one full
#'   \code{montest} object per condition (named by condition, each with its own
#'   \code{minp}), plus a top-level \code{minp} (a single named vector) in which the test-sample
#'   p-values of all cells of all conditions are corrected as one family (Holm, Hochberg, BH, BY and Cauchy
#'   combination), \code{condition} and \code{call}. Pooling or adaptive selection across conditions
#'   (\code{pool}/\code{select = "condition"}) is not available; each condition is fit separately.
#'
#' @seealso montest
#' @export
seqtest <- function(fml, data, condition = NULL, ...) {
  mc <- match.call()
  dots <- list(...)
  if (any(c("condition", "fml", "linearD") %in% names(dots))) {
    stop("`fml`, `condition` and `linearD` are handled by seqtest() itself.", call. = FALSE)
  }

  data <- data.table::as.data.table(data.table::copy(data))

  ################ parse ################
  p <- parse_iv_formula(fml)
  if (is.null(p$iv)) {
    stop("No IV part found. Expected a formula like Y ~ X | FE | D1 + D2 ~ Z.", call. = FALSE)
  }
  if (length(p$rf_parts) == 0L) {
    stop("No reduced-form RHS (covariates) found before the IV part.", call. = FALSE)
  }
  d_terms <- attr(stats::terms(stats::as.formula(call("~", p$iv[[2L]]))), "term.labels")
  if (length(d_terms) != 2L || !all(d_terms %in% names(data))) {
    stop("seqtest() requires exactly two treatment variables in the IV part, ",
         "written D1 + D2 ~ Z, both present in `data`.", call. = FALSE)
  }
  D1 <- d_terms[1L]
  D2 <- d_terms[2L]
  z_expr <- p$iv[[3L]]
  Y <- if (p$has_lhs) all.vars(p$lhs) else NULL

  two_valued <- function(x) data.table::uniqueN(x, na.rm = TRUE) == 2L
  fsd_ok <- two_valued(data[[D1]]) && two_valued(data[[D2]])

  ################ conditions ################
  allowed <- c("KRD", "KRDY", "FSD", "MWD", "MWDY")
  if (is.null(condition)) {
    condition <- c("KRD", if (!is.null(Y)) "KRDY", if (fsd_ok) "FSD")
  } else {
    condition <- match.arg(condition, c(allowed, "all"), several.ok = TRUE)
    if ("all" %in% condition) condition <- allowed
  }
  for (cn in intersect(c("KRDY", "MWDY"), condition)) {
    if (is.null(Y)) {
      stop("Condition ", cn, " requires an outcome Y on the left hand side of `fml`.", call. = FALSE)
    }
  }
  if (any(c("MWD", "MWDY") %in% condition) && !two_valued(data[[D1]])) {
    stop("Conditions MWD and MWDY require a binary D1 (", D1, "), but it has ",
         data.table::uniqueN(data[[D1]], na.rm = TRUE), " distinct values.", call. = FALSE)
  }
  if ("FSD" %in% condition && !fsd_ok) {
    stop("Condition FSD requires both D1 (", D1, ") and D2 (", D2, ") to be binary.", call. = FALSE)
  }

  ################ build formulas ################
  env <- p$env
  make_fml <- function(lhs_vars, d_expr) {
    pipe <- rebuild_pipe(c(p$rf_parts, list(d_expr)))
    left <- if (length(lhs_vars)) {
      call("~", Reduce(function(a, b) call("+", a, b), lapply(lhs_vars, as.name)), pipe)
    } else {
      call("~", pipe)
    }
    stats::as.formula(call("~", left, z_expr), env = env)
  }

  fits <- list()
  for (cond in condition) {
    if (cond %in% c("KRD", "KRDY", "MWD", "MWDY")) {
      ## D2 (and Y for the *DY versions) as outcomes, D1 as the treatment
      f <- make_fml(c(D2, if (grepl("Y$", cond)) Y), as.name(D1))
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = data,
                                              condition = substr(cond, 1L, 2L)), dots))
    } else {
      ## FSD: recode both to {0,1} (order preserving). D1 - D2 takes values in
      ## {-1,0,1} (D2 need not be nested in D1) and is scored linearly, so the
      ## first stage of D1 - D2 is exactly the difference in first stages. The
      ## +1 shift keeps the support at {0,1,2}: if a sample only contains two of
      ## the three values, montest's downgrade to a binary treatment then sees
      ## {0,1} instead of {-1,0}.
      d1 <- as.integer(data[[D1]] == max(data[[D1]], na.rm = TRUE))
      d2 <- as.integer(data[[D2]] == max(data[[D2]], na.rm = TRUE))
      dd <- "Dseq_diff"
      while (dd %in% names(data)) dd <- paste0(dd, "_")
      dat <- data.table::copy(data)
      dat[, (dd) := d1 - d2 + 1L]
      f <- make_fml(NULL, as.name(dd))
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = dat, condition = "simple", linearD = TRUE), dots))
    }
  }

  ## One condition: the montest object itself, plus the call.
  if (length(fits) == 1L) {
    fit <- fits[[1L]]
    names(fit)[names(fit) == "call"] <- "montest_call"
    return(structure(c(fit, list(condition = condition, call = mc)), class = "seqtest"))
  }

  ## Several conditions: one full montest object per condition, and a top-level
  ## minp from correcting the test-sample p-values of ALL cells of ALL conditions
  ## as one family (same recipe as montest's own minp). Holm, Hochberg, BH and BY
  ## are valid under the dependence between conditions (same data), as is CCT.
  praw <- unlist(lapply(fits, function(x) x$results[train == FALSE, p.raw]), use.names = FALSE)
  praw <- replace(praw, is.na(praw), 1)
  minp <- c(
    p.raw = min(praw),
    vapply(c(holm = "holm", hochberg = "hochberg", BH = "BH", BY = "BY"),
           function(m) min(stats::p.adjust(praw, method = m)), numeric(1L)),
    p.CCT = cct_pvalue(praw)
  )
  names(minp)[2:5] <- paste0("p.", c("holm", "hochberg", "BH", "BY"))

  structure(c(fits, list(minp = minp, condition = condition, call = mc)), class = "seqtest")
}

#' @export
print.seqtest <- function(x, ...) {
  if (length(x$condition) == 1L) {
    cat("seqtest, condition", x$condition, ": minimum p-values\n")
    print(signif(x$minp, 4))
  } else {
    cat("seqtest, conditions", paste(x$condition, collapse = ", "),
        ": minimum p-values corrected across all cells of all conditions\n")
    print(signif(x$minp, 4))
    for (cn in x$condition) {
      cat("\n", cn, ":\n", sep = "")
      print(signif(x[[cn]]$minp, 4))
    }
  }
  invisible(x)
}
