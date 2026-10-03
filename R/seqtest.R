#' Tests of sequential-treatment IV conditions with two treatments
#'
#' \code{seqtest()} is a thin wrapper around \code{\link{montest}} for an instrument Z and two
#' (sequential) treatments D1 and D2, written \code{Y ~ X | FE | D1 + D2 ~ Z}. It parses the formula,
#' translates each requested condition into the appropriate \code{montest()} call and passes
#' all other arguments straight through. The three conditions are
#' \describe{
#'   \item{\code{"KRD"}}{The Kwan-Roth conditions with D1 as the treatment and D2 as the only outcome,
#'     i.e. \code{montest(D2 ~ X | FE | D1 ~ Z, condition = "KR")}.}
#'   \item{\code{"KRDY"}}{The Kwan-Roth conditions with D1 as the treatment and (D2, Y) jointly as the
#'     outcome, i.e. \code{montest(D2 + Y ~ X | FE | D1 ~ Z, condition = "KR")}. The sets A range over the
#'     joint support of the binned outcomes; the \code{yval} labels in \code{$results} list the
#'     (D2, Y) tuples in each set and \code{$Ylookup} of the fit decodes them. Requires Y.}
#'   \item{\code{"FSD"}}{The first stage difference condition
#'     \eqn{E[D1|Z=1]-E[D1|Z=0] \ge E[D2|Z=1]-E[D2|Z=0]}, tested as the simple first stage condition for the
#'     constructed treatment D1 - D2, i.e. \code{montest(~ X | FE | D1 - D2 ~ Z, condition = "simple")}.
#'     Rejected unless D1 and D2 are both binary (and D2 <= D1, so that D1 - D2 is a binary treatment).}
#' }
#'
#' @param fml A formula \code{Y ~ X | FE | D1 + D2 ~ Z}. Y may be omitted (\code{~ X | D1 + D2 ~ Z}) unless
#'   \code{"KRDY"} is requested, and may contain several variables joined by \code{+}.
#' @param data A \code{data.frame} or \code{data.table}.
#' @param condition Character vector, any of \code{"KRD"}, \code{"KRDY"}, \code{"FSD"} (or \code{"all"}).
#'   Defaults to \code{"KRD"}, plus \code{"KRDY"} if Y is given and \code{"FSD"} if D1 and D2 are both binary.
#' @param ... Further arguments passed unchanged to \code{\link{montest}}.
#'
#' @return An object of class \code{"seqtest"}: a list with \code{fits} (one \code{montest} result per
#'   condition, named by condition), \code{minp} (a matrix of minimum p-values, one row per condition)
#'   and \code{call}. No correction is made across conditions.
#'
#' @seealso montest
#' @export
seqtest <- function(fml, data, condition = NULL, ...) {
  mc <- match.call()
  dots <- list(...)
  if (any(c("condition", "fml") %in% names(dots))) {
    stop("`fml` and `condition` are handled by seqtest() itself.", call. = FALSE)
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
  allowed <- c("KRD", "KRDY", "FSD")
  if (is.null(condition)) {
    condition <- c("KRD", if (!is.null(Y)) "KRDY", if (fsd_ok) "FSD")
  } else {
    condition <- match.arg(condition, c(allowed, "all"), several.ok = TRUE)
    if ("all" %in% condition) condition <- allowed
  }
  if ("KRDY" %in% condition && is.null(Y)) {
    stop("Condition KRDY requires an outcome Y on the left hand side of `fml`.", call. = FALSE)
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
    if (cond == "KRD") {
      f <- make_fml(D2, as.name(D1))
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = data, condition = "KR"), dots))
    } else if (cond == "KRDY") {
      f <- make_fml(c(D2, Y), as.name(D1))
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = data, condition = "KR"), dots))
    } else {
      ## FSD: recode both to {0,1} (order preserving), require nesting D2 <= D1
      d1 <- as.integer(data[[D1]] == max(data[[D1]], na.rm = TRUE))
      d2 <- as.integer(data[[D2]] == max(data[[D2]], na.rm = TRUE))
      if (any(d2 > d1, na.rm = TRUE)) {
        stop("Condition FSD requires D2 <= D1 for everyone (D2 = 1 implies D1 = 1), ",
             "so that D1 - D2 is a binary treatment.", call. = FALSE)
      }
      dd <- "Dseq_diff"
      while (dd %in% names(data)) dd <- paste0(dd, "_")
      dat <- data.table::copy(data)
      dat[, (dd) := d1 - d2]
      f <- make_fml(NULL, as.name(dd))
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = dat, condition = "simple"), dots))
    }
  }

  minp <- do.call(rbind, lapply(fits, function(x) x$minp))
  structure(list(fits = fits, minp = minp, call = mc), class = "seqtest")
}

#' @export
print.seqtest <- function(x, ...) {
  cat("seqtest: minimum p-values by condition (no correction across conditions)\n")
  print(signif(x$minp, 4))
  invisible(x)
}
