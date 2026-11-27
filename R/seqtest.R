#' Tests of sequential-treatment IV conditions with two treatments
#'
#' \code{seqtest()} is a thin wrapper around \code{\link{montest}} for an instrument Z and two
#' (sequential) treatments D1 and D2, written \code{Y ~ X | FE | D1 + D2 ~ Z}. It parses the formula,
#' translates each requested condition into the appropriate \code{montest()} call and passes
#' all other arguments straight through. The conditions are
#' \describe{
#'   \item{\code{"KRD"}}{The Kwon-Roth conditions with D1 as the treatment and D2 as the only outcome,
#'     i.e. \code{montest(D2 ~ X | FE | D1 ~ Z, condition = "KR")}.}
#'   \item{\code{"KRDY"}}{The Kwon-Roth conditions with D1 as the treatment and (D2, Y) jointly as the
#'     outcome, i.e. \code{montest(D2 + Y ~ X | FE | D1 ~ Z, condition = "KR")}. The sets A range over the
#'     joint support of the binned outcomes; the \code{yval} labels in \code{$results} list the
#'     (D2, Y) tuples in each set and \code{$Ylookup} of the fit decodes them. Requires Y.}
#'   \item{\code{"KRDY2"}}{The intersection of two sets of Kwon-Roth conditions (use case 2): KR with D1 as the
#'     treatment and (D2, Y) jointly as the outcome (as \code{"KRDY"}), and KR with D2 as the treatment and Y as the
#'     outcome. The two problems are stacked as blocks of one \code{montest()} call using its \code{block} argument
#'     (cluster = the unit, or the user's \code{cluster}), so \code{pool = "block"} and \code{select = "block"}
#'     pool or select across them and everything is corrected as one family. Y is binned once (\code{Ysubsets},
#'     \code{gridtypeY}); D1 and D2 are binned by role (see below); \code{block} may not be passed. D1, D2 and Y may
#'     not have missing values; two-valued D1/D2 are recoded to 0/1. \code{$Wlookup} of the fit decodes the outcome
#'     codes appearing in \code{yval}; \code{block} is \code{"KRDY"} for the D1 problem and \code{"KRD2Y"} for the D2 problem. Requires Y and
#'     uses the forest search (\code{testtype = "forest"}).}
#'   \item{\code{"FSD"}}{The first stage difference condition
#'     \eqn{E[D1|Z=1]-E[D1|Z=0] \ge E[D2|Z=1]-E[D2|Z=0]}, tested as the simple first stage condition for the
#'     constructed treatment D1 - D2 (shifted by 1 and scored linearly), i.e.
#'     \code{montest(~ X | FE | D1 - D2 + 1 ~ Z, condition = "simple", linearD = TRUE)}. D2 need not be
#'     nested in D1. Rejected unless D1 and D2 are both binary.}
#'   \item{\code{"FSS"}}{The flipped first stage difference, i.e. the first stage sum
#'     \eqn{E[D1|Z=1]-E[D1|Z=0] + E[D2|Z=1]-E[D2|Z=0] \ge 0}, which is the first stage difference for D1 and
#'     1 - D2 (and equals the sum of the concordant-cell conditions (1,1) and (0,0) of the KR conditions with D2
#'     as outcome, whereas \code{"FSD"} is the sum of the discordant-cell conditions). It is tested as the simple first
#'     stage condition for the constructed treatment D1 + D2 (support \{0,1,2\}, scored linearly), i.e.
#'     \code{montest(~ X | FE | D1 + D2 ~ Z, condition = "simple", linearD = TRUE)}. It has power when D2 moves
#'     against D1 among compliers (e.g. a programme that is rarely repeated), where the complier slack is in the
#'     discordant cells. Together with \code{"FSD"}, the conditions imply \eqn{|\Delta E[D2]| \le \Delta E[D1]}.
#'     Rejected unless D1 and D2 are both binary.}
#'   \item{\code{"FSDY"}}{\code{"FSD"} localized in pre-specified bins of Y. With D1 as the treatment and the joint (D2, binned Y)
#'     as the outcome, as in \code{"KRDY"}, each bin b gives the summed discordant cell
#'     \eqn{1\{D1=1,D2=0,Y\in b\} - 1\{D1=0,D2=1,Y\in b\}}, i.e. the first stage of \eqn{(D1-D2)1\{Y\in b\}} is nonnegative. The bins
#'     are the \code{yval} margin, so \code{pool = "yval"} and \code{select = "yval"} pool or select across bins (pooling all bins gives
#'     \code{"FSD"}). Y is binned with \code{Ysubsets}/\code{gridtypeY}; \code{$Wlookup} decodes the (D2, Y-bin) codes. Requires Y;
#'     rejected unless D1 and D2 are both binary. Implemented by \code{montest()}'s \code{kr_sums} argument.}
#'   \item{\code{"FSSY"}}{\code{"FSS"} localized in bins of Y, as \code{"FSDY"} but with the summed concordant cells
#'     \eqn{1\{D1=1,D2=1,Y\in b\} - 1\{D1=0,D2=0,Y\in b\}}, i.e. the first stage of \eqn{(D1+D2-1)1\{Y\in b\}} is nonnegative.}
#'   \item{\code{"MWD"}}{The Mourifie-Wan conditions with D1 as the treatment and D2 as the outcome,
#'     i.e. \code{montest(D2 ~ X | D1 ~ Z, condition = "MW")}. As for \code{montest}'s \code{"MW"}, Z must be
#'     binary and \code{fml} may not contain fixed effects. D1 may be multivalued and ordered: the endpoint
#'     conditions are then pointwise sign restrictions, and the interior values are tested with outcome sets
#'     learned on the training half (see \code{montest}); the forest search it requires runs automatically.}
#'   \item{\code{"MWDY"}}{As \code{"MWD"} with (D2, Y) as outcomes, i.e.
#'     \code{montest(D2 + Y ~ X | D1 ~ Z, condition = "MW")}; each outcome is residualized separately (see
#'     \code{Y.res}) and all enter the forest. Requires Y.}
#'   \item{\code{"MWDY2"}}{The intersection of two sets of Mourifie-Wan conditions, mirroring \code{"KRDY2"}: MW with
#'     D1 as the treatment and (D2, Y) as outcomes (as \code{"MWDY"}), and MW with D2 as the treatment and Y as the
#'     outcome. The two problems are stacked as blocks of one \code{montest()} call (\code{block = }; cluster = the unit,
#'     or the user's \code{cluster}), so \code{pool = "block"} and \code{select = "block"} pool or select across them
#'     and everything is corrected as one family. The first outcome column is D2 in block 1 and a constant in block 2,
#'     the Y columns are shared; outcomes are residualized and nuisances fitted within each block, and Y is not binned.
#'     \code{block} and \code{Dsubsets} may not be passed. D1 and D2 may be multivalued and ordered (two-valued ones are
#'     recoded to 0/1; as for \code{"MWD"} the forest search it requires runs automatically), and rows with missing values in D1, D2 or Y
#'     are dropped with a message (for all conditions). With a multivalued D2 the second-treatment conditions are the Mourifie-Wan conditions for D2, which are
#'     pointwise only at the endpoints of D2 and learned-set budget bounds in the interior. The intersection of the two
#'     sets of conditions is not sharp when D1 is multivalued. \code{block} is \code{"MWDY"} for the D1 problem and \code{"MWD2Y"} for the D2 problem. Requires Y and the forest search.}
#' }
#'
#' \strong{Binning of D1, D2 and Y.} Pass D1, D2 and Y unbinned; seqtest bins each variable according to its role, with
#' the syntax of \code{montest}'s \code{Dcut}. D1 is always a treatment (\code{D1cut}). D2 is a treatment in the second block of
#' \code{"KRDY2"} and \code{"MWDY2"} (\code{D2Tcut}, by default as D1) and an outcome or conditioning variable W otherwise
#' (\code{D2cut}, by default as \code{D2Tcut} for the KR conditions and \code{"none"} for the MW conditions, whose outcomes need
#' no binning). Y is binned for the KR conditions (\code{Ycut}). Every default is \code{"equisized:4"}, which leaves a variable
#' with at most four values as it is. Binning W only changes the family of outcome events that are tested (power). Binning a
#' treatment is a legitimate choice, but it changes the exclusion restriction that is tested: Z may then affect W only through
#' the binned treatment. The blocks are separate testing problems, so D2 may be coded differently as a treatment than as W,
#' and differently from D1. \code{montest()} is told not to bin again.
#'
#' \strong{How the tests localize.} The Kwon-Roth conditions localize in the conditioning variables (D2, or
#' (D2, Y)) through the margin cells, which are tested with the sample split and selection of cells or globally (see \code{montest}); with covariates in X
#' the search additionally looks for violating subgroups in X within the cells (and can pool or select across them). The Mourifie-Wan conditions have no cells: the conditioning variables (D2, Y) and
#' X enter as features of the forest, so the search itself localizes jointly in (D2, X) or (D2, Y, X).
#'
#' @param fml A formula \code{Y ~ X | FE | D1 + D2 ~ Z}. Y may be omitted (\code{~ X | D1 + D2 ~ Z}) unless
#'   \code{"KRDY"} is requested, and may contain several variables joined by \code{+}.
#' @param data A \code{data.frame} or \code{data.table}.
#' @param condition Character vector, any of \code{"KRD"}, \code{"KRDY"}, \code{"KRDY2"}, \code{"FSD"}, \code{"FSS"},
#'   \code{"FSDY"}, \code{"FSSY"}, \code{"MWD"}, \code{"MWDY"}, \code{"MWDY2"} (or \code{"all"}, which fails if any of them is infeasible for the data).
#'   Defaults to \code{"KRD"}, plus \code{"KRDY"} if Y is given and \code{"FSD"} if D1 and D2 are both binary;
#'   \code{"FSS"}, \code{"FSDY"}, \code{"FSSY"} and the MW conditions are never run by default.
#' @param D1cut,D2cut,D2Tcut,Ycut Binning of D1 (always a treatment), of D2 as an outcome or conditioning variable W
#'   (\code{D2cut}) and as a treatment (\code{D2Tcut}, used in the second block of \code{"KRDY2"} and \code{"MWDY2"}), and of Y
#'   (KR conditions). Each is \code{NULL} (default as described under "Binning of D1, D2 and Y"), \code{"none"} (use the
#'   variable as it is), \code{"equisized:K"} or \code{"equidistant:K"} (K >= 2 bins; no change if the variable has at most K
#'   values), or a strictly increasing numeric vector of cut points giving left-open, right-closed bins, so \code{D1cut = 0}
#'   codes D1 as 0 | above 0 and \code{D2cut = c(0, 1)} codes D2 as 0 | 1 | 2 or more. Ignored by the FS conditions (and
#'   \code{Ycut} by the MW conditions).
#' @param ... Further arguments passed unchanged to \code{\link{montest}} (e.g. \code{stack = FALSE} to run
#'   the margins sequentially, or \code{Zcut}). The deprecated \code{Dsubsets}/\code{gridtypeD} and
#'   \code{Ysubsets}/\code{gridtypeY} are mapped to \code{D1cut} and \code{Ycut} with a warning; \code{Dcut} is not accepted.
#'
#' @return An object of class \code{"seqtest"}. With one condition, the full \code{montest} output for
#'   that condition plus \code{condition} and \code{call}. With several conditions, one full
#'   \code{montest} object per condition (named by condition, each with its own
#'   \code{minp}), plus a top-level \code{minp} (a single named vector) in which the test-sample
#'   p-values of all cells of all conditions are corrected as one family (Holm, Hochberg, BH, BY and Cauchy
#'   combination), \code{condition}, \code{call} and \code{version} (see \code{\link{montest_version}}). Pooling or selection across conditions
#'   (\code{pool}/\code{select = "condition"}) is not available; each condition is fit separately.
#'
#' @seealso montest
#' @export
seqtest <- function(fml, data, condition = NULL, D1cut = NULL, D2cut = NULL, D2Tcut = NULL, Ycut = NULL, ...) {
  mc <- match.call()
  dots <- list(...)
  if (any(c("condition", "fml", "linearD") %in% names(dots))) {
    stop("`fml`, `condition` and `linearD` are handled by seqtest() itself.", call. = FALSE)
  }
  if ("Dcut" %in% names(dots)) {
    stop("seqtest() has no `Dcut`: use `D1cut` for D1 and `D2cut`/`D2Tcut` for D2.", call. = FALSE)
  }
  ## Binning specs. The deprecated Dsubsets/gridtypeD and Ysubsets/gridtypeY are read from `...` (and
  ## removed from it): D1 and Y are binned here, never by montest(). Zcut/Zsubsets pass through to montest().
  spec_D1 <- resolve_bin_spec(D1cut, dots[["Dsubsets"]], if (is.null(dots[["gridtypeD"]])) "equisized" else dots[["gridtypeD"]],
                              !is.null(dots[["Dsubsets"]]), !is.null(dots[["gridtypeD"]]), "D1cut", "`Dsubsets`/`gridtypeD`")
  spec_Y <- resolve_bin_spec(Ycut, dots[["Ysubsets"]], if (is.null(dots[["gridtypeY"]])) "equisized" else dots[["gridtypeY"]],
                             !is.null(dots[["Ysubsets"]]), !is.null(dots[["gridtypeY"]]), "Ycut", "`Ysubsets`/`gridtypeY`")
  dots[c("Dsubsets", "gridtypeD", "Ysubsets", "gridtypeY")] <- NULL
  spec_W <- parse_bin_spec(D2cut, "D2cut")
  spec_T <- parse_bin_spec(D2Tcut, "D2Tcut")

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
  allowed <- c("KRD", "KRDY", "KRDY2", "FSD", "FSS", "FSDY", "FSSY", "MWD", "MWDY", "MWDY2")
  if (is.null(condition)) {
    condition <- c("KRD", if (!is.null(Y)) "KRDY", if (fsd_ok) "FSD")
  } else {
    condition <- match.arg(condition, c(allowed, "all"), several.ok = TRUE)
    if ("all" %in% condition) condition <- allowed
  }
  if ((!is.null(D1cut) || !is.null(D2cut) || !is.null(D2Tcut) || !is.null(Ycut)) && !any(grepl("^(KR|MW)", condition))) {
    message("`D1cut`, `D2cut`, `D2Tcut` and `Ycut` only apply to the KR and MW conditions; ignored.")
  }
  for (cn in intersect(c("KRDY2", "MWDY2"), condition)) {
    bad_dots <- intersect("block", names(dots))
    if (length(bad_dots)) {
      stop("Condition ", cn, " sets `", paste(bad_dots, collapse = "`, `"), "` itself.", call. = FALSE)
    }
    if (!is.null(Y)) {
      miss <- !stats::complete.cases(data[, c(D1, D2, Y), with = FALSE])
      if (any(miss)) {
        message("Condition ", cn, " needs the blocks to be row-aligned: dropping ", sum(miss),
                " rows with missing values in D1, D2 or Y (for all conditions).")
        data <- data[!miss]
        fsd_ok <- two_valued(data[[D1]]) && two_valued(data[[D2]])
      }
    }
  }
  for (cn in intersect(c("KRDY", "KRDY2", "FSDY", "FSSY", "MWDY", "MWDY2"), condition)) {
    if (is.null(Y)) {
      stop("Condition ", cn, " requires an outcome Y on the left hand side of `fml`.", call. = FALSE)
    }
  }
  for (cn in intersect(c("FSD", "FSS", "FSDY", "FSSY"), condition)) {
    if (!fsd_ok) {
      stop("Condition ", cn, " requires both D1 (", D1, ") and D2 (", D2, ") to be binary.", call. = FALSE)
    }
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

  ## Binning of D1, D2 and Y by role. D1 is always a treatment (D1cut). D2 is a treatment (blocks of
  ## KRDY2/MWDY2; D2Tcut, default as D1) and an outcome/conditioning variable W (KRD*, MWD*, block 1 of the
  ## *2 conditions; D2cut, default as D2Tcut for the KR conditions and "none" for the MW conditions, whose
  ## outcomes need no binning). Y is binned for the KR conditions (Ycut). The default is equisized:4.
  ## Everything is binned here and montest() is told not to bin again (Dcut = Ycut = "none").
  wcol <- function() if (is.null(dots$weight)) NULL else data[[dots$weight]]
  make_bins <- function(kr) {
    sT <- if (is.null(spec_T)) spec_D1 else spec_T
    sW <- if (!is.null(spec_W)) spec_W else if (kr) sT else list(type = "none")
    list(d1 = bin_vec(data[[D1]], spec_D1, wcol()),
         d2t = bin_vec(data[[D2]], sT, wcol()),
         d2w = bin_vec(data[[D2]], sW, wcol()))
  }

  fits <- list()
  for (cond in condition) {
    if (cond %in% c("KRD", "KRDY", "MWD", "MWDY")) {
      ## D2 (and Y for the *DY versions) as outcomes, D1 as the treatment
      f <- make_fml(c(D2, if (grepl("Y$", cond)) Y), as.name(D1))
      kr <- cond %in% c("KRD", "KRDY")
      bv <- make_bins(kr)
      dat_c <- data.table::copy(data)
      dots_c <- dots
      data.table::set(dat_c, j = D1, value = bv$d1)
      data.table::set(dat_c, j = D2, value = bv$d2w)
      ## D1 and D2 are binned above and, for the KR conditions, Y here, so montest() must not bin again.
      if (kr && grepl("Y$", cond)) {
        for (v in Y) data.table::set(dat_c, j = v, value = bin_vec(dat_c[[v]], spec_Y, wcol()))
      }
      dots_c$Dcut <- "none"
      if (kr) dots_c$Ycut <- "none"
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = dat_c,
                                              condition = substr(cond, 1L, 2L)), dots_c))
    } else if (cond == "KRDY2") {
      ## Two KR problems stacked as blocks of one montest() call:
      ##   block 1: treatment D1, outcome = joint code of (D2, binned Y)
      ##   block 2: treatment D2, outcome = binned Y
      ## Both copies of a unit share a cluster id, so sample splitting keeps them together.
      st_id <- "seq_id__"; st_blk <- "seq_block__"; st_T <- "seq_T__"; st_W <- "seq_W__"
      dat <- data.table::copy(data)
      data.table::set(dat, j = st_id, value = seq_len(nrow(dat)))
      recode01 <- function(x) if (two_valued(x)) as.integer(x == max(x, na.rm = TRUE)) else x
      ## D1 (treatment, block 1), D2 as a treatment (block 2) and D2 as W (block 1 cell label)
      ## are binned separately; see make_bins().
      bv <- make_bins(TRUE)
      d1v <- recode01(bv$d1); d2tv <- recode01(bv$d2t); d2wv <- recode01(bv$d2w)
      ybins <- paste0(Y, ".seqbin__")
      for (k in seq_along(Y)) data.table::set(dat, j = ybins[k], value = bin_vec(dat[[Y[k]]], spec_Y, wcol()))
      ylab <- do.call(paste, c(lapply(seq_along(Y), function(k) paste0(Y[k], "=", dat[[ybins[k]]])), sep = ","))
      ## Columns are assigned with set()/[[ ]] rather than inside `[`, where a user column
      ## called e.g. D1 or ylab would hide the local variable of that name.
      lab1 <- paste0("(", D2, "=", d2wv, ",", ylab, ")")
      lab2 <- paste0("(", ylab, ")")
      lab_levels <- c(unique(lab1), unique(lab2))
      b1 <- data.table::copy(dat)
      data.table::set(b1, j = st_blk, value = factor("KRDY", levels = c("KRDY", "KRD2Y")))
      data.table::set(b1, j = st_T, value = d1v)
      data.table::set(b1, j = st_W, value = match(lab1, lab_levels) - 1L)
      b2 <- data.table::copy(dat)
      data.table::set(b2, j = st_blk, value = factor("KRD2Y", levels = c("KRDY", "KRD2Y")))
      data.table::set(b2, j = st_T, value = d2tv)
      data.table::set(b2, j = st_W, value = match(lab2, lab_levels) - 1L)
      st <- data.table::rbindlist(list(b1, b2))
      wlookup <- data.table::data.table(code = seq_along(lab_levels) - 1L, label = lab_levels)
      args <- dots
      args$cluster <- if (is.null(dots$cluster)) st_id else dots$cluster
      f <- make_fml(st_W, as.name(st_T))
      ## Everything is binned above; the W code takes one value per (D2, Y-bin) cell, so "none" must allow them all.
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = st, condition = "KR", block = st_blk,
                                              Dcut = "none", Ycut = "none"), args))
      fits[[cond]]$Wlookup <- wlookup
    } else if (cond == "MWDY2") {
      ## Two MW problems stacked as blocks of one montest() call:
      ##   block 1: treatment D1, outcomes (D2, Y)
      ##   block 2: treatment D2, outcomes (0, Y)
      ## The first outcome column is D2 in block 1 and a constant in block 2 (D2 is the
      ## treatment there); the Y columns are shared. Both copies of a unit share a cluster id.
      st_id <- "seq_id__"; st_blk <- "seq_block__"; st_T <- "seq_T__"; st_O <- "seq_O__"
      dat <- data.table::copy(data)
      data.table::set(dat, j = st_id, value = seq_len(nrow(dat)))
      ## Two-valued treatments are recoded to 0/1; multivalued ones are used as they are unless
      ## binned (D1cut / D2Tcut / D2cut, see make_bins()).
      recode01 <- function(x) if (two_valued(x)) as.integer(x == max(x, na.rm = TRUE)) else x
      bv <- make_bins(FALSE)
      d1v <- recode01(bv$d1); d2tv <- recode01(bv$d2t); d2wv <- recode01(bv$d2w)
      b1 <- data.table::copy(dat)
      data.table::set(b1, j = st_blk, value = factor("MWDY", levels = c("MWDY", "MWD2Y")))
      data.table::set(b1, j = st_T, value = d1v)
      data.table::set(b1, j = st_O, value = as.numeric(d2wv))
      b2 <- data.table::copy(dat)
      data.table::set(b2, j = st_blk, value = factor("MWD2Y", levels = c("MWDY", "MWD2Y")))
      data.table::set(b2, j = st_T, value = d2tv)
      data.table::set(b2, j = st_O, value = 0)
      st <- data.table::rbindlist(list(b1, b2))
      args <- dots
      args$cluster <- if (is.null(dots$cluster)) st_id else dots$cluster
      f <- make_fml(c(st_O, Y), as.name(st_T))
      ## With `block`, D must not need binning (see montest()); the treatments are binned above.
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = st, condition = "MW", block = st_blk,
                                              Dcut = "none"), args))
    } else if (cond %in% c("FSDY", "FSSY")) {
      ## FSDY/FSSY: FSD/FSS within pre-specified Y-bins b. KR with treatment D1 and outcome W = (D2, binned Y), but
      ## with one summed cell per bin (montest's `kr_sums`):
      ##   FSDY: 1{D1=1, D2=0, Y in b} - 1{D1=0, D2=1, Y in b}  (= (d1 - d2) 1{Y in b}, the discordant cells)
      ##   FSSY: 1{D1=1, D2=1, Y in b} - 1{D1=0, D2=0, Y in b}  (= (d1 + d2 - 1) 1{Y in b}, the concordant cells)
      ## so that pool/select = "yval" work across bins; pooling all bins gives FSD/FSS.
      st_W <- "seq_W__"
      dat <- data.table::copy(data)
      recode01 <- function(x) as.integer(x == max(x, na.rm = TRUE))
      data.table::set(dat, j = D1, value = recode01(dat[[D1]]))
      data.table::set(dat, j = D2, value = recode01(dat[[D2]]))
      ybins <- paste0(Y, ".seqbin__")
      for (k in seq_along(Y)) data.table::set(dat, j = ybins[k], value = bin_vec(dat[[Y[k]]], spec_Y, wcol()))
      ylab <- do.call(paste, c(lapply(seq_along(Y), function(k) paste0(Y[k], "=", dat[[ybins[k]]])), sep = ","))
      lab <- paste0("(", D2, "=", dat[[D2]], ",", ylab, ")")
      lab_levels <- unique(lab)
      data.table::set(dat, j = st_W, value = match(lab, lab_levels) - 1L)
      wlookup <- data.table::data.table(code = seq_along(lab_levels) - 1L, label = lab_levels)
      ## code of the cell (D2 = d, bin b); -1 (matches nothing) if the cell is empty
      cell_code <- function(d, b) {
        m <- match(paste0("(", D2, "=", d, ",", b, ")"), lab_levels)
        ifelse(is.na(m), -1L, m - 1L)
      }
      bins <- sort(unique(ylab))
      d_hi <- if (cond == "FSDY") 0L else 1L
      sums <- data.table::data.table(label = paste0("(", bins, ")"),
                                     w_hi = cell_code(d_hi, bins), w_lo = cell_code(1L - d_hi, bins))
      args <- dots
      f <- make_fml(st_W, as.name(D1))
      fits[[cond]] <- do.call(montest, c(list(fml = f, data = dat, condition = "KR", kr_sums = sums,
                                              Ycut = "none"), args))
      fits[[cond]]$Wlookup <- wlookup
    } else {
      ## FSD: recode both to {0,1} (order preserving). D1 - D2 takes values in
      ## {-1,0,1} (D2 need not be nested in D1) and is scored linearly, so the
      ## first stage of D1 - D2 is exactly the difference in first stages. The
      ## +1 shift keeps the support at {0,1,2}: if a sample only contains two of
      ## the three values, montest's downgrade to a binary treatment then sees
      ## {0,1} instead of {-1,0}.
      d1 <- as.integer(data[[D1]] == max(data[[D1]], na.rm = TRUE))
      d2 <- as.integer(data[[D2]] == max(data[[D2]], na.rm = TRUE))
      ## FSS (flipped): D1 + D2, which already has support {0,1,2}; its first stage is the sum of
      ## the first stages, i.e. the first stage difference for D1 and 1 - D2.
      dd <- if (cond == "FSS") "Dseq_sum" else "Dseq_diff"
      while (dd %in% names(data)) dd <- paste0(dd, "_")
      dat <- data.table::copy(data)
      data.table::set(dat, j = dd, value = if (cond == "FSS") d1 + d2 else d1 - d2 + 1L)
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

  structure(c(fits, list(minp = minp, condition = condition, call = mc, version = montest_version())),
            class = "seqtest")
}

#' @export
print.seqtest <- function(x, ...) {
  if (length(x$condition) == 1L) {
    cat("seqtest, condition", x$condition, ": minimum p-values\n")
    print(signif(x$minp, 4))
    cat("\nFull montest output is stored in the object; components: ",
        paste(setdiff(names(x), c("condition", "call")), collapse = ", "), "\n", sep = "")
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
