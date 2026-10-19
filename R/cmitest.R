#' Andrews-Shi test of conditional moment inequalities
#'
#' Tests \eqn{H_0: E[m_j(W) | X] \ge 0} almost surely for all \eqn{j = 1,\dots,k}, where \code{X} has one or
#' more columns (continuous and/or discrete; more than three conditioning variables is rarely practical,
#' because the number of boxes grows like \eqn{(2r)^{d_c}}). This follows Andrews and Shi (2013, Econometrica),
#' Section 9 (countable-cube instrument functions, asymptotic GMS critical value), and is intended as a
#' benchmark for the conditional-moment (Mourifie-Wan) tests in \code{\link{montest}}.
#'
#' The conditional inequalities are converted into the unconditional moments \eqn{E[m_j(W) g(X)] \ge 0}
#' for a collection of indicator functions \eqn{g} of "boxes". Continuous columns of \code{X} are
#' transformed to \eqn{(0,1)} (centered, rotated by the inverse square root of their covariance
#' matrix, then passed through the standard normal cdf). At scale \eqn{r = 1, \dots, r_1}
#' the continuous columns are cut into cubes of side \eqn{1/(2r)}, each interacted with every observed
#' value combination of the discrete columns. The number of boxes at scale \eqn{r} is
#' \eqn{L (2r)^{d_c}}, with \eqn{L} the number of discrete cells and \eqn{d_c} the number of continuous
#' columns. With only discrete columns there is a single scale (the value cells).
#'
#' For each box \eqn{g} the studentized moments are
#' \eqn{t_j(g) = \sqrt{n}\,\bar m_j(g)/\tilde\sigma_j(g)}, with
#' \eqn{\tilde\sigma_j^2(g) = \hat\Sigma_{jj}(g) + \varepsilon\,\hat\Sigma_{jj}(1)} (Andrews-Shi eq. 3.5).
#' The statistic is the Kolmogorov-Smirnov type \eqn{\max_g S(t(g))} or the Cramer-von Mises type
#' \eqn{\sum_r w(r) \, N_r^{-1} \sum_{g \in r} S(t(g))} with \eqn{w(r) = (r^2 + 100)^{-1}}, where \eqn{S} is
#' the Sum function \eqn{S_1 = \sum_j [t_j]_-^2} or the Max function \eqn{S_3 = \max_j [t_j]_-^2}.
#' The critical value is the \eqn{1-\alpha+\eta} quantile (plus \eqn{\eta}) of the statistic simulated from
#' the Gaussian process with the sample covariance kernel, shifted by the GMS slackness term
#' \eqn{\varphi_j(g) = B_n \mathbf{1}\{t_j(g)/\kappa_n > 1\}}.
#'
#' Not implemented: equality moments, the QLR function \eqn{S_2}, and the bootstrap critical value.
#' Weights across discrete cells in the CvM statistic are uniform (an extension; the paper only treats
#' continuous regressors with such interactions informally).
#'
#' @param m Numeric vector or \eqn{n \times k} matrix of moment values with \eqn{E[m|X] \ge 0} under the null.
#' @param X Vector, matrix or data frame with the conditioning variables (any number of columns; a warning is
#'   issued above three).
#' @param stat \code{"KS"} (default) or \code{"CvM"}.
#' @param S \code{"S1"} (sum, default) or \code{"S3"} (max); identical when \eqn{k = 1}.
#' @param discrete Logical vector, one entry per column of \code{X}. Default: a column is discrete if it has at
#'   most \code{maxlev} distinct values.
#' @param maxlev See \code{discrete}.
#' @param r1 Finest scale (cube side \eqn{1/(2 r_1)}). Default: the largest \eqn{r} such that the expected
#'   number of observations in the smallest cube is at least \code{min_obs}.
#' @param min_obs Target expected number of observations in the smallest cube (Andrews-Shi recommend 10-20).
#' @param maxboxes Upper bound on the total number of boxes; \code{r1} is lowered until it holds.
#' @param include_all If \code{TRUE}, add the whole-sample moment (the instrument \eqn{g \equiv 1}) as an extra box,
#'   counted as scale \eqn{r = 0} (CvM weight \eqn{(0^2+100)^{-1}}). Default \code{FALSE} (the paper's cube class
#'   contains only cubes of side \eqn{\le 1/2}). It helps against violations spread thinly over the whole sample.
#' @param eps Studentization constant \eqn{\varepsilon}.
#' @param kappa,B GMS constants; defaults \eqn{(0.3 \ln n)^{1/2}} and \eqn{(0.4 \ln n / \ln\ln n)^{1/2}}.
#' @param reps Number of simulation draws for the critical value.
#' @param alpha Nominal level.
#' @param eta Infinitesimal uniformity constant.
#' @param seed Optional seed (restores the previous RNG state on exit).
#' @return An object of class \code{cmitest}: \code{stat}, \code{cv}, \code{pval}, \code{reject}, \code{nboxes},
#'   \code{r1}, \code{boxes_per_scale}, \code{worst} (the most violating box: scale, cell, its studentized
#'   moments and a logical membership vector \code{idx}), and \code{settings}.
#' @references Andrews, D. W. K. and X. Shi (2013), Inference based on conditional moment inequalities,
#'   Econometrica 81(2), 609-666.
#' @export
cmitest <- function(m, X, stat = c("KS", "CvM"), S = c("S1", "S3"), discrete = NULL, maxlev = 10,
                    r1 = NULL, min_obs = 15, maxboxes = Inf, include_all = FALSE, eps = 0.05, kappa = NULL, B = NULL,
                    reps = 1000, alpha = 0.05, eta = 1e-6, seed = NULL) {
  stat <- match.arg(stat)
  S <- match.arg(S)
  m <- as.matrix(m)
  if (!is.numeric(m)) stop("'m' must be numeric.")
  n <- nrow(m)
  k <- ncol(m)
  X <- as.data.frame(X)
  if (nrow(X) != n) stop("'m' and 'X' must have the same number of rows.")
  dX <- ncol(X)
  if (dX < 1L) stop("'X' must have at least one column.")
  if (dX > 3L) warning("'X' has ", dX, " columns; the number of boxes grows exponentially in the number of continuous columns.")
  if (anyNA(m) || anyNA(X)) stop("Missing values in 'm' or 'X' are not allowed.")
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  if (is.null(discrete)) discrete <- vapply(X, function(x) length(unique(x)) <= maxlev, NA)
  if (length(discrete) != dX) stop("'discrete' must have one entry per column of 'X'.")
  if (is.null(kappa)) kappa <- sqrt(0.3 * log(n))
  if (is.null(B)) B <- sqrt(0.4 * log(n) / log(log(n)))

  ## ---- conditioning variables: discrete cell id and (0,1)-transformed continuous columns
  dc <- sum(!discrete)
  if (any(discrete)) {
    cellf <- interaction(X[discrete], drop = TRUE)
    cell <- as.integer(cellf)
    L <- nlevels(cellf)
    cell_labels <- levels(cellf)
  } else {
    cell <- rep(1L, n); L <- 1L; cell_labels <- "all"
  }
  if (dc > 0L) {
    Xc <- as.matrix(X[!discrete])
    mu <- colMeans(Xc)
    Xc <- sweep(Xc, 2L, mu)
    if (dc == 1L) {
      U <- stats::pnorm(Xc / stats::sd(Xc[, 1L]))
    } else {
      es <- eigen(stats::cov(Xc), symmetric = TRUE)
      if (min(es$values) <= 1e-12) stop("The continuous columns of 'X' are collinear.")
      Sih <- es$vectors %*% diag(1 / sqrt(es$values), dc) %*% t(es$vectors)
      U <- stats::pnorm(Xc %*% Sih)
    }
  }

  ## ---- scales and number of boxes
  if (dc == 0L) {
    scales <- NA_integer_
    nbox_r <- L
    if (L > maxboxes) warning("Number of discrete cells (", L, ") exceeds 'maxboxes'; maxboxes is ignored.")
  } else {
    if (is.null(r1)) {
      r1 <- floor(0.5 * (n / (L * min_obs))^(1 / dc))
      if (r1 < 1L) {
        warning("Sample too small for min_obs = ", min_obs, " in the coarsest cubes; using r1 = 1.")
        r1 <- 1L
      }
    }
    r1 <- as.integer(r1)
    if (r1 < 1L) stop("'r1' must be at least 1.")
    tot <- function(r) sum(L * (2 * seq_len(r))^dc)
    while (r1 > 1L && tot(r1) > maxboxes) r1 <- r1 - 1L
    if (tot(r1) > maxboxes)
      warning("Even r1 = 1 gives ", tot(1L), " boxes, more than 'maxboxes' = ", maxboxes, "; using r1 = 1.")
    scales <- seq_len(r1)
    nbox_r <- L * (2 * scales)^dc
  }
  sc <- scales                                         # scale codes (NA: discrete cells only; 0: whole sample)
  nb <- nbox_r
  wts <- if (dc == 0L) 1 else 1 / (scales^2 + 100)
  if (include_all) {
    sc <- c(0L, sc); nb <- c(1, nb); wts <- c(if (dc == 0L) 1 else 1 / 100, wts)
  }
  ns <- length(sc)

  ## cube id per scale (one cube per unit and scale)
  cube_id <- function(r) {
    if (isTRUE(r == 0)) return(rep(1, n))
    if (dc == 0L) return(cell)
    idx <- pmin(pmax(ceiling(U * 2 * r), 1), 2 * r)
    cub <- 1 + as.vector((idx - 1) %*% (2 * r)^(seq_len(dc) - 1L))    # mixed-radix cube index (double)
    (cell - 1) * (2 * r)^dc + cub
  }

  sig1 <- pmax(colMeans(m^2) - colMeans(m)^2, 0)
  if (any(sig1 <= 0)) warning("A moment has zero sample variance.")

  ## ---- observed statistic
  Sfun <- function(neg2) if (S == "S1") rowSums(neg2) else do.call(pmax, as.data.frame(neg2))
  obs <- vector("list", ns)
  Tn <- if (stat == "KS") -Inf else 0
  worst <- list(S = -Inf)
  for (si in seq_len(ns)) {
    r <- sc[si]
    cid0 <- cube_id(r)
    uids <- sort(unique(cid0))
    cid <- match(cid0, uids)                          # dense integer ids of the non-empty boxes
    s1 <- rowsum(m, cid, reorder = TRUE)
    s2 <- rowsum(m^2, cid, reorder = TRUE)
    mbar <- s1 / n
    sigt <- sqrt(pmax(s2 / n - mbar^2, 0) + sweep(matrix(eps, nrow(mbar), k), 2L, sig1, "*"))
    tt <- sqrt(n) * mbar / sigt
    neg2 <- pmax(-tt, 0)^2
    Sg <- Sfun(neg2)
    obs[[si]] <- list(cid = cid, ids = uids, mbar = mbar, sigt = sigt, tt = tt)
    if (stat == "KS") {
      Tn <- max(Tn, Sg)
    } else {
      Tn <- Tn + wts[si] * sum(Sg) / nb[si]
    }
    if (max(Sg) > worst$S) {
      jj <- which.max(Sg)
      worst <- list(S = max(Sg), scale = si, id = obs[[si]]$ids[jj], t = tt[jj, ])
    }
  }
  wr <- sc[worst$scale]
  widx <- cube_id(wr) == worst$id
  wcell_lab <- if (isTRUE(wr == 0)) "all" else if (dc == 0L) cell_labels[worst$id] else
    cell_labels[floor((worst$id - 1) / (2 * wr)^dc) + 1]
  worst_out <- list(scale = wr, cell = wcell_lab, t = worst$t,
                    n_in_box = sum(widx), idx = widx)

  ## ---- critical value: GMS, asymptotic (Gaussian process with the sample kernel)
  chunk <- max(1L, min(reps, floor(2e7 / n)))
  Tb <- numeric(0)
  done <- 0L
  while (done < reps) {
    rb <- min(chunk, reps - done)
    xi <- matrix(stats::rnorm(n * rb), n, rb)
    sxi <- colSums(xi)
    Tsim <- if (stat == "KS") rep(-Inf, rb) else numeric(rb)
    for (si in seq_len(ns)) {
      o <- obs[[si]]
      Ssim <- NULL
      for (j in seq_len(k)) {
        A <- rowsum(xi * m[, j], o$cid, reorder = TRUE)
        nu <- (A - outer(o$mbar[, j], sxi)) / sqrt(n) / o$sigt[, j]
        phi <- B * (o$tt[, j] / kappa > 1)
        neg2 <- pmax(-(nu + phi), 0)^2
        Ssim <- if (is.null(Ssim)) neg2 else if (S == "S1") Ssim + neg2 else pmax(Ssim, neg2)
      }
      if (stat == "KS") {
        Tsim <- pmax(Tsim, apply(Ssim, 2L, max))
      } else {
        Tsim <- Tsim + wts[si] * colSums(Ssim) / nb[si]
      }
    }
    Tb <- c(Tb, Tsim)
    done <- done + rb
  }
  cv <- unname(stats::quantile(Tb, min(1, 1 - alpha + eta), type = 1)) + eta
  pval <- mean(Tb >= Tn)

  structure(list(stat = Tn, cv = cv, pval = pval, reject = Tn > cv,
                 nboxes = sum(nb), r1 = if (dc == 0L) NA_integer_ else r1,
                 boxes_per_scale = stats::setNames(nb, c(if (include_all) "all", if (dc == 0L) "cells" else paste0("r", scales))),
                 worst = worst_out,
                 settings = list(stat = stat, S = S, discrete = discrete, min_obs = min_obs, maxboxes = maxboxes,
                                 include_all = include_all, eps = eps, kappa = kappa, B = B, reps = reps, alpha = alpha, n = n, k = k)),
            class = "cmitest")
}

#' @export
print.cmitest <- function(x, ...) {
  cat("Andrews-Shi conditional moment inequality test (", x$settings$stat, ", ", x$settings$S, ")\n", sep = "")
  cat(sprintf("  statistic = %.4f, critical value (%g%%) = %.4f, p-value = %.4f\n",
              x$stat, 100 * x$settings$alpha, x$cv, x$pval))
  cat(sprintf("  %d boxes", x$nboxes),
      if (!is.na(x$r1)) sprintf(", finest scale r1 = %d", x$r1), "\n", sep = "")
  cat(sprintf("  most violating box: cell '%s'%s, %d obs, studentized moment(s): %s\n",
              x$worst$cell,
              if (!is.na(x$worst$scale)) sprintf(", scale r = %d", x$worst$scale) else "",
              x$worst$n_in_box, paste(round(x$worst$t, 2), collapse = ", ")))
  invisible(x)
}
