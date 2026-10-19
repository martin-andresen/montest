#' Chernozhukov-Lee-Rosen test of conditional moment inequalities (local linear)
#'
#' R port of the Stata command \code{clrtest, method("local")} (Chernozhukov, Kim, Lee and Rosen, 2015,
#' Stata Journal 15(1); based on Chernozhukov, Lee and Rosen, 2013, Econometrica). It tests
#' \eqn{H_0: E[y_j | x_j = x] \le 0} for all \eqn{x} in a grid and all inequalities \eqn{j = 1,\dots,J},
#' by estimating each conditional mean with a local linear regression, forming a uniform confidence
#' band with a multiplier simulation of the supremum of the studentized process, and rejecting if the
#' lower intersection bound \eqn{\max_x(\hat g(x) - \hat s(x) c_{1-\alpha})} is non-negative.
#' Because the standard errors are positive this is equivalent to \eqn{\max_x \hat g(x)/\hat s(x) \ge c}.
#'
#' The steps follow \code{clr_k.ado}: (1) plug-in bandwidth from a quartic polynomial fit
#' (\eqn{h = 2.036\,sd(x)\,(\sigma^2/\sum \hat g''^2 w)^{1/5}}, with weight \eqn{w} = 1 between the 10th and 90th
#' percentiles of \eqn{x}), undersmoothed by \eqn{n^{1/5} n^{-1/3.5}}; (2) biweight kernel
#' \eqn{(15/16)(1-u^2)^2}; (3) local linear estimate and pointwise standard error
#' \eqn{\hat s(z) = \sqrt{\overline{(K v)^2}/h}\,/\,(\hat f(z)\sqrt{nh})} with residuals \eqn{v}; (4) critical values from
#' the \eqn{(1-\alpha)}-quantile of the maximum over the grid (and inequalities) of the simulated
#' studentized process, using one matrix of standard normal draws shared by all inequalities and both
#' passes; (5) adaptive inequality selection (AIS): a first pass at level \eqn{\gamma_n = 1 - 0.1/\log n} keeps
#' the grid points with \eqn{\hat g \ge \max(\hat g - \hat s c_\gamma) - 2\hat s c_\gamma}, and the critical value is recomputed on
#' them with the same draws.
#'
#' The local linear method is local in \emph{one} conditioning variable per inequality (as in the Stata
#' command). The Stata random-number stream cannot be reproduced, so estimates, standard errors,
#' bandwidths and the selected grid agree with Stata, while critical values and p-values agree up to
#' simulation error.
#'
#' @param y Numeric vector or \eqn{n \times J} matrix; one column per inequality, \eqn{E[y_j|x_j] \le 0} under the null.
#' @param x Conditioning variable: a vector (used for every inequality) or an \eqn{n \times J} matrix.
#' @param grid Evaluation points: a numeric vector (used for every inequality) or a list of \eqn{J} vectors. Default:
#'   \code{seq} from the 2.5th to the 97.5th percentile of \code{x} in steps of \code{grid_step}
#'   (the grid used in the Farbmacher et al. MW replication code), or \code{ngrid} equally spaced points if \code{ngrid} is given.
#' @param ngrid,grid_step See \code{grid}.
#' @param alpha Significance level for \code{reject} and \code{cv}.
#' @param levels Optional confidence levels (e.g. \code{c(0.5, 0.9, 0.95, 0.99)}) for which the Stata outputs
#'   \code{bd} (lower intersection bound) and \code{det} (1 if the null is not rejected) are returned.
#' @param ais Apply adaptive inequality selection (Stata default).
#' @param undersmooth Undersmooth the bandwidth (Stata default).
#' @param bandwidth Optional bandwidth, scalar or one per inequality; default the plug-in rule above.
#' @param rnd Number of normal draws for the critical value (Stata default 10000).
#' @param seed Seed for the draws (Stata \code{seed(0)}); the previous RNG state is restored on exit. \code{NULL} for no seeding.
#' @return An object of class \code{clrtest}: \code{stat} (maximum studentized estimate over the selected grid), \code{cv},
#'   \code{pval} (share of simulated suprema at least \code{stat}), \code{reject}, \code{bounds}, \code{det}, \code{bandwidth},
#'   \code{n_grid} (points used), \code{n_dropped} (grid points removed because the standard error is zero or undefined,
#'   e.g. when the moment is identically zero in the window), \code{n_selected}, \code{theta} and \code{se} (per-inequality estimates on the full grid), \code{settings}.
#' @references Chernozhukov, V., S. Lee and A. M. Rosen (2013), Intersection bounds: estimation and inference, Econometrica 81(2), 667-737.
#'   Chernozhukov, V., W. Kim, S. Lee and A. M. Rosen (2015), Implementing intersection bounds in Stata, Stata Journal 15(1), 21-44.
#' @export
clrtest <- function(y, x, grid = NULL, ngrid = NULL, grid_step = 0.01, alpha = 0.05, levels = NULL,
                    ais = TRUE, undersmooth = TRUE, bandwidth = NULL, rnd = 10000L, seed = 0) {
  y <- as.matrix(y)
  if (!is.numeric(y)) stop("'y' must be numeric.")
  J <- ncol(y)
  x <- as.matrix(x)
  if (ncol(x) == 1L && J > 1L) x <- x[, rep(1L, J), drop = FALSE]
  if (nrow(x) != nrow(y) || ncol(x) != J) stop("'x' must be a vector or an n x J matrix conforming with 'y'.")
  keep <- stats::complete.cases(y, x)
  y <- y[keep, , drop = FALSE]
  x <- x[keep, , drop = FALSE]
  n <- nrow(y)
  if (n < 20L) stop("Too few complete observations.")
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
      on.exit(if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  ## ---- grids
  if (is.null(grid)) {
    grid <- lapply(seq_len(J), function(i) {
      q <- stats::quantile(x[, i], c(0.025, 0.975), names = FALSE)
      if (is.null(ngrid)) seq(q[1L], q[2L], by = grid_step) else seq(q[1L], q[2L], length.out = ngrid)
    })
  } else if (!is.list(grid)) {
    grid <- rep(list(as.numeric(grid)), J)
  }
  if (length(grid) != J) stop("'grid' must be a vector or a list with one vector per inequality.")
  grid <- lapply(grid, function(g) g[!is.na(g)])
  ng <- vapply(grid, length, 1L)
  if (any(ng == 0L)) stop("Empty grid.")

  ## ---- bandwidths
  h <- if (is.null(bandwidth)) {
    vapply(seq_len(J), function(i) .clr_bandwidth(y[, i], x[, i], undersmooth), 0)
  } else rep_len(as.numeric(bandwidth), J)
  if (any(!is.finite(h) | h <= 0)) stop("Bandwidth could not be computed; supply 'bandwidth'.")

  ## ---- local linear estimates, standard errors, and standardized simulated process
  est <- lapply(seq_len(J), function(i) .clr_estimate(y[, i], x[, i], grid[[i]], h[i]))
  ## grid points with no usable local variation (e.g. the moment is exactly zero for every unit in the window)
  ## give 0/0 in the studentization; no inference is possible there, so they are dropped and counted
  est <- lapply(est, function(e) {
    ok <- is.finite(e$ghat) & is.finite(e$se) & e$se > 0 & is.finite(e$nrm) & e$nrm > 0
    e$n_dropped <- sum(!ok)
    if (e$n_dropped > 0L) {
      e$ghat <- e$ghat[ok]; e$se <- e$se[ok]; e$Sv <- e$Sv[ok, , drop = FALSE]; e$nrm <- e$nrm[ok]
    }
    e
  })
  n_dropped <- vapply(est, `[[`, 0L, "n_dropped")
  ng <- ng - n_dropped
  if (any(ng == 0L)) stop("No grid point with positive standard error.")
  rnd <- as.integer(rnd)
  chunk <- max(1L, min(rnd, floor(2e7 / n)))
  W <- lapply(ng, function(m) matrix(0, m, rnd))       # standardized process, rows = grid points
  done <- 0L
  while (done < rnd) {
    rb <- min(chunk, rnd - done)
    xi <- matrix(stats::rnorm(n * rb), n, rb)
    for (i in seq_len(J)) W[[i]][, done + seq_len(rb)] <- est[[i]]$Sv %*% xi / est[[i]]$nrm
    done <- done + rb
  }
  supcol <- function(idx_list) {                         # max over selected rows of all inequalities, per draw
    mx <- rep(-Inf, rnd)
    for (i in seq_len(J)) if (length(idx_list[[i]]))
      mx <- pmax(mx, apply(W[[i]][idx_list[[i]], , drop = FALSE], 2L, max))
    mx
  }
  theta <- lapply(est, `[[`, "ghat")
  se <- lapply(est, `[[`, "se")

  ## ---- adaptive inequality selection
  all_idx <- lapply(ng, seq_len)
  if (ais) {
    gam <- 1 - 0.1 / log(n)
    cv_g <- stats::quantile(supcol(all_idx), gam, type = 1, names = FALSE)
    lb <- max(unlist(Map(function(t, s) t - s * cv_g, theta, se)))
    sel <- Map(function(t, s) which(t >= lb - 2 * s * cv_g), theta, se)
  } else sel <- all_idx
  nsel <- vapply(sel, length, 1L)
  if (sum(nsel) == 0L) stop("Adaptive inequality selection removed every grid point.")

  ## ---- test
  sup2 <- supcol(sel)
  stat <- max(unlist(Map(function(t, s, k) if (length(k)) (t / s)[k] else -Inf, theta, se, sel)))
  cvf <- function(L) stats::quantile(sup2, L, type = 1, names = FALSE)
  bound <- function(L) {
    cv <- cvf(L)
    max(unlist(Map(function(t, s, k) if (length(k)) (t - s * cv)[k] else -Inf, theta, se, sel)))
  }
  bounds <- det <- NULL
  if (!is.null(levels)) {
    bounds <- stats::setNames(vapply(levels, bound, 0), levels)
    det <- stats::setNames(as.integer(bounds < 0), levels)
  }
  structure(list(stat = stat, cv = cvf(1 - alpha), pval = mean(sup2 >= stat), reject = stat >= cvf(1 - alpha),
                 bounds = bounds, det = det, bandwidth = h, n_grid = ng, n_dropped = n_dropped, n_selected = nsel,
                 theta = theta, se = se,
                 settings = list(alpha = alpha, ais = ais, undersmooth = undersmooth, rnd = rnd, n = n, J = J)),
            class = "clrtest")
}

#' @export
print.clrtest <- function(x, ...) {
  cat("CLR test of conditional moment inequalities (local linear), ", x$settings$J, " inequalit",
      if (x$settings$J == 1L) "y" else "ies", ", n = ", x$settings$n, "\n", sep = "")
  cat(sprintf("  max studentized estimate = %.4f, critical value (%g%%) = %.4f, p-value = %.4f\n",
              x$stat, 100 * (1 - x$settings$alpha), x$cv, x$pval))
  cat("  grid points: ", paste(x$n_grid, collapse = ", "), "; selected by AIS: ", paste(x$n_selected, collapse = ", "),
      "; bandwidth: ", paste(signif(x$bandwidth, 4), collapse = ", "), "\n", sep = "")
  invisible(x)
}

## biweight kernel (Stata kl.mata)
.clr_kl <- function(u) {
  out <- (15 / 16) * (1 - u^2)^2
  out[abs(u) > 1] <- 0
  out
}

## plug-in bandwidth (Stata bsp: quartic polynomial fit)
.clr_bandwidth <- function(y, x, undersmooth) {
  n <- length(y)
  sdx <- stats::sd(x)
  xx <- (x - mean(x)) / sdx
  xreg <- cbind(1, xx, xx^2, xx^3, xx^4)
  para <- qr.solve(xreg, y)
  resi <- y - xreg %*% para
  sig2 <- mean(resi^2)
  drv2 <- 2 * para[3] + 6 * para[4] * xx + 12 * para[5] * xx^2
  q <- stats::quantile(xx, c(0.1, 0.9), names = FALSE)
  w0h <- as.numeric(xx >= q[1L] & xx <= q[2L])
  den <- sum(drv2^2 * w0h)
  num <- sig2 * (max(w0h) - min(w0h))
  h <- 2.036 * sdx * (num / den)^(1 / 5)
  if (undersmooth) h <- h * n^(1 / 5) * n^(-1 / 3.5)
  h
}

## local linear estimate at points z (Stata llre) given a precomputed kernel matrix K (|z| x n), D = x - z
.clr_llre_from <- function(K, D, y) {
  n <- ncol(K)
  W1 <- K * D
  S0 <- rowMeans(K)
  S1 <- rowMeans(W1)
  S2 <- rowMeans(W1 * D)
  T0 <- drop(K %*% y) / n
  T1 <- drop(W1 %*% y) / n
  den <- S2 * S0 - S1^2
  den <- den + 1e-8 * (den == 0)
  (S2 * T0 - S1 * T1) / den
}

## llre evaluated at all observations, in column chunks
.clr_llre_obs <- function(y, x, h) {
  n <- length(x)
  step <- max(1L, floor(5e6 / n))
  out <- numeric(n)
  for (s in seq(1L, n, by = step)) {
    idx <- s:min(n, s + step - 1L)
    D <- outer(x[idx], x, function(a, b) b - a)           # |idx| x n, rows = evaluation points
    out[idx] <- .clr_llre_from(.clr_kl(-D / h), D, y)
  }
  out
}

## estimates, standard errors and the weighted kernel matrix for one inequality
.clr_estimate <- function(y, x, z, h) {
  n <- length(y)
  D <- outer(z, x, function(a, b) b - a)                  # |z| x n
  K <- .clr_kl(D / h)
  ghat <- .clr_llre_from(K, D, y)
  fhat <- rowMeans(K) / h
  v <- y - .clr_llre_obs(y, x, h)
  Sv <- K * rep(v, each = length(z))                      # row z, column i: K(z,x_i) v_i
  nrm <- sqrt(rowSums(Sv^2))                              # = sqrt(n * mean((K v)^2))
  se <- sqrt(rowMeans(Sv^2) / h) / fhat / sqrt(n * h)
  list(ghat = ghat, se = se, Sv = Sv, nrm = nrm)
}
