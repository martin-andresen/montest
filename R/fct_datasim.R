## Generalization of the DGP in Farbmacher et al. (2020)
## to multivalued treatments and instruments.
##
## Convention:
##   J = number of treatment thresholds, so D in {0, ..., J}
##   K = number of instrument margins, so Z in {0, ..., K}
##   mono_bad has columns:
##     j = treatment margin D >= j
##     k = instrument margin Z: k - 1 -> k
##  Note: alpha reparametrized in probability space.

fct_datasim <- function(
    setup, n,
    J = 1, K = 1,
    condition = NULL,
    mono_bad = NULL,
    excl_bad = NULL,
    alpha_good = 0.132,
    alpha_bad  = -0.273,
    gap = 1.0,
    gamma_bad = 1.25,
    tau = rep(1, J),
    eps = 1e-6
) {

  p <- 3
  betaXY <- c(0.3, 0.3, 0.3)

  if (length(tau) != J) {
    stop("tau must have length J")
  }

  cov <- matrix(c(1, 0.3,
                  0.3, 1), 2, 2)

  errors <- MASS::mvrnorm(n, rep(0, 2), cov)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("Xvar", 1:p)

  Xdf <- as.data.frame(X)

  eval_inside <- function(x, name) {
    if (is.character(x) && length(x) == 1L) {
      out <- eval(parse(text = x), envir = Xdf)
    } else {
      out <- x
    }

    if (!is.numeric(out)) {
      stop(name, " must be numeric or a string evaluating to numeric")
    }

    out
  }

  make_alpha_array <- function(x, name) {
    x <- eval_inside(x, name)

    if (length(x) == 1L) {
      arr <- array(x, dim = c(n, J, K))
    } else if (length(x) == n) {
      arr <- array(rep(x, each = J * K), dim = c(J, K, n))
      arr <- aperm(arr, c(3, 1, 2))
    } else if (length(x) == J * K) {
      mat <- matrix(x, nrow = J, ncol = K)
      arr <- array(rep(mat, each = n), dim = c(n, J, K))
    } else if (is.matrix(x) && all(dim(x) == c(J, K))) {
      arr <- array(rep(x, each = n), dim = c(n, J, K))
    } else if (length(x) == n * J * K) {
      arr <- array(x, dim = c(n, J, K))
    } else {
      stop(
        name, " must have length 1, n, J*K, or n*J*K, ",
        "or be a J by K matrix"
      )
    }

    if (any(!is.finite(arr))) {
      stop(name, " contains non-finite values")
    }

    arr
  }

  make_n_vector <- function(x, name) {
    x <- eval_inside(x, name)

    if (length(x) == 1L) {
      x <- rep(x, n)
    } else if (length(x) != n) {
      stop(name, " must have length 1 or n")
    }

    if (any(!is.finite(x))) {
      stop(name, " contains non-finite values")
    }

    as.numeric(x)
  }

  if (is.null(condition)) {
    viol <- rep(TRUE, n)
  } else {
    viol <- eval(parse(text = condition), envir = Xdf)
    if (!is.logical(viol) || length(viol) != n) {
      stop("condition must evaluate to a logical vector of length n")
    }
    viol[is.na(viol)] <- FALSE
  }

  if (setup == "A") {
    Z <- sample(0:K, n, replace = TRUE)
    b <- rep(0, n)

  } else if (setup == "B") {
    z_score <- 0.2 * rowSums(X) + rnorm(n)

    br <- unique(quantile(
      z_score,
      probs = seq(0, 1, length.out = K + 2),
      na.rm = TRUE
    ))

    Z <- as.integer(cut(
      z_score,
      breaks = br,
      include.lowest = TRUE,
      labels = FALSE
    )) - 1L

    b <- 0.5 * log1p(exp(rowSums(X)))

  } else {
    stop("invalid choice of setup")
  }

  base <- b + errors[, 1]

  alpha_good_arr <- make_alpha_array(alpha_good, "alpha_good")
  alpha_bad_arr  <- alpha_good_arr

  if (!is.null(mono_bad)) {
    if (!all(c("j", "k") %in% names(mono_bad))) {
      stop("mono_bad must contain columns j and k")
    }

    alpha_bad_eval <- make_alpha_array(alpha_bad, "alpha_bad")

    for (r in seq_len(nrow(mono_bad))) {
      j <- mono_bad$j[r]
      k <- mono_bad$k[r]

      if (j < 1 || j > J || k < 1 || k > K) {
        stop("mono_bad contains invalid j or k")
      }

      alpha_bad_arr[, j, k] <- alpha_bad_eval[, j, k]
    }
  }

  make_shift <- function(alpha_arr, base, gap, eps) {
    c_j <- (seq_len(J) - 1L) * gap
    shift <- array(0, dim = c(n, J, K + 1L))

    for (j in seq_len(J)) {
      p0 <- mean(base > c_j[j])

      pz <- matrix(NA_real_, nrow = n, ncol = K + 1L)
      pz[, 1L] <- p0

      for (k in seq_len(K)) {
        pz[, k + 1L] <- pz[, k] + alpha_arr[, j, k]
      }

      pz <- pmin(pmax(pz, eps), 1 - eps)

      q <- stats::quantile(
        base,
        probs = as.vector(1 - pz),
        type = 8,
        names = FALSE
      )

      shift[, j, ] <- c_j[j] - matrix(q, nrow = n, ncol = K + 1L)
    }

    shift
  }

  shift_good <- make_shift(alpha_good_arr, base, gap, eps)
  shift_bad  <- make_shift(alpha_bad_arr,  base, gap, eps)

  L <- matrix(NA_real_, n, J)

  for (j in seq_len(J)) {
    sh <- shift_good[cbind(seq_len(n), j, Z + 1L)]

    if (!is.null(mono_bad)) {
      sh[viol] <- shift_bad[cbind(which(viol), j, Z[viol] + 1L)]
    }

    raw_j <- base + sh - (j - 1L) * gap

    if (j == 1L) {
      L[, j] <- raw_j
    } else {
      L[, j] <- pmin(raw_j, L[, j - 1L] - 1e-8)
    }
  }

  D <- rowSums(L > 0)

  eta_good <- matrix(0, nrow = n, ncol = K)
  eta_bad  <- eta_good

  if (!is.null(excl_bad)) {
    if (any(excl_bad < 1 | excl_bad > K)) {
      stop("excl_bad contains invalid instrument margin k")
    }

    gamma_bad_vec <- make_n_vector(gamma_bad, "gamma_bad")

    for (k in excl_bad) {
      eta_bad[, k] <- gamma_bad_vec
    }
  }

  gamma_good_z <- cbind(0, t(apply(eta_good, 1L, cumsum)))
  gamma_bad_z  <- cbind(0, t(apply(eta_bad,  1L, cumsum)))

  gamma_z <- gamma_good_z[cbind(seq_len(n), Z + 1L)]

  if (!is.null(excl_bad)) {
    gamma_z[viol] <- gamma_bad_z[cbind(which(viol), Z[viol] + 1L)]
  }

  tau_D <- rep(0, n)
  for (j in seq_len(J)) {
    tau_D <- tau_D + tau[j] * as.numeric(D >= j)
  }

  Y <- as.vector(tau_D + gamma_z + X %*% betaXY + errors[, 2])

  data.frame(Y, D, Z, X)
}
