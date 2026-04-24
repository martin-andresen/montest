##This is a generalization of the DGP in Farbmacher et. al (2020) to multivalued treatments and instruments.
fct_datasim <- function(
    setup, n,
    J = 1, K = 1,
    condition = NULL,
    mono_bad = NULL,          # data.frame(k, m), k = D >= k, m = Z margin m-1 -> m
    excl_bad = NULL,          # vector of Z margins m
    alpha_good = 0.40,
    alpha_bad  = -0.75,
    gap = 1.0,
    gamma_bad = 1.25,
    tau = rep(1, K)           # tau[k] effect of crossing D >= k
) {

  p <- 3
  betaXY <- c(0.3, 0.3, 0.3)

  if (length(tau) != K) {
    stop("tau must have length K")
  }

  cov <- matrix(c(1, 0.3,
                  0.3, 1), 2, 2)

  errors <- MASS::mvrnorm(n, rep(0, 2), cov)

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("Xvar", 1:p)

  Xdf <- as.data.frame(X)

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
    Z <- sample(0:J, n, replace = TRUE)
    b <- rep(0, n)

  } else if (setup == "B") {
    z_score <- 0.2 * rowSums(X) + rnorm(n)

    br <- unique(quantile(
      z_score,
      probs = seq(0, 1, length.out = J + 2),
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

  alpha_good_mat <- matrix(alpha_good, nrow = K, ncol = J)
  alpha_bad_mat  <- alpha_good_mat

  if (!is.null(mono_bad)) {
    stopifnot(all(c("k", "m") %in% names(mono_bad)))

    for (r in seq_len(nrow(mono_bad))) {
      k <- mono_bad$k[r]
      m <- mono_bad$m[r]

      if (k < 1 || k > K || m < 1 || m > J) {
        stop("mono_bad contains invalid k or m")
      }

      alpha_bad_mat[k, m] <- alpha_bad
    }
  }

  shift_good <- cbind(0, t(apply(alpha_good_mat, 1, cumsum)))
  shift_bad  <- cbind(0, t(apply(alpha_bad_mat,  1, cumsum)))

  base <- b + errors[, 1]

  L <- matrix(NA_real_, n, K)

  for (k in seq_len(K)) {
    sh <- shift_good[k, Z + 1L]

    if (!is.null(mono_bad)) {
      sh[viol] <- shift_bad[k, Z[viol] + 1L]
    }

    raw_k <- base + sh - (k - 1) * gap

    if (k == 1L) {
      L[, k] <- raw_k
    } else {
      L[, k] <- pmin(raw_k, L[, k - 1L] - 1e-8)
    }
  }

  D <- rowSums(L > 0)

  eta_good <- rep(0, J)
  eta_bad  <- eta_good

  if (!is.null(excl_bad)) {
    if (any(excl_bad < 1 | excl_bad > J)) {
      stop("excl_bad contains invalid margin m")
    }

    eta_bad[excl_bad] <- gamma_bad
  }

  gamma_good_z <- c(0, cumsum(eta_good))
  gamma_bad_z  <- c(0, cumsum(eta_bad))

  gamma_z <- gamma_good_z[Z + 1L]

  if (!is.null(excl_bad)) {
    gamma_z[viol] <- gamma_bad_z[Z[viol] + 1L]
  }

  tau_D <- rep(0, n)
  for (k in seq_len(K)) {
    tau_D <- tau_D + tau[k] * as.numeric(D >= k)
  }

  Y <- as.vector(tau_D + gamma_z + X %*% betaXY + errors[, 2])

  data.frame(Y, D, Z, X)
}
