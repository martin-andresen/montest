## Generalization of the DGP in Farbmacher et al. (2020)
## to multivalued treatments and instruments.
##
## Convention:
##   J = number of treatment thresholds, so D in {0, ..., J}
##   K = number of instrument margins, so Z in {0, ..., Kop}
##   mono_bad has columns:
##     j = treatment margin D >= j
##     k = instrument margin Z: k - 1 -> k
##  Note: alpha reparametrized in probability space.
##  nFE1/nFE2: optional number of levels for one/two grouping (fixed-effect)
##    variables FE1/FE2. When set, Z's distribution shifts by FE group
##    (strength fe_z_strength) and Y gets an additive FE main effect (sd
##    fe_y_sd) -- an FE-correlated instrument plus an FE-correlated outcome
##    confound, i.e. exactly the setting where absorbing FE is load-bearing
##    for validity. NULL (default) reproduces the previous FE-free behavior.
##  fe_d_strength/fe_zd_corr: setup "A". fe_d_strength is the sd of a
##    group-level shift added to D's latent index; fe_zd_corr is its
##    correlation with the group's own Z-shift. It is recycled to length 2
##    -- element 1 drives FE1, element 2 drives FE2 -- so a single FE
##    dimension can be made adversarial while the other stays benign (e.g.
##    fe_d_strength = c(0.6, 0) hits FE1 only). With fe_zd_corr < 0, FE
##    groups where Z = 1 is more common get a systematically LOWER baseline
##    P(D = 1): the first stage stays strictly positive within every group,
##    but the FE-marginal first stage E[D|Z=1] - E[D|Z=0] can turn negative
##    (a Simpson-type artifact that absorbing that FE removes).
##    fe_d_strength = 0 (default) = no D-side FE shift, i.e. previous
##    behavior.
##  z_score_expr: optional string overriding setup "B"'s Z-generating
##    formula, evaluated the same way alpha_good is -- lets a caller plant
##    genuine nonlinearity in Z's true dependence on X (e.g. to test
##    parametric=TRUE under a misspecified fml.Z). For K == 1 (binary Z),
##    it is interpreted as P(Z=1|X) directly (clipped to [eps,1-eps], drawn
##    via rbinom) -- the same probability-space convention alpha_good
##    already uses for D -- rather than a latent index perturbed with noise
##    and quantile-cut, so a linear-probability-model fit of the matching
##    functional form can recover it exactly. For K > 1 it still falls back
##    to the original noisy-index/quantile-cut mechanism (0.2*rowSums(X) by
##    default), since a single probability doesn't generalize to more than
##    two categories.
fct_datasim <- function(
    setup, n,
    J = 1, K = 1,
    condition = NULL,
    mono_bad = NULL,
    excl_bad = NULL,
    alpha_good = 0.132,
    alpha_bad  = -0.273,
    gamma_bad = 1.25,
    tau = rep(1, J),
    eps = 1e-6,
    min_cell_prob = 0.02,
    nFE1 = NULL,
    nFE2 = NULL,
    fe1_balanced = FALSE,
    fe2_balanced = FALSE,
    fe_z_strength = 0.75,
    fe_y_sd = 1,
    fe_d_strength = 0,
    fe_zd_corr = 0,
    z_score_expr = NULL,
    return_design = FALSE
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

  has_FE1 <- !is.null(nFE1)
  has_FE2 <- !is.null(nFE2)

  ## Group sizes follow a mild power-law (rather than a uniform draw) so
  ## that, with a moderate number of groups, a realistic mix of large and
  ## singleton/near-singleton groups occurs -- giving montest's
  ## drop_singletons/drop_novar_Z options something real to do, rather than
  ## requiring an unrealistically large nFE1/nFE2 to ever see a singleton.
  ## `fe1_balanced`/`fe2_balanced = TRUE` overrides this with a uniform draw
  ## (roughly equal group sizes) -- appropriate when the FE stands in for an
  ## ordinal grouping like cohort or year, where levels are of similar size.
  FE1 <- if (has_FE1) {
    pr <- if (isTRUE(fe1_balanced)) NULL else (seq_len(nFE1))^(-1)
    sample.int(nFE1, n, replace = TRUE, prob = pr)
  } else {
    NULL
  }

  FE2 <- if (has_FE2) {
    pr <- if (isTRUE(fe2_balanced)) NULL else (seq_len(nFE2))^(-1)
    sample.int(nFE2, n, replace = TRUE, prob = pr)
  } else {
    NULL
  }

  if (has_FE1) Xdf$FE1 <- FE1
  if (has_FE2) Xdf$FE2 <- FE2

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
    if (has_FE1) {
      ## Z's marginal distribution shifts by FE group (and, if present, a
      ## second FE group) -- unconditionally correlated with FE, but still
      ## exogenous *within* FE, since nothing else about the DGP depends on
      ## the FE-group shift used here. This is what makes controlling for
      ## FE load-bearing for validity: an instrument "as good as randomly
      ## assigned" only within group, the standard motivation for a
      ## fixed-effects IV design (e.g. judge/region/cohort designs).
      fe1_z_shift <- stats::rnorm(nFE1, 0, fe_z_strength)
      z_score <- fe1_z_shift[FE1]

      if (has_FE2) {
        fe2_z_shift <- stats::rnorm(nFE2, 0, fe_z_strength)
        z_score <- z_score + fe2_z_shift[FE2]
      }

      z_score <- z_score + rnorm(n)

      ## Group-level shift in D's latent index, correlated (via fe_zd_corr)
      ## with the group's own Z-shift -- see the fe_d_strength / fe_zd_corr
      ## note in the header. fe_d_strength is recycled to length 2:
      ## element 1 drives FE1, element 2 drives FE2 (so a single dimension
      ## can be made adversarial while the other stays benign). All-zero =>
      ## `b` below is exactly the previous rep(0, n).
      fed <- rep_len(fe_d_strength, 2L)

      shift_from <- function(zshift, nlev, strength) {
        if (strength == 0) return(numeric(nlev))
        zstd <- zshift / fe_z_strength
        strength * (fe_zd_corr * zstd +
          sqrt(max(0, 1 - fe_zd_corr^2)) * stats::rnorm(nlev))
      }

      fe1_d_shift <- shift_from(fe1_z_shift, nFE1, fed[1])
      fe2_d_shift <- if (has_FE2) shift_from(fe2_z_shift, nFE2, fed[2]) else NULL

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

    } else {
      Z <- sample(0:K, n, replace = TRUE)
    }

    b <- if (has_FE1) {
      bb <- fe1_d_shift[FE1]
      if (has_FE2) bb <- bb + fe2_d_shift[FE2]
      bb
    } else {
      rep(0, n)
    }

  } else if (setup == "B") {
    if (!is.null(z_score_expr) && K == 1L) {
      ## Direct probability-space specification (K == 1, binary Z only):
      ## z_score_expr is evaluated and used *as* P(Z=1|X) itself, the same
      ## convention alpha_good already uses for D -- not as a latent index
      ## to be perturbed with noise and quantile-cut. That makes the true
      ## P(Z=1|X) exactly whatever function of X is written in
      ## z_score_expr, with no threshold-crossing/probit link between the
      ## two, so a linear-probability-model fit (parametric = TRUE) using
      ## the matching functional form can recover it exactly -- letting a
      ## "correct vs misspecified functional form" comparison isolate
      ## regressor completeness alone, uncontaminated by a link-function
      ## mismatch the estimator has no way to get right regardless.
      p_z <- eval_inside(z_score_expr, "z_score_expr")
      p_z <- pmin(pmax(p_z, eps), 1 - eps)

      Z <- stats::rbinom(n, size = 1L, prob = p_z)

    } else {
      z_score <- if (!is.null(z_score_expr)) {
        eval_inside(z_score_expr, "z_score_expr")
      } else {
        0.2 * rowSums(X)
      }
      z_score <- z_score + rnorm(n)

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
    }

    b <- 0.5 * log1p(exp(rowSums(X)))

  } else {
    stop("invalid choice of setup")
  }

  base <- b + errors[, 1]

  alpha_eff_arr <- make_alpha_array(alpha_good, "alpha_good")

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

      alpha_eff_arr[viol, j, k] <- alpha_bad_eval[viol, j, k]
    }
  }

  choose_p0_feasible <- function(
    alpha_arr, J, K,
    eps = 1e-6,
    min_cell_prob = 0.02
  ) {
    n <- dim(alpha_arr)[1]

    delta <- array(0, dim = c(n, J, K + 1L))

    for (k in seq_len(K)) {
      delta[, , k + 1L] <- delta[, , k] + alpha_arr[, , k]
    }

    lower <- matrix(NA_real_, nrow = n, ncol = J)
    upper <- matrix(NA_real_, nrow = n, ncol = J)

    for (j in seq_len(J)) {
      dmin <- apply(delta[, j, , drop = FALSE], 1, min)
      dmax <- apply(delta[, j, , drop = FALSE], 1, max)

      lower[, j] <- eps - dmin
      upper[, j] <- 1 - eps - dmax
    }

    if (any(lower >= upper)) {
      bad <- which(lower >= upper, arr.ind = TRUE)[1, ]

      stop(
        "Alpha values are infeasible for treatment margin j = ",
        bad[2], ". No baseline p0 can keep probabilities inside (0, 1)."
      )
    }

    ## Required baseline gaps between adjacent treatment margins.
    ## Need:
    ##   p0_j + delta_jz >= p0_{j+1} + delta_{j+1,z} + min_cell_prob
    ## for every z.
    ##
    ## Equivalently:
    ##   p0_j - p0_{j+1} >=
    ##   min_cell_prob + max_z(delta_{j+1,z} - delta_{jz})
    req_gap <- NULL

    if (J >= 2L) {
      req_gap <- matrix(NA_real_, nrow = n, ncol = J - 1L)

      for (j in seq_len(J - 1L)) {
        diff_delta <- delta[, j + 1L, , drop = FALSE] -
          delta[, j, , drop = FALSE]

        req_gap[, j] <- min_cell_prob + apply(diff_delta, 1, max)
      }
    }

    p0 <- matrix(NA_real_, nrow = n, ncol = J)

    if (J == 1L) {
      p0[, 1L] <- 0.5 * (lower[, 1L] + upper[, 1L])

    } else {
      ## Deterministic constructive solution.
      ## Start from the lowest feasible final margin, then move upward.
      p0[, J] <- lower[, J]

      for (j in (J - 1L):1L) {
        p0[, j] <- pmax(
          lower[, j],
          p0[, j + 1L] + req_gap[, j]
        )
      }

      ## Check upper bounds.
      too_high <- p0 > upper

      if (any(too_high)) {
        bad <- which(too_high, arr.ind = TRUE)[1, ]

        stop(
          "Could not choose feasible baseline probabilities. ",
          "Requested alpha values require too much probability mass between ",
          "treatment margins. First failure at unit ", bad[1],
          ", treatment margin j = ", bad[2], "."
        )
      }

      ## Optional: add common upward slack where possible to avoid sitting exactly
      ## at lower bounds. This preserves all required gaps.
      slack_up <- apply(upper - p0, 1, min)
      p0 <- p0 + 0.5 * slack_up
    }

    pz <- array(NA_real_, dim = c(n, J, K + 1L))

    for (j in seq_len(J)) {
      for (z in 0:K) {
        pz[, j, z + 1L] <- p0[, j] + delta[, j, z + 1L]
      }
    }

    if (any(pz <= eps | pz >= 1 - eps)) {
      bad_range <- range(pz, finite = TRUE)

      stop(
        "Final probabilities are outside (eps, 1 - eps). ",
        "Range was [", signif(bad_range[1], 4), ", ",
        signif(bad_range[2], 4), "]."
      )
    }

    if (J >= 2L) {
      for (j in seq_len(J - 1L)) {
        gap_j <- pz[, j, ] - pz[, j + 1L, ]

        if (any(gap_j < min_cell_prob - 1e-10)) {
          stop(
            "Final probabilities violate nesting for margins j = ",
            j, " and j + 1 = ", j + 1, "."
          )
        }
      }
    }

    list(
      p0 = p0,
      pz = pz,
      delta = delta,
      lower = lower,
      upper = upper,
      req_gap = req_gap
    )
  }

  make_shift_from_alpha_targets <- function(
    alpha_arr, base, J, K,
    eps = 1e-6,
    min_cell_prob = 0.02
  ) {
    n <- length(base)

    design <- choose_p0_feasible(
      alpha_arr = alpha_arr,
      J = J,
      K = K,
      eps = eps,
      min_cell_prob = min_cell_prob
    )

    p0 <- design$p0
    pz <- design$pz

    cutoffs <- matrix(NA_real_, nrow = n, ncol = J)
    shift <- array(0, dim = c(n, J, K + 1L))

    for (j in seq_len(J)) {
      cutoffs[, j] <- stats::quantile(
        base,
        probs = 1 - p0[, j],
        type = 8,
        names = FALSE
      )

      q <- stats::quantile(
        base,
        probs = as.vector(1 - pz[, j, ]),
        type = 8,
        names = FALSE
      )

      q <- matrix(q, nrow = n, ncol = K + 1L)

      shift[, j, ] <- cutoffs[, j] - q
    }

    list(
      shift = shift,
      cutoffs = cutoffs,
      p0 = p0,
      pz = pz,
      delta = design$delta,
      lower = design$lower,
      upper = design$upper
    )
  }

  shift_obj <- make_shift_from_alpha_targets(
    alpha_arr = alpha_eff_arr,
    base = base,
    J = J,
    K = K,
    eps = eps,
    min_cell_prob = min_cell_prob
  )

  shift <- shift_obj$shift
  cutoffs <- shift_obj$cutoffs

  L <- matrix(NA_real_, nrow = n, ncol = J)

  for (j in seq_len(J)) {
    sh <- shift[cbind(seq_len(n), j, Z + 1L)]

    raw_j <- base + sh - cutoffs[, j]

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

  row_cumsum_with_zero <- function(M) {
    M <- as.matrix(M)

    if (ncol(M) == 0L) {
      return(matrix(0, nrow = nrow(M), ncol = 1L))
    }

    out <- M

    if (ncol(M) >= 2L) {
      for (kk in 2:ncol(M)) {
        out[, kk] <- out[, kk - 1L] + M[, kk]
      }
    }

    cbind(0, out)
  }

  gamma_good_z <- row_cumsum_with_zero(eta_good)
  gamma_bad_z  <- row_cumsum_with_zero(eta_bad)

  gamma_z <- gamma_good_z[cbind(seq_len(n), Z + 1L)]

  if (!is.null(excl_bad) && any(viol)) {
    gamma_z[viol] <- gamma_bad_z[cbind(which(viol), Z[viol] + 1L)]
  }

  tau_D <- rep(0, n)

  for (j in seq_len(J)) {
    tau_D <- tau_D + tau[j] * as.numeric(D >= j)
  }

  ## FE main effects on Y: a generic outcome confounder correlated with FE
  ## and, via Z's own FE-dependence above, unconditionally correlated with
  ## the instrument -- the classic omitted-fixed-effect confound that
  ## controlling for FE in `fml` is meant to absorb.
  fe1_y <- if (has_FE1) stats::rnorm(nFE1, 0, fe_y_sd)[FE1] else 0
  fe2_y <- if (has_FE2) stats::rnorm(nFE2, 0, fe_y_sd)[FE2] else 0

  Y <- as.vector(tau_D + gamma_z + X %*% betaXY + fe1_y + fe2_y + errors[, 2])

  out <- data.frame(Y, D, Z, X)

  if (has_FE1) out$FE1 <- FE1
  if (has_FE2) out$FE2 <- FE2

  if (return_design) {
    expected_first_stage <- matrix(NA_real_, nrow = J, ncol = K)
    rownames(expected_first_stage) <- paste0("D_ge_", seq_len(J))
    colnames(expected_first_stage) <- paste0("Z_margin_", seq_len(K))

    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        expected_first_stage[j, k] <- mean(alpha_eff_arr[, j, k])
      }
    }

    expected_first_stage_by_condition <- NULL

    if (!is.null(condition)) {
      expected_first_stage_by_condition <- list(
        good = matrix(NA_real_, nrow = J, ncol = K),
        bad  = matrix(NA_real_, nrow = J, ncol = K)
      )

      rownames(expected_first_stage_by_condition$good) <- paste0("D_ge_", seq_len(J))
      colnames(expected_first_stage_by_condition$good) <- paste0("Z_margin_", seq_len(K))

      rownames(expected_first_stage_by_condition$bad) <- paste0("D_ge_", seq_len(J))
      colnames(expected_first_stage_by_condition$bad) <- paste0("Z_margin_", seq_len(K))

      for (j in seq_len(J)) {
        for (k in seq_len(K)) {
          expected_first_stage_by_condition$good[j, k] <-
            if (any(!viol)) mean(alpha_eff_arr[!viol, j, k]) else NA_real_

          expected_first_stage_by_condition$bad[j, k] <-
            if (any(viol)) mean(alpha_eff_arr[viol, j, k]) else NA_real_
        }
      }
    }

    realized_first_stage <- matrix(NA_real_, nrow = J, ncol = K)
    rownames(realized_first_stage) <- paste0("D_ge_", seq_len(J))
    colnames(realized_first_stage) <- paste0("Z_margin_", seq_len(K))

    for (j in seq_len(J)) {
      Dj <- as.numeric(D >= j)

      for (k in seq_len(K)) {
        left  <- Z == k - 1L
        right <- Z == k

        realized_first_stage[j, k] <-
          if (any(left) && any(right)) {
            mean(Dj[right]) - mean(Dj[left])
          } else {
            NA_real_
          }
      }
    }

    realized_first_stage_by_condition <- NULL

    if (!is.null(condition)) {
      realized_first_stage_by_condition <- list(
        good = matrix(NA_real_, nrow = J, ncol = K),
        bad  = matrix(NA_real_, nrow = J, ncol = K)
      )

      rownames(realized_first_stage_by_condition$good) <- paste0("D_ge_", seq_len(J))
      colnames(realized_first_stage_by_condition$good) <- paste0("Z_margin_", seq_len(K))

      rownames(realized_first_stage_by_condition$bad) <- paste0("D_ge_", seq_len(J))
      colnames(realized_first_stage_by_condition$bad) <- paste0("Z_margin_", seq_len(K))

      for (j in seq_len(J)) {
        Dj <- as.numeric(D >= j)

        for (k in seq_len(K)) {
          left_good  <- !viol & Z == k - 1L
          right_good <- !viol & Z == k

          left_bad  <- viol & Z == k - 1L
          right_bad <- viol & Z == k

          realized_first_stage_by_condition$good[j, k] <-
            if (any(left_good) && any(right_good)) {
              mean(Dj[right_good]) - mean(Dj[left_good])
            } else {
              NA_real_
            }

          realized_first_stage_by_condition$bad[j, k] <-
            if (any(left_bad) && any(right_bad)) {
              mean(Dj[right_bad]) - mean(Dj[left_bad])
            } else {
              NA_real_
            }
        }
      }
    }

    attr(out, "design") <- list(
      expected_first_stage = expected_first_stage,
      expected_first_stage_by_condition = expected_first_stage_by_condition,
      realized_first_stage = realized_first_stage,
      realized_first_stage_by_condition = realized_first_stage_by_condition,
      p0 = shift_obj$p0,
      pz = shift_obj$pz,
      delta = shift_obj$delta,
      lower = shift_obj$lower,
      upper = shift_obj$upper,
      alpha_eff_arr = alpha_eff_arr,
      cutoffs = shift_obj$cutoffs,
      realized_D_dist = prop.table(table(factor(D, levels = 0:J))),
      realized_Z_dist = prop.table(table(factor(Z, levels = 0:K))),
      condition_share = mean(viol),
      min_cell_prob = min_cell_prob
    )
  }

  out
}
