## ARCHIVED prototype (from ../montest_wts, never in the package): train-side power-optimal weighted
## Bonferroni over the candidate cells of all pool/select designs, per sample direction (sample tested
## separately, budget alpha/2 per direction), p.wbonf = min(1, p.raw / w). See NOTE.md.
##
## Needs: forest_test_core()/forest_test() run with defer_selection = TRUE and return_u = TRUE (the
## influence vectors `train_u` now live in R/functions.R together with cell_corr_from_u() and
## power_union(), which screen = "power" uses), and the hook below inserted into the legacy wrapper
## (dev/adaptive/forest_test_adaptive.R) after the design grid is built. montest() step 7 then needs a
## branch that sets minp = c(p.raw, p.wbonf, p.wbonf_hi, p.equal) from the p.wbonf columns.

## ---- hook inside the legacy wrapper ----
  if (identical(getOption("montest.adaptive_method", "legacy"), "weights")) {
    ## Train-side power-optimal weighted Bonferroni over all candidate cells of all designs,
    ## per sample direction; sample is tested separately (no cross-half pooling).
    names_w <- character(nrow(design_grid))
    fits_w <- vector("list", nrow(design_grid))
    for (ii in seq_len(nrow(design_grid))) {
      eff <- make_eff(design_grid[ii])
      names_w[ii] <- design_name(design_grid[ii])
      fits_w[[ii]] <- run_core(setdiff(eff$pool, sample_col), eff$select, defer = TRUE, ret_u = TRUE)
    }
    out <- adaptive_weights_fit(fits_w, names_w, alpha = alpha, sample_col = sample_col)
    out$adaptive_dims <- adaptive_dims
    return(out)
  }

## ---- the fit ----
## ---------------------------------------------------------------------------------------
## Train-side power-optimal weighted Bonferroni over candidate cells (all designs pooled).
## `fits`: forest_test_core() outputs run with defer_selection = TRUE and return_u = TRUE,
## sample tested separately. Per sample direction s (train on half s, test on the other):
## candidates = all train cells of all designs; predictive held-out statistics
## U = -T ~ N(a, (1 + phi) R), a = -t_train, R from cluster-level influence vectors;
## weights w' (sum 1) maximize P(some U_k >= qnorm(1 - (alpha/2) w'_k)); final weights are
## w = w'/2 (the budget is split equally across directions, so no cross-half leakage).
adaptive_weights_fit <- function(fits, design_names, alpha = 0.05, sample_col = "sample",
                                 phi = getOption("montest.phi", 1), phi_hi = getOption("montest.phi_hi", 4), maxcand = 8L) {
  stat_cols <- c("cell_id", "train_row_id", "train", "relevant", "G", "N", "dof_rank", "df",
                 "coef", "stderr", "t", "tau_cutoff", "p.raw", ".uidx", "design_id", sample_col)
  cand <- data.table::rbindlist(lapply(seq_along(fits), function(d) {
    tc <- data.table::copy(fits[[d]]$train_cells)
    tc[, `:=`(design_id = d, .uidx = seq_len(.N))]
    tc
  }), use.names = TRUE, fill = TRUE)
  key_cols <- setdiff(names(cand), stat_cols)

  power_fun <- function(w, a, Sig, alpha_dir) {
    on <- w > 1e-9
    cc <- stats::qnorm(1 - alpha_dir * pmin(w[on], 1))
    if (sum(on) == 1L) return(1 - stats::pnorm(cc, a[on], sqrt(Sig[on, on])))
    ## Genz-Bretz is randomized: fix the seed per call (common random numbers) so the
    ## objective is smooth for the optimizer; the caller's RNG state is restored below.
    set.seed(20240101L)
    1 - mvtnorm::pmvnorm(upper = cc, mean = a[on], sigma = Sig[on, on, drop = FALSE],
                         algorithm = mvtnorm::GenzBretz(abseps = 1e-4, maxpts = 20000L))[1L]
  }
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  }

  sel_list <- list(); R_list <- list()
  for (s_dir in c(1L, 2L)) {
    cs <- cand[get(sample_col) == s_dir & is.finite(t)]
    if (nrow(cs) == 0L) next
    cs[, a := -t]
    data.table::setorder(cs, -a)
    cs <- cs[seq_len(min(.N, 3L * maxcand))]
    ## influence matrix over clusters
    ul <- lapply(seq_len(nrow(cs)), function(i) fits[[cs$design_id[i]]]$train_u[[cs$.uidx[i]]])
    all_cl <- sort(unique(unlist(lapply(ul, function(z) names(z$u)))))
    U <- matrix(0, length(all_cl), length(ul))
    for (i in seq_along(ul)) U[match(names(ul[[i]]$u), all_cl), i] <- ul[[i]]$u
    sc <- vapply(ul, function(z) sqrt(z$G / z$df), numeric(1))
    V <- (crossprod(U) * outer(sc, sc))
    se_chk <- vapply(ul, function(z) z$se, numeric(1))
    bad <- abs(sqrt(diag(V)) / se_chk - 1) > 1e-4
    if (any(bad, na.rm = TRUE)) {
      warning("adaptive weights: influence-based SE differs from reported SE for ",
              sum(bad, na.rm = TRUE), " candidate(s); using reported SE scale.", call. = FALSE)
    }
    R <- stats::cov2cor(V + diag(1e-8, nrow(V)))
    ## drop near-duplicate candidates (identical cells), keep the stronger one
    keep <- rep(TRUE, nrow(cs))
    for (i in seq_len(nrow(cs))) if (keep[i]) {
      dup <- which(seq_len(nrow(cs)) > i & abs(R[i, ]) > 0.9999)
      keep[dup] <- FALSE
    }
    keep[cumsum(keep) > maxcand] <- FALSE
    cs <- cs[keep]; R <- R[keep, keep, drop = FALSE]
    K <- nrow(cs)
    a <- cs$a
    alpha_dir <- alpha / 2
    solve_w <- function(Sig) {
      if (K == 1L) return(1)
      starts <- c(list(rep(0, K)), lapply(seq_len(min(K, 2L)), function(k) replace(rep(-2, K), k, 2)))
      best <- NULL
      for (st in starts) {
        o <- stats::optim(st, function(par) -power_fun(exp(par) / sum(exp(par)), a, Sig, alpha_dir),
                          method = "Nelder-Mead", control = list(maxit = 200))
        if (is.null(best) || o$value < best$value) best <- o
      }
      w <- exp(best$par) / sum(exp(best$par))
      w[w < 0.005] <- 0
      w / sum(w)
    }
    Sig <- (1 + phi) * R
    w <- solve_w(Sig)
    w_hi <- solve_w((1 + phi_hi) * R)
    w_eq <- rep(1 / K, K)
    cs[, `:=`(w = w / 2, w_hi = w_hi / 2, w_eq = w_eq / 2,
              pred_power = power_fun(w, a, Sig, alpha_dir), direction = s_dir)]
    R_list[[s_dir]] <- list(cells = cs[, c(key_cols, "design_id", "t"), with = FALSE], R = R)
    sel_list[[length(sel_list) + 1L]] <- cs
  }
  sel <- data.table::rbindlist(sel_list, use.names = TRUE, fill = TRUE)
  if (nrow(sel) == 0L) stop("adaptive weights: no candidate cells.")
  chosen <- sel[w > 0 | w_hi > 0 | w_eq > 0]

  ## held-out rows of the chosen candidates: same design, same cell key, opposite sample
  tests <- data.table::rbindlist(lapply(seq_len(nrow(chosen)), function(i) {
    r <- chosen[i]
    tr <- data.table::as.data.table(fits[[r$design_id]]$results)[train == FALSE]
    ok <- tr[[sample_col]] == (3L - r$direction)
    for (cc in intersect(key_cols, names(tr))) {
      v <- r[[cc]]
      if (!is.na(v)) ok <- ok & (tr[[cc]] == v)
    }
    hit <- tr[ok]
    if (nrow(hit) != 1L) {
      stop("adaptive weights: expected one held-out row per chosen cell, found ", nrow(hit), ".")
    }
    hit[, `:=`(design = design_names[r$design_id], w = r$w, w_hi = r$w_hi, w_eq = r$w_eq, direction = r$direction)]
    hit
  }), use.names = TRUE, fill = TRUE)
  tests[, p.wbonf := pmin(1, p.raw / w)]
  tests[, p.wbonf_hi := pmin(1, p.raw / w_hi)]
  tests[, p.equal := pmin(1, p.raw / w_eq)]
  tests[, relevant := as.integer(w > 0)]

  trains <- sel[, c(setdiff(names(sel), c(".uidx", "a", "pred_power", "direction", "w", "w_hi", "w_eq", "design_id")),
                    "w", "design_id"), with = FALSE]
  trains[, `:=`(train = TRUE, relevant = as.integer(w > 0), design = design_names[design_id], design_id = NULL)]
  results <- data.table::rbindlist(list(trains, tests), use.names = TRUE, fill = TRUE)

  list(
    results = results, grid = NULL, global = NULL, Xmeans = NULL, Xmeans_all = NULL, XSD = NULL,
    sample_design = "separate_weighted", sample_pool_after_margin_select = FALSE,
    adaptive_design = "weights", adaptive_method = "weights", adaptive_R = R_list,
    adaptive_weights = sel[, c(key_cols, "design_id", "direction", "t", "w", "pred_power"), with = FALSE]
  )
}

