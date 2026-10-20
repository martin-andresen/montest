## ARCHIVED (removed from R/functions.R). Adaptive choice between pooling and selecting over a margin:
## the legacy forest_test() wrapper that fitted all 2^k pool/select designs with forest_test_core()
## and kept the winner by a train-side score. See NOTE.md. Last commit containing it in the package:
## git tag archive/adaptive-pool-select. It calls forest_test_core() (renamed forest_test() after the removal).

####### FIND OPTIMAL CUTOFF AND TEST #############
forest_test <- function(
    data,
    cluster = NULL,
    weight  = NULL,
    sample  = "sample",
    pred    = "pred",
    pred_o  = "pred_o",
    sort_pred = NULL,
    sort_pred_o = NULL,
    scores  = "scores",
    scores_test = NULL,
    x_names = NULL,
    minsize = 50L,
    arm = NULL,
    minarm = 10L,
    minarm_clusters = 5L,
    margins = NULL,
    pool = NULL,
    select = NULL,
    gridpoints = NULL,
    store_grid = TRUE,
    verbose = FALSE,
    screen = c("stepdown", "none", "minimum", "negative", "nonpositive", "fgk_relevant", "power"),
    alpha = 0.05,
    fe_expr = NULL,
    fe_rank_adj = !is.null(fe_expr),
    fe_rank_conservative = TRUE,
    x_rank_vars = character(0),
    sandwich = NULL,
    center = NULL,
    resid_treat = NULL,
    resid_outcome = NULL,
    sample_weight = NULL,
    recenter_propensity = FALSE,
    recenter_binary = FALSE,
    v = NULL,
    center_inv_v = FALSE,
    defer_selection = FALSE
) {
  screen <- match.arg(screen)

  if (is.null(pool)) pool <- character()
  if (is.null(select)) select <- character()
  if (is.null(margins)) margins <- character()

  pool <- unique(as.character(pool))
  select <- unique(as.character(select))
  margins <- as.character(margins)

  sample_col <- as.character(sample)
  stopifnot(length(sample_col) == 1L)

  adaptive_dims <- intersect(pool, select)

  if (isTRUE(defer_selection) && length(adaptive_dims) > 0L) {
    stop("Internal error: adaptive pool/select designs are incompatible with deferred selection.")
  }

  if (sample_col %in% adaptive_dims) {
    stop(
      "`sample` cannot appear in both `pool` and `select`.\n",
      "This would adaptively choose between sample-pooled and sample-separated ",
      "testing using both sample halves. Put `sample` in either `pool` or `select`, ",
      "but not both."
    )
  }

  design_score <- function(fit, pool_eff, select_eff, design_row) {
    ## All candidate train cells (not just screen survivors): K = number of cells the
    ## design searched over, so a design that selects over more margins pays for it, and an
    ## empty screen is not scored as an uncorrected p. Cutoff-grid search is common to all
    ## designs and cancels in the ranking.
    tr <- data.table::as.data.table(fit$results)[train == TRUE]

    if (nrow(tr) == 0L) return(Inf)
    if (!("t" %in% names(tr))) return(Inf)

    ## Log scale: pnorm(t) underflows to 0 for t < ~-38, which would tie all designs.
    ## The smallest Holm-adjusted p is K * min(p), i.e. log(K) + min(log p), capped at 1.
    tt <- tr$t
    tt <- tt[is.finite(tt)]

    if (!length(tt)) return(Inf)

    min(0, log(length(tt)) + min(stats::pnorm(tt, log.p = TRUE)))
  }

  run_core <- function(pool_eff, select_eff) {
    forest_test_core(
      data = data,
      cluster = cluster,
      weight = weight,
      sample = sample,
      pred = pred,
      pred_o = pred_o,
      sort_pred = sort_pred,
      sort_pred_o = sort_pred_o,
      scores = scores,
      scores_test = scores_test,
      x_names = x_names,
      minsize = minsize,
      arm = arm,
      minarm = minarm,
      minarm_clusters = minarm_clusters,
      margins = margins,
      pool = pool_eff,
      select = select_eff,
      gridpoints = gridpoints,
      store_grid = store_grid,
      verbose = verbose,
      screen = screen,
      alpha = alpha,
      fe_expr = fe_expr,
      fe_rank_adj = fe_rank_adj,
      fe_rank_conservative = fe_rank_conservative,
      x_rank_vars = x_rank_vars,
      sandwich = sandwich,
      center = center,
      resid_treat = resid_treat,
      resid_outcome = resid_outcome,
      sample_weight = sample_weight,
      recenter_propensity = recenter_propensity,
      recenter_binary = recenter_binary,
      v = v,
      center_inv_v = center_inv_v,
      defer_selection = defer_selection
    )
  }

  if (length(adaptive_dims) == 0L) {
    out <- run_core(pool, select)
    out$adaptive_dims <- character()
    out$adaptive_design <- NA_character_
    out$adaptive_scores <- NULL
    return(out)
  }

  pool_base <- setdiff(pool, adaptive_dims)
  select_base <- setdiff(select, adaptive_dims)

  design_opts_by_dim <- lapply(adaptive_dims, function(dd) {
      c(
        "pooled_before_search",
        "selected_then_separate"
      )
  })

  names(design_opts_by_dim) <- adaptive_dims

  design_grid <- data.table::as.data.table(
    expand.grid(design_opts_by_dim, stringsAsFactors = FALSE)
  )

  make_eff <- function(design_row) {
    pool_eff <- pool_base
    select_eff <- select_base

    for (dd in adaptive_dims) {
      choice <- as.character(design_row[[dd]])



        if (choice == "pooled_before_search") {
          ## Pool this non-sample margin before the forest/grid search.
          pool_eff <- c(pool_eff, dd)

        } else if (choice == "selected_then_separate") {
          ## Select over this non-sample margin and keep final tests separate.
          select_eff <- c(select_eff, dd)

        } else if (choice == "selected_then_pooled") {
          stop(
            "Internal error: selected_then_pooled is no longer allowed for non-sample margins."
          )

        } else {
          stop("Unknown adaptive design choice: ", choice)
        }
    }

    pool_eff <- unique(pool_eff)
    select_eff <- unique(select_eff)

    ## Non-sample overlap should never be passed to the core.
    bad_overlap <- setdiff(intersect(pool_eff, select_eff), sample_col)
    if (length(bad_overlap) > 0L) {
      stop(
        "Internal error: adaptive design produced non-sample overlap: ",
        paste(bad_overlap, collapse = ", ")
      )
    }

    ## Sample overlap should also never be passed to the core.
    if (sample_col %in% pool_eff && sample_col %in% select_eff) {
      stop("Internal error: adaptive design produced sample in both pool and select.")
    }

    list(
      pool = pool_eff,
      select = select_eff
    )
  }

  design_name <- function(design_row) {
    paste(
      paste(names(design_row), as.character(unlist(design_row)), sep = "="),
      collapse = "; "
    )
  }

  fits <- vector("list", nrow(design_grid))
  scores_design <- rep(Inf, nrow(design_grid))
  names_design <- character(nrow(design_grid))
  pool_designs <- vector("list", nrow(design_grid))
  select_designs <- vector("list", nrow(design_grid))

  for (ii in seq_len(nrow(design_grid))) {
    eff <- make_eff(design_grid[ii])
    names_design[ii] <- design_name(design_grid[ii])
    pool_designs[[ii]] <- eff$pool
    select_designs[[ii]] <- eff$select

    fits[[ii]] <- run_core(eff$pool, eff$select)

    scores_design[ii] <- design_score(
      fit = fits[[ii]],
      pool_eff = eff$pool,
      select_eff = eff$select,
      design_row = design_grid[ii]
    )
  }

  flex_score <- vapply(names_design, function(x) {
    parts <- strsplit(x, "; ", fixed = TRUE)[[1]]

    sum(vapply(parts, function(p) {
      dim <- sub("=.*$", "", p)
      val <- sub("^[^=]+=", "", p)

      if (dim == sample_col && val == "selected_then_pooled") return(0L)
      if (dim == sample_col && val == "selected_then_separate") return(1L)

      if (val == "pooled_before_search") return(0L)
      if (val == "selected_then_separate") return(1L)

      99L
    }, integer(1)))
  }, integer(1))

  winner_i <- which.min(scores_design + 1e-12 * flex_score)

  out <- fits[[winner_i]]

  out$adaptive_dims <- adaptive_dims
  out$adaptive_design <- names_design[winner_i]
  out$adaptive_scores <- data.table::data.table(
    design_id = seq_along(names_design),
    design = names_design,
    score = scores_design,
    winner = seq_along(names_design) == winner_i,
    pool = vapply(pool_designs, function(z) paste(z, collapse = ","), character(1L)),
    select = vapply(select_designs, function(z) paste(z, collapse = ","), character(1L))
  )[order(score)]

  out$adaptive_pool <- pool_designs[[winner_i]]
  out$adaptive_select <- select_designs[[winner_i]]

  out
}
