

##################################################### HELPER FUNCTIONS ######################################
shrink_te_crossfit <- function(data,
                               pred,
                               pred_var,
                               pred_out,
                               pred_out_var,
                               margins = NULL,
                               sample = "sample",
                               weight = NULL,
                               gamma = 1,
                               out_pred = pred,
                               out_pred_out = pred_out) {

  stopifnot(data.table::is.data.table(data))

  pred_col         <- as.character(pred)
  pred_var_col     <- as.character(pred_var)
  pred_out_col     <- as.character(pred_out)
  pred_out_var_col <- as.character(pred_out_var)
  sample_col       <- as.character(sample)
  margins          <- as.character(margins)
  out_pred_col     <- as.character(out_pred)
  out_pred_out_col <- as.character(out_pred_out)
  if (!is.null(weight)) weight_col <- as.character(weight) else weight_col <- NULL

  stopifnot(length(pred_col) == 1L, pred_col %in% names(data))
  stopifnot(length(pred_var_col) == 1L, pred_var_col %in% names(data))
  stopifnot(length(pred_out_col) == 1L, pred_out_col %in% names(data))
  stopifnot(length(pred_out_var_col) == 1L, pred_out_var_col %in% names(data))
  stopifnot(length(sample_col) == 1L, sample_col %in% names(data))
  if (length(margins) > 0L) stopifnot(all(margins %in% names(data)))
  if (!is.null(weight_col)) stopifnot(length(weight_col) == 1L, weight_col %in% names(data))

  stopifnot(is.numeric(gamma), length(gamma) == 1L, is.finite(gamma), gamma >= 0, gamma <= 1)

  if (!all(data[[sample_col]] %in% c(1L, 2L, 1, 2))) {
    stop("sample indicator must take values 1 and 2.")
  }

  DT <- data
  byvars <- c(margins, sample_col)

  wmean <- function(x, w) sum(w * x) / sum(w)

  wvar_pop <- function(x, w) {
    mu <- wmean(x, w)
    sum(w * (x - mu)^2) / sum(w)
  }

  eb_pars <- DT[, {
    y <- get(pred_col)
    v <- get(pred_var_col)
    ww <- if (is.null(weight_col)) rep(1, .N) else get(weight_col)

    ok <- is.finite(y) & is.finite(v) & v >= 0 & is.finite(ww) & ww > 0
    y  <- y[ok]
    v  <- v[ok]
    ww <- ww[ok]

    if (length(y) == 0L) {
      .(eb_mu = NA_real_, eb_tau2 = NA_real_)
    } else if (length(y) == 1L) {
      .(eb_mu = y[1L], eb_tau2 = 0)
    } else {
      mu_hat <- wmean(y, ww)
      v_y    <- wvar_pop(y, ww)
      m_v    <- wmean(v, ww)
      tau2   <- max(0, v_y - m_v)
      .(eb_mu = mu_hat, eb_tau2 = tau2)
    }
  }, by = byvars]

  DT[eb_pars, c("..eb_mu", "..eb_tau2") := .(i.eb_mu, i.eb_tau2), on = byvars]

  DT[, (out_pred_col) := {
    y  <- get(pred_col)
    v  <- get(pred_var_col)
    mu <- ..eb_mu
    t2 <- ..eb_tau2

    out <- y
    ok <- is.finite(y) & is.finite(v) & v >= 0 & is.finite(mu) & is.finite(t2)

    if (any(ok)) {
      B  <- t2[ok] / (t2[ok] + v[ok])
      eb <- mu[ok] + B * (y[ok] - mu[ok])
      out[ok] <- (1 - gamma) * y[ok] + gamma * eb
    }
    out
  }]

  eb_pars_out <- data.table::copy(eb_pars)
  eb_pars_out[, (sample_col) := ifelse(get(sample_col) == 1L, 2L, 1L)]

  DT[eb_pars_out, c("..eb_mu_out", "..eb_tau2_out") := .(i.eb_mu, i.eb_tau2), on = byvars]

  DT[, (out_pred_out_col) := {
    y  <- get(pred_out_col)
    v  <- get(pred_out_var_col)
    mu <- ..eb_mu_out
    t2 <- ..eb_tau2_out

    out <- y
    ok <- is.finite(y) & is.finite(v) & v >= 0 & is.finite(mu) & is.finite(t2)

    if (any(ok)) {
      B  <- t2[ok] / (t2[ok] + v[ok])
      eb <- mu[ok] + B * (y[ok] - mu[ok])
      out[ok] <- (1 - gamma) * y[ok] + gamma * eb
    }
    out
  }]

  DT[, c("..eb_mu", "..eb_tau2", "..eb_mu_out", "..eb_tau2_out") := NULL]

  invisible(DT)
}

CART_test <- function(
    data,
    sample  = "sample",
    scores  = "scores",
    x_names,
    margins = NULL,
    weight  = NULL,
    cluster = NULL,
    select = NULL,            # subset of margins; may also include "sample"
    cp = 0.0,
    maxrankcp = 10L,
    alpha = 0.05,
    prune = TRUE,
    minsize = 50L,
    preselect = c("none", "minimum", "negative", "nonpositive", "fgk_relevant"),
    store_trees = FALSE,
    verbose = FALSE,
    xval=10,
    rpart_options=NULL
) {
  stopifnot(data.table::is.data.table(data))
  preselect <- match.arg(preselect)

  sample_col <- as.character(sample)
  scores_col <- as.character(scores)
  stopifnot(length(sample_col) == 1L, sample_col %chin% names(data))
  stopifnot(length(scores_col) == 1L, scores_col %chin% names(data))

  x_names <- as.character(x_names)
  stopifnot(length(x_names) >= 1L, all(x_names %chin% names(data)))

  if (is.null(margins)) margins <- character()
  margins <- unique(as.character(margins))
  if (length(margins) > 0L) stopifnot(all(margins %chin% names(data)))

  weight_col <- if (is.null(weight)) NULL else as.character(weight)
  if (!is.null(weight_col)) stopifnot(length(weight_col) == 1L, weight_col %chin% names(data))

  cluster_col <- if (is.null(cluster)) NULL else as.character(cluster)
  if (!is.null(cluster_col)) stopifnot(length(cluster_col) == 1L, cluster_col %chin% names(data))

  stopifnot(is.numeric(minsize), length(minsize) == 1L, is.finite(minsize), minsize >= 1)
  minsize <- as.integer(minsize)

  stopifnot(is.numeric(alpha), length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  stopifnot(is.numeric(cp), length(cp) == 1L, is.finite(cp), cp >= 0)
  stopifnot(is.numeric(maxrankcp), length(maxrankcp) == 1L, is.finite(maxrankcp), maxrankcp >= 1)
  maxrankcp <- as.integer(maxrankcp)

  svec <- as.integer(data[[sample_col]])
  stopifnot(all(svec %in% c(1L, 2L)))

  # ---- selection setup ----
  if (is.null(select)) select <- character()
  select <- unique(as.character(select))
  select[select == "sample"] <- sample_col

  allowed_select <- unique(c(margins, sample_col))
  bad_select <- setdiff(select, allowed_select)
  if (length(bad_select) > 0L) {
    stop("select contains invalid names: ", paste(bad_select, collapse = ", "))
  }

  select_sample  <- sample_col %chin% select
  select_margins <- intersect(select, margins)
  cell_cols      <- setdiff(margins, select_margins)   # margins still kept separate

  # ---------------- helpers ----------------

  mean_test_crv1 <- function(y, w = NULL, cl = NULL) {
    y <- as.numeric(y)
    n <- length(y)
    if (n == 0L) return(list(coef=NA_real_, se=NA_real_, t=NA_real_, N=0L, G=0L))

    if (is.null(w)) w <- rep(1.0, n) else w <- as.numeric(w)
    ok <- is.finite(y) & is.finite(w) & (w >= 0)
    if (!any(ok)) return(list(coef=NA_real_, se=NA_real_, t=NA_real_, N=0L, G=0L))
    y <- y[ok]
    w <- w[ok]

    if (is.null(cl)) {
      cl <- seq_len(length(y))
    } else {
      cl <- cl[ok]
      cl <- as.integer(factor(cl, exclude = NULL))
    }

    dt0 <- data.table::data.table(y = y, w = w, cl = cl)
    gb <- dt0[, .(U = sum(w * y), W = sum(w)), by = cl]
    G <- nrow(gb)
    U <- sum(gb$U)
    W <- sum(gb$W)
    if (!is.finite(W) || W == 0) {
      return(list(coef=NA_real_, se=NA_real_, t=NA_real_, N=length(y), G=G))
    }

    theta <- U / W
    ug <- gb$U - theta * gb$W
    se <- sqrt((G / pmax.int(G - 1L, 1L)) * sum(ug^2)) / abs(W)
    if (G < 2L || !is.finite(se)) se <- NA_real_
    list(coef = theta, se = se, t = theta / se, N = length(y), G = G)
  }

  select_leaves_local <- function(dt_train_leaf) {
    L <- nrow(dt_train_leaf)
    if (L == 0L) return(integer())

    if (preselect == "none") {
      return(as.integer(dt_train_leaf$leaf))
    }

    if (preselect == "minimum") {
      return(as.integer(dt_train_leaf[which.min(t), leaf]))
    }

    if (preselect == "negative") {
      thr <- stats::qnorm(alpha / L)
      keep <- dt_train_leaf[t <= thr, leaf]
      if (length(keep) == 0L) keep <- dt_train_leaf[which.min(t), leaf]
      return(as.integer(keep))
    }

    if (preselect == "nonpositive") {
      thr <- stats::qnorm(1 - alpha / L)
      keep <- dt_train_leaf[t < thr, leaf]
      if (length(keep) == 0L) keep <- dt_train_leaf[which.min(t), leaf]
      return(as.integer(keep))
    }

    if (preselect == "fgk_relevant") {
      # Mirror of FGK's relevance screening for left-tailed tests:
      # keep leaves that are not strongly in the wrong direction.
      thr <- stats::qnorm(1 - alpha / L)
      keep <- dt_train_leaf[t < thr, leaf]
      return(as.integer(keep))
    }

    as.integer(dt_train_leaf$leaf)
  }

  # Choose ONE cell (across selected margins / maybe sample), then choose leaves within it.
  choose_group_selection <- function(dt_grp) {
    dtf <- dt_grp[is.finite(t)]
    if (nrow(dtf) == 0L) {
      return(list(
        best_cell = NA_integer_,
        best_train_s = NA_integer_,
        sel_leaves = integer()
      ))
    }

    if (preselect == "minimum") {
      j <- which.min(dtf$t)
      return(list(
        best_cell = as.integer(dtf$cell_id[j]),
        best_train_s = as.integer(dtf$train_s[j]),
        sel_leaves = as.integer(dtf$leaf[j])
      ))
    }

    if (preselect == "none") {
      j <- which.min(dtf$t)
      bc <- as.integer(dtf$cell_id[j])
      bt <- as.integer(dtf$train_s[j])
      return(list(
        best_cell = bc,
        best_train_s = bt,
        sel_leaves = as.integer(dtf[cell_id == bc & train_s == bt, leaf])
      ))
    }

    if (preselect == "negative") {
      thr <- stats::qnorm(alpha / nrow(dtf))
      keep <- dtf[t <= thr]
      if (nrow(keep) == 0L) {
        j <- which.min(dtf$t)
        return(list(
          best_cell = as.integer(dtf$cell_id[j]),
          best_train_s = as.integer(dtf$train_s[j]),
          sel_leaves = as.integer(dtf$leaf[j])
        ))
      } else {
        j <- which.min(keep$t)
        bc <- as.integer(keep$cell_id[j])
        bt <- as.integer(keep$train_s[j])
        return(list(
          best_cell = bc,
          best_train_s = bt,
          sel_leaves = as.integer(keep[cell_id == bc & train_s == bt, leaf])
        ))
      }
    }

    if (preselect == "nonpositive") {
      thr <- stats::qnorm(1 - alpha / nrow(dtf))
      keep <- dtf[t < thr]
      if (nrow(keep) == 0L) {
        j <- which.min(dtf$t)
        return(list(
          best_cell = as.integer(dtf$cell_id[j]),
          best_train_s = as.integer(dtf$train_s[j]),
          sel_leaves = as.integer(dtf$leaf[j])
        ))
      } else {
        j <- which.min(keep$t)
        bc <- as.integer(keep$cell_id[j])
        bt <- as.integer(keep$train_s[j])
        return(list(
          best_cell = bc,
          best_train_s = bt,
          sel_leaves = as.integer(keep[cell_id == bc & train_s == bt, leaf])
        ))
      }
    }

    if (preselect == "fgk_relevant") {
      # After local FGK screening, choose the most negative remaining cell/sample,
      # and test all locally-screened leaves in that chosen cell/sample.
      j <- which.min(dtf$t)
      bc <- as.integer(dtf$cell_id[j])
      bt <- as.integer(dtf$train_s[j])
      return(list(
        best_cell = bc,
        best_train_s = bt,
        sel_leaves = as.integer(dtf[cell_id == bc & train_s == bt, leaf])
      ))
    }

    j <- which.min(dtf$t)
    list(
      best_cell = as.integer(dtf$cell_id[j]),
      best_train_s = as.integer(dtf$train_s[j]),
      sel_leaves = as.integer(dtf$leaf[j])
    )
  }

  fit_one_tree_cell <- function(df_cell, train_sample) {
    dtr <- df_cell[df_cell[[sample_col]] == train_sample]
    if (nrow(dtr) < 2L) return(NULL)

    cols_tr <- c(scores_col, x_names)
    df_tr <- as.data.frame(dtr[, ..cols_tr])

    w_tr <- if (!is.null(weight_col)) as.numeric(dtr[[weight_col]]) else NULL
    fml <- stats::as.formula(paste0(scores_col, " ~ ", paste(x_names, collapse = " + ")))


    tree <- rpart::rpart(
      formula = fml,
      data = df_tr,
      method = "anova",
      weights = w_tr,
      control = do.call(rpart::rpart.control, modifyList(rpart_options %||% list(),
                                                         list(cp = cp, minbucket = minsize)))
    )

    if (!is.null(tree$cp) && nrow(tree$cp) > 0L) {
      maxrankcp2 <- min(maxrankcp, nrow(tree$cptable))
      maxcp <- tree$cp[maxrankcp2, 1]
      if (isTRUE(prune)) {
        opcpid <- which.min(tree$cptable[, 4])
        opcp <- tree$cptable[opcpid, 1]
        tree <- rpart::prune(tree, cp = opcp)
      } else {
        tree <- rpart::prune(tree, cp = maxcp)
      }
    }

    df_all <- as.data.frame(df_cell[, ..x_names])
    leaf_all <- as.integer(rpart:::pred.rpart(tree, as.matrix(df_all)))
    list(tree = tree, leaf_all = leaf_all)
  }

  global_means_one_cell <- function(df_cell, key_dt) {
    dtg <- data.table::data.table(
      sample = as.integer(df_cell[[sample_col]]),
      score  = as.numeric(df_cell[[scores_col]]),
      w      = if (!is.null(weight_col)) as.numeric(df_cell[[weight_col]]) else rep(1.0, nrow(df_cell))
    )
    dtg$w[!is.finite(dtg$w)] <- 0
    if (!is.null(cluster_col)) {
      dtg[, cl := as.integer(factor(df_cell[[cluster_col]], exclude = NULL))]
    } else {
      dtg[, cl := .I]
    }

    out <- dtg[, {
      o <- mean_test_crv1(score, w, cl)
      data.table::data.table(
        train = FALSE,
        G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
        p.raw = stats::pnorm(o$t)
      )
    }, by = "sample"]

    if (!is.null(key_dt) && ncol(key_dt) > 0L) {
      for (cc in names(key_dt)) out[, (cc) := key_dt[[cc]][1]]
      data.table::setcolorder(out, c(names(key_dt), "train", "sample"))
    }
    out
  }

  # ---------------- build full-margin cells ----------------
  n <- nrow(data)
  gid <- if (length(margins) == 0L) rep(1L, n) else data.table::frankv(data[, ..margins], ties.method = "dense")
  idx_list <- split(seq_len(n), gid)
  n_cells <- length(idx_list)

  cell_objs <- vector("list", n_cells)
  global_list <- vector("list", n_cells)
  trees_list <- if (store_trees) vector("list", n_cells) else NULL

  for (g in seq_len(n_cells)) {
    idx <- idx_list[[g]]
    key_dt <- if (length(margins) > 0L) data[idx[1L], ..margins] else NULL
    df_cell <- data[idx]

    global_list[[g]] <- global_means_one_cell(df_cell, key_dt)
    if (store_trees) trees_list[[g]] <- vector("list", 2L)

    dirs <- vector("list", 2L)

    for (train_s in c(1L, 2L)) {
      est_s <- if (train_s == 1L) 2L else 1L

      fit <- fit_one_tree_cell(df_cell, train_s)
      if (is.null(fit)) {
        dirs[[train_s]] <- NULL
        next
      }
      if (store_trees) trees_list[[g]][[train_s]] <- fit$tree

      leaf_all <- fit$leaf_all
      s_cell <- as.integer(df_cell[[sample_col]])
      score_cell <- as.numeric(df_cell[[scores_col]])
      w_cell <- if (!is.null(weight_col)) as.numeric(df_cell[[weight_col]]) else rep(1.0, nrow(df_cell))
      w_cell[!is.finite(w_cell)] <- 0
      cl_cell <- if (!is.null(cluster_col)) df_cell[[cluster_col]] else NULL

      dt0 <- data.table::data.table(
        leaf   = leaf_all,
        sample = s_cell,
        score  = score_cell,
        w      = w_cell
      )
      if (!is.null(cluster_col)) {
        dt0[, cl := as.integer(factor(cl_cell, exclude = NULL))]
      } else {
        dt0[, cl := .I]
      }

      dt_train <- dt0[sample == train_s & !is.na(leaf)]
      train_leaf <- dt_train[, {
        G_here <- data.table::uniqueN(cl)
        N_here <- .N
        if (G_here < minsize) {
          list(G = G_here, N = N_here, coef = NA_real_, stderr = NA_real_, t = NA_real_)
        } else {
          o <- mean_test_crv1(score, w, cl)
          list(G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t)
        }
      }, by = leaf]

      train_leaf[, `:=`(train = TRUE, sample = train_s)]
      train_leaf[, p.raw := stats::pnorm(t)]

      prom <- select_leaves_local(train_leaf[is.finite(t), .(leaf, t)])

      if (!is.null(key_dt) && ncol(key_dt) > 0L) {
        for (cc in names(key_dt)) train_leaf[, (cc) := key_dt[[cc]][1]]
      }

      dirs[[train_s]] <- list(
        idx = idx,
        key_dt = key_dt,
        leaf_all = leaf_all,
        train_leaf = train_leaf,
        prom = as.integer(prom),
        train_s = train_s,
        est_s = est_s
      )
    }

    cell_objs[[g]] <- dirs

    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("CART_test: processed %d/%d margin cells", g, n_cells))
    }
  }

  global_out <- data.table::rbindlist(global_list, use.names = TRUE, fill = TRUE)

  # ---------------- collect candidate leaves ----------------
  leaf_rows <- list()
  lr <- 0L

  for (g in seq_len(n_cells)) {
    for (train_s in c(1L, 2L)) {
      obj <- cell_objs[[g]][[train_s]]
      if (is.null(obj)) next

      tl <- data.table::copy(obj$train_leaf)
      tl[, `:=`(
        cell_id = g,
        train_s = train_s,
        est_s = obj$est_s,
        selected_local = as.integer(leaf %in% obj$prom)
      )]
      leaf_rows[[lr <- lr + 1L]] <- tl
    }
  }

  leaf_tbl <- data.table::rbindlist(leaf_rows, use.names = TRUE, fill = TRUE)

  # Restrict candidate set to locally preselected ("promising") leaves.
  leaf_tbl_use <- leaf_tbl[selected_local == 1L & is.finite(t)]

  if (nrow(leaf_tbl_use) == 0L) {
    stop("Testing subset is empty (after local preselection).")
  }

  # group across selected dimensions; keep non-selected dimensions separate
  grp_cols <- c(cell_cols, if (!select_sample) "train_s" else character())
  grp_cols <- unique(grp_cols)

  if (length(grp_cols) == 0L) {
    sel_obj <- choose_group_selection(leaf_tbl_use[, .(cell_id, train_s, leaf, t)])
    sel_tbl <- data.table::data.table(
      best_cell = sel_obj$best_cell,
      best_train_s = sel_obj$best_train_s
    )
    sel_tbl[, sel_leaves := list(sel_obj$sel_leaves)]
  } else {
    sel_tbl <- leaf_tbl_use[, {
      z <- choose_group_selection(.SD[, .(cell_id, train_s, leaf, t)])
      data.table::data.table(
        best_cell = z$best_cell,
        best_train_s = z$best_train_s,
        sel_leaves = list(z$sel_leaves)
      )
    }, by = grp_cols]
  }

  # ---------------- mark relevant train rows + build testing jobs ----------------
  # ---------------- mark relevant train rows + build testing jobs ----------------
  train_rows <- list()
  tr_k <- 0L
  test_jobs <- list()
  tj_k <- 0L

  for (g in seq_len(n_cells)) {
    for (train_s in c(1L, 2L)) {
      obj <- cell_objs[[g]][[train_s]]
      if (is.null(obj)) next

      key_row <- data.table::data.table()
      if (length(cell_cols) > 0L && !is.null(obj$key_dt) && ncol(obj$key_dt) > 0L) {
        for (cc in cell_cols) key_row[, (cc) := obj$key_dt[[cc]][1]]
      }
      if (!select_sample) key_row[, train_s := train_s]

      if (nrow(sel_tbl) == 1L && length(grp_cols) == 0L) {
        sel_here <- sel_tbl
      } else {
        sel_here <- merge(key_row, sel_tbl, by = grp_cols, all.x = TRUE)
      }

      sel_leaves <- integer()
      chosen <- FALSE
      if (nrow(sel_here) == 1L &&
          is.finite(sel_here$best_cell[1]) &&
          is.finite(sel_here$best_train_s[1]) &&
          g == sel_here$best_cell[1] &&
          train_s == sel_here$best_train_s[1]) {
        sel_leaves <- as.integer(unlist(sel_here$sel_leaves[[1]], use.names = FALSE))
        chosen <- length(sel_leaves) > 0L
      }

      tl <- data.table::copy(obj$train_leaf)
      tl[, relevant := as.integer(leaf %in% sel_leaves)]
      train_rows[[tr_k <- tr_k + 1L]] <- tl

      if (chosen) {
        for (lf in sel_leaves) {
          test_jobs[[tj_k <- tj_k + 1L]] <- list(
            cell_id = g,
            train_s = train_s,
            est_s = obj$est_s,
            leaf = as.integer(lf)
          )
        }
      }
    }
  }

  train_out <- data.table::rbindlist(train_rows, use.names = TRUE, fill = TRUE)

  if (length(test_jobs) == 0L) stop("Testing subset is empty (after selection).")

  # ---------------- testing: one test per selected leaf ----------------
  test_rows <- list()
  tt_k <- 0L

  for (j in seq_along(test_jobs)) {
    job <- test_jobs[[j]]
    obj <- cell_objs[[job$cell_id]][[job$train_s]]

    idx <- obj$idx
    leaf_all <- obj$leaf_all
    est_s <- obj$est_s

    keep_local <- which(
      svec[idx] == est_s &
        !is.na(leaf_all) &
        leaf_all == job$leaf
    )
    if (length(keep_local) == 0L) next

    idx_keep <- idx[keep_local]

    y <- as.numeric(data[[scores_col]])[idx_keep]
    w <- if (!is.null(weight_col)) as.numeric(data[[weight_col]])[idx_keep] else rep(1.0, length(idx_keep))
    w[!is.finite(w)] <- 0
    cl <- if (!is.null(cluster_col)) data[[cluster_col]][idx_keep] else NULL

    o <- mean_test_crv1(y, w, cl)

    rr <- data.table::data.table(
      train = FALSE,
      sample = est_s,
      leaf = as.integer(job$leaf),
      G = o$G,
      N = o$N,
      coef = o$coef,
      stderr = o$se,
      t = o$t,
      p.raw = stats::pnorm(o$t),
      relevant = NA_integer_
    )

    if (!is.null(obj$key_dt) && ncol(obj$key_dt) > 0L) {
      for (cc in names(obj$key_dt)) rr[, (cc) := obj$key_dt[[cc]][1]]
    }

    test_rows[[tt_k <- tt_k + 1L]] <- rr
  }

  if (length(test_rows) == 0L) stop("Testing subset is empty (after selection).")
  test_out <- data.table::rbindlist(test_rows, use.names = TRUE, fill = TRUE)

  # combine results
  results_out <- data.table::rbindlist(list(train_out, test_out), use.names = TRUE, fill = TRUE)
  if (!("relevant" %chin% names(results_out))) results_out[, relevant := NA_integer_]
  results_out[train == FALSE, relevant := NA_integer_]
  data.table::setorder(results_out, train)

  out <- list(
    results = results_out,
    global  = global_out
  )
  if (store_trees) out$trees <- trees_list
  out
}

weighted_quantile <- function(x, w = NULL, probs = c(0.25, 0.5, 0.75),
                              na.rm = TRUE,
                              interpolation = c("linear", "left", "right", "midpoint")) {
  interpolation <- match.arg(interpolation)

  # default: no weights -> equal weights
  if (is.null(w)) w <- rep(1, length(x))

  # basic checks
  if (length(x) != length(w)) stop("x and w must have the same length.")
  if (na.rm) {
    keep <- is.finite(x) & is.finite(w)
    x <- x[keep]; w <- w[keep]
  }
  if (any(w < 0)) stop("weights must be nonnegative.")
  if (!length(x)) return(rep(NA_real_, length(probs)))
  s <- sum(w)
  if (s == 0) return(rep(NA_real_, length(probs)))

  # sort by x
  o <- order(x, na.last = NA)
  x <- x[o]; w <- w[o]

  # cumulative weight share in [0,1]
  cw <- cumsum(w) / s

  # helper to get quantile for a single p
  get_q <- function(p) {
    if (p <= 0) return(x[1L])
    if (p >= 1) return(x[length(x)])

    # first index where cw >= p
    k <- which.max(cw >= p)  # same as min(which(cw >= p)) but faster
    if (cw[k] == p) {
      # exact hit: handle by interpolation rule
      if (interpolation == "left")      return(x[k])
      if (interpolation == "right")     return(x[k])
      if (interpolation == "midpoint") {
        j <- k
        while (j > 1L && cw[j - 1L] == cw[k]) j <- j - 1L
        i <- k
        while (i < length(x) && cw[i + 1L] == cw[k]) i <- i + 1L
        return((x[j] + x[i]) / 2)
      }
      return(x[k]) # linear or default on exact hit
    } else {
      # cw[k-1] < p < cw[k]
      if (k == 1L) return(x[1L])
      if (interpolation == "left")      return(x[k - 1L])
      if (interpolation == "right")     return(x[k])
      if (interpolation == "midpoint")  return((x[k - 1L] + x[k]) / 2)

      # linear interpolation in cumulative weight space
      w_below <- cw[k - 1L]
      w_above <- cw[k]
      t <- (p - w_below) / (w_above - w_below)
      return(x[k - 1L] + t * (x[k] - x[k - 1L]))
    }
  }

  # vectorize over probs
  res <- vapply(probs, get_q, numeric(1))
  names(res) <- paste0(probs)
  res
}

w_sd <- function(x, w = NULL, na.rm = TRUE, unbiased = FALSE) {
  if (is.null(w)) return(sd(x, na.rm = na.rm))

  if (na.rm) {
    keep <- is.finite(x) & is.finite(w)
    x <- x[keep]; w <- w[keep]
  }
  if (!length(x)) return(NA_real_)

  mu <- weighted.mean(x, w)
  v  <- sum(w * (x - mu)^2) / sum(w)

  if (unbiased && sum(w) > 1)
    v <- v * sum(w) / (sum(w) - 1)

  sqrt(v)
}




######### MAIN FOREST HELPERS ###############

# ========= Helper: Randomize into K folds =========
make_group_folds <- function(DT,
                             K = 10L,
                             cluster_name = NULL,
                             fold_col = "fold",
                             by_col = NULL,        # NULL = global; else within groups
                             diag_prefix = "cf_",
                             verbose = TRUE) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.numeric(K) && K >= 2)
  stopifnot(is.character(fold_col), length(fold_col) == 1L)

  if (!is.null(by_col)) {
    stopifnot(is.character(by_col), length(by_col) == 1L, by_col %chin% names(DT))
  }
  if (!is.null(cluster_name)) {
    stopifnot(is.character(cluster_name), length(cluster_name) == 1L, cluster_name %chin% names(DT))
  }

  has_diag <- !is.null(diag_prefix)
  if (has_diag) {
    k_col <- paste0(diag_prefix, "K")
    g_col <- paste0(diag_prefix, "G")
    n_col <- paste0(diag_prefix, "n_ok")
  }

  # helper: create folds for one group
  make_one <- function(n, cl = NULL) {
    ok <- rep(TRUE, n)
    n_ok <- n

    fold <- rep(NA_integer_, n)

    if (!is.null(cl)) {
      ok <- ok & !is.na(cl)
      n_ok <- sum(ok)

      if (n_ok == 0L) {
        if (has_diag) return(list(fold, 0L, 0L, 0L))
        return(fold)
      }

      cl_ok <- cl[ok]
      ucl <- unique(cl_ok)
      G <- length(ucl)
      Kg <- as.integer(min(K, G))

      if (Kg < 2L) {
        fold[ok] <- 1L
        if (has_diag) return(list(fold, Kg, G, n_ok))
        return(fold)
      }

      fold_cl <- sample(rep_len(seq_len(Kg), G))
      fold[ok] <- as.integer(fold_cl[match(cl_ok, ucl)])

      if (has_diag) return(list(fold, Kg, G, n_ok))
      return(fold)
    }

    # no clustering
    if (n_ok == 0L) {
      if (has_diag) return(list(fold, 0L, NA_integer_, 0L))
      return(fold)
    }

    Kg <- as.integer(min(K, n_ok))
    if (Kg < 2L) {
      fold[ok] <- 1L
      if (has_diag) return(list(fold, Kg, NA_integer_, n_ok))
      return(fold)
    }

    fold[ok] <- as.integer(sample(rep_len(seq_len(Kg), n_ok)))

    if (has_diag) return(list(fold, Kg, NA_integer_, n_ok))
    fold
  }

  if (is.null(by_col)) {
    # ---------- GLOBAL SPLIT ----------
    n <- nrow(DT)
    cl <- if (!is.null(cluster_name)) DT[[cluster_name]] else NULL

    if (has_diag) {
      DT[, c(fold_col, k_col, g_col, n_col) := make_one(n, cl)]
    } else {
      DT[, (fold_col) := make_one(n, cl)]
    }

    if (verbose) {
      message(sprintf(
        "Cross-fitting folds: created '%s' with up to %d folds (global).",
        fold_col, K
      ))
    }

  } else {
    # ---------- WITHIN-GROUP SPLIT ----------
    if (has_diag) {
      DT[, c(fold_col, k_col, g_col, n_col) := {
        cl <- if (!is.null(cluster_name)) get(cluster_name) else NULL
        make_one(.N, cl)
      }, by = by_col]
    } else {
      DT[, (fold_col) := {
        cl <- if (!is.null(cluster_name)) get(cluster_name) else NULL
        make_one(.N, cl)
      }, by = by_col]
    }

    if (verbose) {
      message(sprintf(
        "Cross-fitting folds: created '%s' with up to %d folds within '%s'.",
        fold_col, K, by_col
      ))
    }
  }

  invisible(DT)
}





# ========= Helper: Predict nuissance using regression-forest - cross-fit across sample parts =========
crossfit_hat <- function(DT,
                         i = NULL,
                         y_name,
                         x_names,
                         margins = NULL,
                         weight_name = NULL,
                         folds = NULL,              # NULL => OOB if mode = "within"
                         forest_opts = list(),
                         hat_suffix = ".hat",
                         mode = c("within", "across"),
                         verbose = FALSE) {

  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.character(y_name), length(y_name) == 1L)
  stopifnot(is.character(x_names), length(x_names) >= 1L)
  stopifnot("sample" %chin% names(DT))
  stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

  mode <- match.arg(mode)

  if (!is.null(folds)) {
    stopifnot(is.character(folds), length(folds) == 1L, folds %chin% names(DT))
  }

  if (is.null(i)) i <- DT[, .I]
  stopifnot(is.integer(i) || is.numeric(i))
  i <- as.integer(i)
  n_i <- length(i)
  if (n_i == 0L) return(invisible(DT))

  if (is.null(margins)) margins <- character()
  grp_cols <- as.character(margins)

  missing_cols <- setdiff(
    c("sample", y_name, x_names, grp_cols,
      if (!is.null(weight_name)) weight_name,
      if (!is.null(folds)) folds),
    names(DT)
  )
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  hat_col <- paste0(y_name, hat_suffix)
  if (!(hat_col %chin% names(DT))) DT[, (hat_col) := NA_real_]

  X_all <- as.matrix(DT[i, ..x_names])

  rid_col <- ".__rid__"
  DT[i, (rid_col) := seq_len(n_i)]
  rid <- DT[[rid_col]]

  sample_all <- DT[["sample"]]
  y_all <- as.numeric(DT[[y_name]])
  w_all <- if (!is.null(weight_name)) as.numeric(DT[[weight_name]]) else NULL
  f_all <- if (!is.null(folds)) DT[[folds]] else NULL

  preds_buf <- rep(NA_real_, n_i)

  fit_predict_idx <- function(idx_tr, idx_te) {
    n_te <- length(idx_te)
    preds <- rep(NA_real_, n_te)
    if (n_te == 0L) return(preds)

    n_tr <- length(idx_tr)
    if (n_tr == 0L) return(preds)

    y_tr <- y_all[idx_tr]
    w_tr <- if (is.null(w_all)) NULL else w_all[idx_tr]

    if (n_tr < 2L) {
      mu <- if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
      preds[] <- mu
      return(preds)
    }

    y_fin <- y_tr[is.finite(y_tr)]
    if (length(y_fin) == 0L) return(preds)
    if (!any(y_fin != y_fin[1L])) {
      preds[] <- y_fin[1L]
      return(preds)
    }

    X_tr <- X_all[rid[idx_tr], , drop = FALSE]
    rf <- do.call(
      grf::regression_forest,
      c(list(
        X = X_tr,
        Y = y_tr,
        sample.weights = w_tr,
        compute.oob.predictions = FALSE
      ), forest_opts)
    )

    X_te <- X_all[rid[idx_te], , drop = FALSE]
    preds[] <- as.numeric(predict(rf, X_te)$predictions)
    preds
  }

  within_pred_idx <- function(idx_s) {
    n <- length(idx_s)
    if (n == 0L) return(numeric())

    f <- f_all[idx_s]
    if (data.table::uniqueN(f) < 2L) {
      y_s <- y_all[idx_s]
      w_s <- if (is.null(w_all)) NULL else w_all[idx_s]

      if (n < 2L) {
        mu <- if (is.null(w_s)) mean(y_s) else stats::weighted.mean(y_s, w_s)
        return(rep(mu, n))
      }

      y_fin <- y_s[is.finite(y_s)]
      if (length(y_fin) == 0L) return(rep(NA_real_, n))
      if (!any(y_fin != y_fin[1L])) return(rep(y_fin[1L], n))

      fit <- do.call(
        grf::regression_forest,
        c(list(
          X = X_all[rid[idx_s], , drop = FALSE],
          Y = y_s,
          sample.weights = w_s,
          compute.oob.predictions = TRUE
        ), forest_opts)
      )
      return(as.numeric(predict(fit)$predictions))
    }

    pos_by_fold <- split(seq_len(n), f)

    p <- rep(NA_real_, n)
    tr_mask <- rep.int(TRUE, n)

    for (k in names(pos_by_fold)) {
      te_pos <- pos_by_fold[[k]]
      tr_mask[te_pos] <- FALSE
      tr_pos <- which(tr_mask)
      tr_mask[te_pos] <- TRUE

      idx_tr <- idx_s[tr_pos]
      idx_te <- idx_s[te_pos]

      if (length(idx_tr) < 2L) {
        y_tr <- y_all[idx_tr]
        w_tr <- if (is.null(w_all)) NULL else w_all[idx_tr]
        mu <- if (length(y_tr) == 0L) NA_real_ else {
          if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
        }
        p[te_pos] <- mu
      } else {
        p[te_pos] <- fit_predict_idx(idx_tr, idx_te)
      }
    }
    p
  }

  within_oob_idx <- function(idx_s) {
    n <- length(idx_s)
    p <- rep(NA_real_, n)
    if (n == 0L) return(p)

    y_s <- y_all[idx_s]
    w_s <- if (is.null(w_all)) NULL else w_all[idx_s]

    if (n < 2L) {
      mu <- if (is.null(w_s)) mean(y_s) else stats::weighted.mean(y_s, w_s)
      p[] <- mu
      return(p)
    }

    y_fin <- y_s[is.finite(y_s)]
    if (length(y_fin) == 0L) return(p)
    if (!any(y_fin != y_fin[1L])) {
      p[] <- y_fin[1L]
      return(p)
    }

    fit <- do.call(
      grf::regression_forest,
      c(list(
        X = X_all[rid[idx_s], , drop = FALSE],
        Y = y_s,
        sample.weights = w_s,
        compute.oob.predictions = TRUE
      ), forest_opts)
    )

    p[] <- as.numeric(predict(fit)$predictions)
    p
  }

  idx_all <- i
  if (length(grp_cols) == 0L) {
    group_list <- list(idx_all)
  } else {
    group_dt <- DT[i, ..grp_cols]
    gid <- data.table::frankv(group_dt, ties.method = "dense")
    group_list <- split(idx_all, gid)
  }

  for (g in seq_along(group_list)) {
    idx_g <- group_list[[g]]
    if (length(idx_g) == 0L) next

    idx1 <- idx_g[sample_all[idx_g] == 1L]
    idx2 <- idx_g[sample_all[idx_g] == 2L]

    if (mode == "across") {
      if (length(idx1) > 1L && length(idx2) > 0L) {
        preds_buf[rid[idx2]] <- fit_predict_idx(idx1, idx2)
      } else if (length(idx1) == 1L && length(idx2) > 0L) {
        preds_buf[rid[idx2]] <- y_all[idx1]
      }

      if (length(idx2) > 1L && length(idx1) > 0L) {
        preds_buf[rid[idx1]] <- fit_predict_idx(idx2, idx1)
      } else if (length(idx2) == 1L && length(idx1) > 0L) {
        preds_buf[rid[idx1]] <- y_all[idx2]
      }

    } else {
      if (is.null(folds)) {
        if (length(idx1) > 0L) preds_buf[rid[idx1]] <- within_oob_idx(idx1)
        if (length(idx2) > 0L) preds_buf[rid[idx2]] <- within_oob_idx(idx2)
      } else {
        if (length(idx1) > 0L) preds_buf[rid[idx1]] <- within_pred_idx(idx1)
        if (length(idx2) > 0L) preds_buf[rid[idx2]] <- within_pred_idx(idx2)
      }
    }

    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("crossfit_hat: processed %d/%d groups", g, length(group_list)))
    }
  }

  DT[i, (hat_col) := preds_buf]
  DT[, (rid_col) := NULL]

  invisible(DT)
}

# ========= Helper: Estimate causal/regression/instrumental forests =========
fit_models <- function(DT,
                       i = NULL,
                       forest_type = c("causal", "regression", "instrumental"),
                       y_name,
                       x_names,
                       margins = NULL,
                       w_name = NULL,
                       z_name = NULL,
                       folds = NULL,              # NULL => within-sample OOB for pred
                       weight_name = NULL,
                       cluster_name = NULL,
                       forest_opts = NULL,
                       aipw.clip = 1e-3,
                       shrink = FALSE,
                       verbose = FALSE) {

  stopifnot(data.table::is.data.table(DT))
  forest_type <- match.arg(forest_type)

  stopifnot("sample" %chin% names(DT))
  stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

  if (!is.null(folds)) stopifnot(is.character(folds), length(folds) == 1L, folds %chin% names(DT))
  if (!is.null(aipw.clip)) {
    stopifnot(is.numeric(aipw.clip), length(aipw.clip) == 1L, is.finite(aipw.clip), aipw.clip > 0, aipw.clip < 1)
  }
  stopifnot(is.logical(shrink), length(shrink) == 1L, !is.na(shrink))

  if (is.null(i)) i <- DT[, .I]
  i <- as.integer(i)
  n_i <- length(i)
  if (n_i == 0L) return(invisible(DT))

  if (is.null(margins)) margins <- character()
  grp_cols <- as.character(margins)

  if (is.null(forest_opts)) forest_opts <- list()

  if (!("pred"   %chin% names(DT))) DT[, pred   := NA_real_]
  if (!("pred_o" %chin% names(DT))) DT[, pred_o := NA_real_]
  if (!("scores" %chin% names(DT))) DT[, scores := NA_real_]
  if (shrink) {
    if (!("pred_var"   %chin% names(DT))) DT[, pred_var   := NA_real_]
    if (!("pred_o_var" %chin% names(DT))) DT[, pred_o_var := NA_real_]
  }

  y_hat <- paste0(y_name, ".hat")
  w_hat <- if (!is.null(w_name)) paste0(w_name, ".hat") else NULL
  z_hat <- if (!is.null(z_name)) paste0(z_name, ".hat") else NULL

  if (forest_type %in% c("causal", "instrumental")) {
    stopifnot(!is.null(w_name), w_name %chin% names(DT))
    stopifnot(y_hat %chin% names(DT), w_hat %chin% names(DT))
  }
  if (forest_type == "instrumental") {
    stopifnot(!is.null(z_name), z_name %chin% names(DT))
    stopifnot(z_hat %chin% names(DT))
  }

  X_all <- as.matrix(DT[i, ..x_names])

  rid_col <- ".__rid__"
  DT[i, (rid_col) := seq_len(n_i)]
  rid <- DT[[rid_col]]

  sample_all <- DT[["sample"]]
  folds_all  <- if (!is.null(folds)) DT[[folds]] else NULL
  wgt_all    <- if (!is.null(weight_name)) as.numeric(DT[[weight_name]]) else NULL
  cl_all     <- if (!is.null(cluster_name)) DT[[cluster_name]] else NULL

  y_all <- as.numeric(DT[[y_name]])

  if (forest_type %in% c("causal", "instrumental")) {
    w_all    <- as.numeric(DT[[w_name]])
    yhat_all <- as.numeric(DT[[y_hat]])
    what_all <- as.numeric(DT[[w_hat]])
  }
  if (forest_type == "instrumental") {
    z_all    <- as.numeric(DT[[z_name]])
    zhat_all <- as.numeric(DT[[z_hat]])
  }

  build_forest_idx <- function(type, idx, compute.oob.predictions = FALSE) {
    X <- X_all[rid[idx], , drop = FALSE]
    Y <- y_all[idx]
    sw <- if (is.null(wgt_all)) NULL else wgt_all[idx]
    cl <- if (is.null(cl_all)) NULL else cl_all[idx]

    if (type == "regression") {
      do.call(grf::regression_forest,
              c(list(X = X, Y = Y, sample.weights = sw, clusters = cl,
                     compute.oob.predictions = compute.oob.predictions), forest_opts))
    } else if (type == "causal") {
      W <- w_all[idx]
      do.call(grf::causal_forest,
              c(list(X = X, Y = Y, W = W,
                     Y.hat = yhat_all[idx],
                     W.hat = what_all[idx],
                     sample.weights = sw,
                     clusters = cl,
                     compute.oob.predictions = compute.oob.predictions), forest_opts))
    } else {
      W <- w_all[idx]
      Z <- z_all[idx]
      do.call(grf::instrumental_forest,
              c(list(X = X, Y = Y, W = W, Z = Z,
                     Y.hat = yhat_all[idx],
                     W.hat = what_all[idx],
                     Z.hat = zhat_all[idx],
                     sample.weights = sw,
                     clusters = cl,
                     compute.oob.predictions = compute.oob.predictions), forest_opts))
    }
  }

  predict_with_optional_var <- function(fit, X) {
    if (!shrink) {
      p <- predict(fit, X)
      return(list(pred = as.numeric(p$predictions), var = NULL))
    }

    p <- predict(fit, X, estimate.variance = TRUE)
    v <- if (!is.null(p$variance.estimates)) {
      pmax(as.numeric(p$variance.estimates), 0)
    } else {
      rep(NA_real_, nrow(X))
    }

    list(pred = as.numeric(p$predictions), var = v)
  }

  predict_oob_with_optional_var <- function(fit, n) {
    if (!shrink) {
      p <- predict(fit)
      return(list(pred = as.numeric(p$predictions), var = NULL))
    }

    p <- predict(fit, estimate.variance = TRUE)
    v <- if (!is.null(p$variance.estimates)) {
      pmax(as.numeric(p$variance.estimates), 0)
    } else {
      rep(NA_real_, n)
    }

    list(pred = as.numeric(p$predictions), var = v)
  }

  within_pred_idx <- function(idx_s) {
    n <- length(idx_s)
    if (n == 0L) {
      return(list(pred = numeric(), var = if (shrink) numeric() else NULL))
    }

    Xs <- X_all[rid[idx_s], , drop = FALSE]

    f <- folds_all[idx_s]
    if (data.table::uniqueN(f) < 2L) {
      fit <- build_forest_idx(forest_type, idx_s, compute.oob.predictions = TRUE)
      return(predict_oob_with_optional_var(fit, n))
    }

    pos_by_fold <- split(seq_len(n), f)

    p  <- rep(NA_real_, n)
    vv <- if (shrink) rep(NA_real_, n) else NULL
    tr_mask <- rep.int(TRUE, n)

    for (k in names(pos_by_fold)) {
      te_pos <- pos_by_fold[[k]]
      tr_mask[te_pos] <- FALSE
      tr_pos <- which(tr_mask)
      tr_mask[te_pos] <- TRUE

      if (length(tr_pos) < 2L) {
        idx_tr <- idx_s[tr_pos]
        y_tr <- y_all[idx_tr]
        w_tr <- if (is.null(wgt_all)) NULL else wgt_all[idx_tr]
        mu <- if (length(y_tr) == 0L) {
          NA_real_
        } else {
          if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
        }
        p[te_pos] <- mu
        if (shrink) vv[te_pos] <- NA_real_
      } else {
        fit_k <- build_forest_idx(forest_type, idx_s[tr_pos], compute.oob.predictions = FALSE)
        out_k <- predict_with_optional_var(fit_k, Xs[te_pos, , drop = FALSE])
        p[te_pos] <- out_k$pred
        if (shrink) vv[te_pos] <- out_k$var
      }
    }

    list(pred = p, var = vv)
  }

  within_oob_idx <- function(idx_s) {
    n <- length(idx_s)
    if (n == 0L) {
      return(list(pred = numeric(), var = if (shrink) numeric() else NULL))
    }
    fit <- build_forest_idx(forest_type, idx_s, compute.oob.predictions = TRUE)
    predict_oob_with_optional_var(fit, n)
  }

  compute_aipw_vec <- function(idx_s, tau_vec) {
    n <- length(idx_s)
    if (n == 0L) return(numeric())

    if (forest_type == "regression") return(y_all[idx_s])
    if (forest_type != "causal")     return(rep(NA_real_, n))

    Y <- y_all[idx_s]
    m <- yhat_all[idx_s]
    e <- what_all[idx_s]
    W <- w_all[idx_s]

    W <- as.numeric(W > 0.5)

    if (any(e < 0 | e > 1, na.rm = TRUE)) e <- plogis(e)
    if (!is.null(aipw.clip)) e <- pmin(pmax(e, aipw.clip), 1 - aipw.clip)

    tau <- as.numeric(tau_vec)
    stopifnot(length(tau) == n)

    m1 <- m + (1 - e) * tau
    m0 <- m - e * tau

    tau + (W / e) * (Y - m1) - ((1 - W) / (1 - e)) * (Y - m0)
  }

  idx_all <- i
  if (length(grp_cols) == 0L) {
    group_list <- list(idx_all)
  } else {
    group_dt <- DT[i, ..grp_cols]
    gid <- data.table::frankv(group_dt, ties.method = "dense")
    group_list <- split(idx_all, gid)
  }

  pred_buf   <- rep(NA_real_, n_i)
  pred_o_buf <- rep(NA_real_, n_i)
  score_buf  <- rep(NA_real_, n_i)

  if (shrink) {
    pred_var_buf   <- rep(NA_real_, n_i)
    pred_o_var_buf <- rep(NA_real_, n_i)
  }

  for (g in seq_along(group_list)) {
    idx_g <- group_list[[g]]
    if (length(idx_g) == 0L) next

    idx1 <- idx_g[sample_all[idx_g] == 1L]
    idx2 <- idx_g[sample_all[idx_g] == 2L]

    fit1 <- if (length(idx1)) build_forest_idx(
      forest_type, idx1,
      compute.oob.predictions = is.null(folds_all)
    ) else NULL

    fit2 <- if (length(idx2)) build_forest_idx(
      forest_type, idx2,
      compute.oob.predictions = is.null(folds_all)
    ) else NULL

    if (is.null(folds_all)) {
      res1 <- if (length(idx1)) within_oob_idx(idx1) else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }
      res2 <- if (length(idx2)) within_oob_idx(idx2) else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }
    } else {
      res1 <- if (length(idx1)) within_pred_idx(idx1) else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }
      res2 <- if (length(idx2)) within_pred_idx(idx2) else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }
    }

    res1_o <- if (!is.null(fit2) && length(idx1)) {
      predict_with_optional_var(fit2, X_all[rid[idx1], , drop = FALSE])
    } else list(pred = numeric(), var = if (shrink) numeric() else NULL)

    res2_o <- if (!is.null(fit1) && length(idx2)) {
      predict_with_optional_var(fit1, X_all[rid[idx2], , drop = FALSE])
    } else list(pred = numeric(), var = if (shrink) numeric() else NULL)

    p1   <- res1$pred
    p2   <- res2$pred
    p1_o <- res1_o$pred
    p2_o <- res2_o$pred

    if (shrink) {
      v1   <- res1$var
      v2   <- res2$var
      v1_o <- res1_o$var
      v2_o <- res2_o$var
    }

    sc1 <- if (length(idx1)) compute_aipw_vec(idx1, p1) else numeric()
    sc2 <- if (length(idx2)) compute_aipw_vec(idx2, p2) else numeric()

    if (length(idx1)) {
      r1 <- rid[idx1]
      pred_buf[r1]   <- p1
      pred_o_buf[r1] <- p1_o
      score_buf[r1]  <- sc1
      if (shrink) {
        pred_var_buf[r1]   <- v1
        pred_o_var_buf[r1] <- v1_o
      }
    }
    if (length(idx2)) {
      r2 <- rid[idx2]
      pred_buf[r2]   <- p2
      pred_o_buf[r2] <- p2_o
      score_buf[r2]  <- sc2
      if (shrink) {
        pred_var_buf[r2]   <- v2
        pred_o_var_buf[r2] <- v2_o
      }
    }

    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("fit_models: processed %d/%d groups", g, length(group_list)))
    }
  }

  if (shrink) {
    DT[i, `:=`(pred = pred_buf,
               pred_var = pred_var_buf,
               pred_o = pred_o_buf,
               pred_o_var = pred_o_var_buf,
               scores = score_buf)]
  } else {
    DT[i, `:=`(pred = pred_buf,
               pred_o = pred_o_buf,
               scores = score_buf)]
  }

  DT[, (rid_col) := NULL]

  invisible(DT)
}

####### FIND OPTIMAL CUTOFF AND TEST #############
forest_test <- function(
    data,
    cluster = NULL,
    weight  = NULL,
    sample  = "sample",
    pred    = "pred",
    pred_o  = "pred_o",
    scores  = "scores",
    x_names = NULL,
    minsize = 50L,
    margins = NULL,
    pool = NULL,
    select = NULL,
    gridpoints = NULL,
    store_grid = TRUE,
    verbose = FALSE
) {
  stopifnot(data.table::is.data.table(data))

  sample_col <- as.character(sample)
  pred_col   <- as.character(pred)
  pred_o_col <- as.character(pred_o)
  scores_col <- as.character(scores)
  weight_col <- if (is.null(weight)) NULL else as.character(weight)

  stopifnot(length(sample_col) == 1L, sample_col %in% names(data))
  stopifnot(length(pred_col)   == 1L, pred_col   %in% names(data))
  stopifnot(length(pred_o_col) == 1L, pred_o_col %in% names(data))
  stopifnot(length(scores_col) == 1L, scores_col %in% names(data))
  stopifnot(is.numeric(minsize), length(minsize) == 1L, is.finite(minsize), minsize >= 1L)

  if (!is.null(weight_col)) stopifnot(weight_col %in% names(data))

  if (!is.null(cluster)) {
    cluster_col <- as.character(cluster)
    stopifnot(length(cluster_col) == 1L, cluster_col %in% names(data))
  } else {
    cluster_col <- NULL
  }

  if (is.null(margins)) margins <- character()
  margins <- as.character(margins)
  if (length(margins) > 0L) stopifnot(all(margins %in% names(data)))

  if (!is.null(x_names)) {
    x_names <- as.character(x_names)
    stopifnot(length(x_names) >= 1L, all(x_names %in% names(data)))
  }

  allowed_dim <- unique(c(margins, sample_col))

  if (is.null(pool)) pool <- character()
  pool <- unique(as.character(pool))
  bad_pool <- setdiff(pool, allowed_dim)
  if (length(bad_pool) > 0L) {
    stop("pool contains invalid names: ", paste(bad_pool, collapse = ", "))
  }

  if (is.null(select)) select <- character()
  select <- unique(as.character(select))
  bad_select <- setdiff(select, allowed_dim)
  if (length(bad_select) > 0L) {
    stop("select contains invalid names: ", paste(bad_select, collapse = ", "))
  }

  overlap <- intersect(pool, select)
  if (length(overlap) > 0L) {
    stop("The following names appear in both pool and select: ",
         paste(overlap, collapse = ", "))
  }

  pool_sample   <- sample_col %in% pool
  select_sample <- sample_col %in% select
  if (pool_sample && select_sample) {
    stop("sample cannot appear in both pool and select.")
  }

  pool_margins   <- intersect(pool, margins)
  select_margins <- intersect(select, margins)

  cell_cols <- setdiff(margins, pool_margins)
  adjust_cols <- setdiff(cell_cols, select_margins)

  test_by <- c(
    adjust_cols,
    if (!pool_sample && !select_sample) sample_col else character()
  )

  do_shares <- (length(pool_margins) > 0L)

  search_grid_on <- !is.null(gridpoints)
  if (search_grid_on) {
    stopifnot(is.numeric(gridpoints), length(gridpoints) == 1L,
              is.finite(gridpoints), gridpoints >= 2)
    gridpoints <- as.integer(gridpoints)
  }

  n <- nrow(data)
  samp <- as.integer(data[[sample_col]])
  stopifnot(all(samp %in% c(1L, 2L)))

  predv   <- as.numeric(data[[pred_col]])
  predov  <- as.numeric(data[[pred_o_col]])
  scorev  <- as.numeric(data[[scores_col]])

  wv <- if (!is.null(weight_col)) as.numeric(data[[weight_col]]) else rep(1.0, n)
  wv[!is.finite(wv)] <- 0

  clv <- if (!is.null(cluster_col)) {
    as.integer(factor(data[[cluster_col]], exclude = NULL))
  } else {
    seq_len(n)
  }

  crv1_mean <- function(score, w, cl) {
    dt0 <- data.table::data.table(score = score, w = w, cl = cl)
    gb <- dt0[, .(U = sum(w * score), W = sum(w)), by = cl]
    G <- nrow(gb)
    U <- sum(gb$U); W <- sum(gb$W)
    theta <- U / W
    ug <- gb$U - theta * gb$W
    se <- sqrt((G / pmax.int(G - 1L, 1L)) * sum(ug^2)) / abs(W)
    if (G < 2L || !is.finite(se)) se <- NA_real_
    list(coef = theta, se = se, t = theta / se, G = G, N = length(score))
  }

  percentile_indices <- function(w, k) {
    n0 <- length(w)
    if (n0 == 0L) return(integer())
    if (n0 <= k) return(seq_len(n0))
    w2 <- w
    w2[!is.finite(w2)] <- 0
    w2[w2 < 0] <- 0
    totw <- sum(w2)
    if (!is.finite(totw) || totw <= 0) {
      return(unique(as.integer(round(seq(1, n0, length.out = k)))))
    }
    cw <- cumsum(w2)
    probs <- (seq_len(k)) / (k + 1)
    targets <- probs * cw[n0]
    idx <- vapply(targets, function(tt) which(cw >= tt)[1L], integer(1))
    unique(pmin(pmax(idx, 1L), n0))
  }

  w_mean <- function(x, w) {
    ok <- is.finite(x) & is.finite(w)
    if (!any(ok)) return(NA_real_)
    sum(w[ok] * x[ok]) / sum(w[ok])
  }

  w_sd <- function(x, w) {
    ok <- is.finite(x) & is.finite(w)
    if (sum(ok) < 2L) return(NA_real_)
    mu <- sum(w[ok] * x[ok]) / sum(w[ok])
    sqrt(sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok]))
  }

  cutoff_index <- function(pred_local, cutoff) {
    if (!is.finite(cutoff) || length(pred_local) == 0L) return(NA_integer_)
    j <- which(pred_local >= cutoff)[1L]
    if (length(j) == 0L) length(pred_local) else j
  }

  if (length(cell_cols) == 0L) {
    idx_list <- list(seq_len(n))
    keys_list <- list(NULL)
  } else {
    gid <- data.table::frankv(data[, ..cell_cols], ties.method = "dense")
    idx_list <- split(seq_len(n), gid)
    keys_list <- lapply(idx_list, function(idx) data[idx[1L], ..cell_cols])
  }

  run_one_cell_idx <- function(idx, key_dt = NULL, cell_id = NA_integer_) {
    s  <- samp[idx]
    pr <- predv[idx]
    po <- predov[idx]
    sc <- scorev[idx]
    w  <- wv[idx]
    cl <- clv[idx]

    dt <- data.table::data.table(
      sample = s,
      pred   = pr,
      pred_o = po,
      score  = sc,
      w      = w,
      cl     = cl,
      rid    = seq_along(idx)
    )

    data.table::setorder(dt, sample, pred)

    dt[, N := seq_len(.N), by = sample]
    dt[, `:=`(a = w * score, b = w)]
    dt[, `:=`(WgY = cumsum(a), Wg = cumsum(b)), by = .(sample, cl)]
    dt[, `:=`(SW = cumsum(b), SWY = cumsum(a)), by = sample]
    dt[, m := SWY / SW, by = sample]

    dt[, `:=`(
      dTA2 =  WgY^2 - (WgY - a)^2,
      dTB2 =  Wg^2  - (Wg  - b)^2,
      dTAB =  WgY * Wg - (WgY - a) * (Wg - b)
    )]
    dt[, `:=`(TA2 = cumsum(dTA2), TB2 = cumsum(dTB2), TAB = cumsum(dTAB)), by = sample]
    dt[, G := cumsum(!duplicated(cl)), by = sample]

    dt[, sumS2 := (TA2 - 2*m*TAB + (m^2)*TB2) / (SW^2)]
    dt[, se := sqrt((G / pmax.int(G - 1L, 1L)) * sumS2)]
    dt[G < 2L, se := NA_real_]
    dt[, t_stat := m / se]

    tau_tr <- dt[G >= minsize, .(tau_tr = min(pred, na.rm = TRUE)), by = sample]

    dt_est <- dt[, .(sample, pred_o, cl, w)]
    data.table::setorder(dt_est, sample, pred_o)
    dt_est[, G_est := cumsum(!duplicated(cl)), by = sample]
    tau_est <- dt_est[G_est >= minsize, .(tau_est = min(pred_o, na.rm = TRUE)), by = sample]

    if (nrow(tau_tr) < 2L || nrow(tau_est) < 2L) {
      msg <- "minsize too large: cannot reach minsize clusters in train or est order for both samples."
      if (!is.null(key_dt) && ncol(key_dt) > 0L) {
        msg <- paste0(
          msg, " Cell: ",
          paste(paste(names(key_dt), as.character(key_dt[1, ]), sep = "="), collapse = ", ")
        )
      }
      stop(msg)
    }

    tau_vec <- c(
      max(as.numeric(c(tau_tr[sample == 1L, tau_tr], tau_est[sample == 2L, tau_est]))),
      max(as.numeric(c(tau_tr[sample == 2L, tau_tr], tau_est[sample == 1L, tau_est])))
    )
    dt[, tau_constraint := data.table::fifelse(sample == 1L, tau_vec[1], tau_vec[2])]

    if (!search_grid_on) {
      dt[, cand_search := TRUE]
    } else {
      dt[, cand_search := {
        idxq <- percentile_indices(w, gridpoints)
        cand <- rep(FALSE, .N)
        cand[idxq] <- TRUE
        cand
      }, by = sample]
    }

    elig <- (dt$pred >= dt$tau_constraint) & (dt$G >= minsize) & (dt$cand_search)

    res <- dt[elig, .SD[which.min(t_stat)], by = sample,
              .SDcols = c("G", "N", "m", "se", "t_stat", "pred")]

    cut1 <- res[sample == 1L, pred]
    cut2 <- res[sample == 2L, pred]
    if (length(cut1) == 0L || !is.finite(cut1)) cut1 <- tau_vec[1]
    if (length(cut2) == 0L || !is.finite(cut2)) cut2 <- tau_vec[2]

    grid_out <- NULL
    if (store_grid) {
      k_store <- 100L
      grid_out <- dt[, {
        idx_q <- percentile_indices(w, k_store)
        cut_s <- if (.BY$sample == 1L) cut1 else cut2
        idx_c <- cutoff_index(pred, cut_s)
        idx_all <- unique(c(idx_q, idx_c))
        idx_all <- idx_all[is.finite(idx_all)]
        idx_all <- sort(idx_all)
        out <- data.table::data.table(tau = pred[idx_all], t = t_stat[idx_all], chosen = FALSE)
        if (is.finite(idx_c) && nrow(out) > 0L) {
          j <- match(pred[idx_c], out$tau)
          if (!is.na(j)) out[j, chosen := TRUE]
        }
        out
      }, by = "sample"]

      if (!is.null(key_dt) && ncol(key_dt) > 0L) {
        for (cc in names(key_dt)) grid_out[, (cc) := key_dt[[cc]][1]]
        data.table::setcolorder(grid_out, c(names(key_dt), "sample", "t", "tau", "chosen"))
      }
    }

    idx_test_by_sample <- list(
      `1` = idx[dt$rid[dt$sample == 2L & dt$pred_o <= cut1]],
      `2` = idx[dt$rid[dt$sample == 1L & dt$pred_o <= cut2]]
    )

    train_out <- data.table::data.table(
      cell_id = cell_id,
      train = TRUE,
      relevant = 0L,
      sample = res$sample,
      G = res$G,
      N = res$N,
      coef = res$m,
      stderr = res$se,
      t = res$t_stat,
      tau_cutoff = res$pred,
      p.raw = stats::pnorm(res$t_stat)
    )

    if (!is.null(key_dt) && ncol(key_dt) > 0L) {
      for (cc in names(key_dt)) train_out[, (cc) := key_dt[[cc]][1]]
      data.table::setcolorder(train_out, c("cell_id", names(key_dt), "train", "relevant", "sample"))
    }

    list(train = train_out, grid = grid_out, idx_test = idx_test_by_sample)
  }

  n_cells <- length(idx_list)
  train_list <- vector("list", n_cells)
  grid_list  <- if (store_grid) vector("list", n_cells) else NULL
  idx_test_list <- vector("list", n_cells)

  for (g in seq_len(n_cells)) {
    ans <- run_one_cell_idx(idx_list[[g]], keys_list[[g]], cell_id = g)
    train_list[[g]] <- ans$train
    idx_test_list[[g]] <- ans$idx_test
    if (store_grid) grid_list[[g]] <- ans$grid
    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("forest_test: processed %d/%d cells", g, n_cells))
    }
  }

  train_all <- data.table::rbindlist(train_list, use.names = TRUE, fill = TRUE)
  grid_out  <- if (store_grid) data.table::rbindlist(grid_list, use.names = TRUE, fill = TRUE) else NULL

  selected_train <- train_all

  has_selection <- (length(select_margins) > 0L) || select_sample

  if (has_selection) {
    if (pool_sample) {
      agg_by <- c(adjust_cols, select_margins, "cell_id")
      cand <- train_all[, .(sel_t = mean(t, na.rm = TRUE)), by = agg_by]
      choose_by <- adjust_cols
      id_by <- c(adjust_cols, select_margins, "cell_id")
      keep_ids <- cand[, .SD[which.min(sel_t)], by = choose_by][, ..id_by]

      selected_train <- merge(train_all, keep_ids,
                              by = c(adjust_cols, select_margins, "cell_id"))
    } else if (select_sample) {
      agg_by <- c(adjust_cols, select_margins, sample_col, "cell_id")
      cand <- train_all[, .(sel_t = mean(t, na.rm = TRUE)), by = agg_by]
      choose_by <- adjust_cols
      id_by <- c(adjust_cols, select_margins, sample_col, "cell_id")
      keep_ids <- cand[, .SD[which.min(sel_t)], by = choose_by][, ..id_by]

      selected_train <- merge(train_all, keep_ids,
                              by = c(adjust_cols, select_margins, sample_col, "cell_id"))
    } else {
      agg_by <- c(adjust_cols, select_margins, sample_col, "cell_id")
      cand <- train_all[, .(sel_t = mean(t, na.rm = TRUE)), by = agg_by]
      choose_by <- c(adjust_cols, sample_col)
      id_by <- c(adjust_cols, sample_col, select_margins, "cell_id")
      keep_ids <- cand[, .SD[which.min(sel_t)], by = choose_by][, ..id_by]

      selected_train <- merge(train_all, keep_ids,
                              by = c(adjust_cols, sample_col, select_margins, "cell_id"))
    }
  }

  train_out <- data.table::copy(train_all)
  win_keys <- unique(selected_train[, .(cell_id, sample)])
  train_out[win_keys, on = .(cell_id, sample), relevant := 1L]
  train_out[, cell_id := NULL]

  idx_test_all <- integer()
  if (nrow(selected_train) > 0L) {
    for (i in seq_len(nrow(selected_train))) {
      cid <- selected_train$cell_id[i]
      s_i <- as.character(selected_train[[sample_col]][i])
      idx_test_all <- c(idx_test_all, idx_test_list[[cid]][[s_i]])
    }
  }
  idx_test_all <- unique(idx_test_all)

  in_test <- rep(FALSE, n)
  if (length(idx_test_all) > 0L) in_test[idx_test_all] <- TRUE
  if (!any(in_test)) stop("Testing subset is empty.")

  # keep selected margin values for final test rows
  selected_test_keys <- NULL
  if (length(select_margins) > 0L) {
    selected_test_keys <- unique(selected_train[, c(adjust_cols, select_margins), with = FALSE])
  }

  dt_test <- data.table::data.table(
    score = scorev,
    w     = wv,
    cl    = clv
  )
  if (length(test_by) > 0L) dt_test[, (test_by) := data[, ..test_by]]
  dt_test <- dt_test[in_test]

  if (length(test_by) == 0L) {
    o <- crv1_mean(dt_test$score, dt_test$w, dt_test$cl)
    test_out <- data.table::data.table(
      train = FALSE,
      relevant = 1L,
      sample = NA_integer_,
      G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
      tau_cutoff = NA_real_,
      p.raw = stats::pnorm(o$t)
    )
    if (!is.null(selected_test_keys)) {
      for (cc in select_margins) test_out[, (cc) := selected_test_keys[[cc]][1L]]
    }
  } else {
    test_out <- dt_test[, {
      o <- crv1_mean(score, w, cl)
      data.table::data.table(
        train = FALSE,
        relevant = 1L,
        G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
        tau_cutoff = NA_real_,
        p.raw = stats::pnorm(o$t)
      )
    }, by = test_by]

    if (!(sample_col %in% names(test_out))) test_out[, (sample_col) := NA_integer_]

    if (!is.null(selected_test_keys)) {
      # selected margins are constant within each adjust_cols group
      merge_by <- adjust_cols
      if (length(merge_by) == 0L) {
        for (cc in select_margins) test_out[, (cc) := selected_test_keys[[cc]][1L]]
      } else {
        test_out <- merge(
          test_out,
          selected_test_keys,
          by = merge_by,
          all.x = TRUE,
          sort = FALSE
        )
      }
    }
  }

  results_out <- data.table::rbindlist(list(train_out, test_out), use.names = TRUE, fill = TRUE)

  global <- global_means_crv1(
    data = data,
    scorev = scorev,
    wv = wv,
    clv = clv,
    by_cols = test_by,
    sample_col = sample_col,
    crv1_mean_fun = crv1_mean
  )

  Xmeans <- Xmeans_all <- XSD <- NULL
  if (!is.null(x_names)) {
    cols_need <- unique(c(test_by, x_names))
    df_all <- data[, ..cols_need]
    df_tst <- data[in_test, ..cols_need]

    wmeans_dt <- function(df, w, by_cols) {
      DTtmp <- data.table::as.data.table(df)
      DTtmp[, w__ := w]
      if (length(by_cols) > 0L) {
        DTtmp[, lapply(.SD, function(x) w_mean(as.numeric(x), w__)),
              by = by_cols, .SDcols = x_names]
      } else {
        DTtmp[, lapply(.SD, function(x) w_mean(as.numeric(x), w__)),
              .SDcols = x_names]
      }
    }

    wsds_dt <- function(df, w, by_cols) {
      DTtmp <- data.table::as.data.table(df)
      DTtmp[, w__ := w]
      if (length(by_cols) > 0L) {
        DTtmp[, lapply(.SD, function(x) w_sd(as.numeric(x), w__)),
              by = by_cols, .SDcols = x_names]
      } else {
        DTtmp[, lapply(.SD, function(x) w_sd(as.numeric(x), w__)),
              .SDcols = x_names]
      }
    }

    Xmeans <- wmeans_dt(df_tst, wv[in_test], test_by)
    Xmeans_all <- wmeans_dt(df_all, wv, test_by)
    XSD <- wsds_dt(df_all, wv, test_by)
  }

  shares <- NULL
  if (do_shares) {
    denom_by <- c(
      adjust_cols,
      if (!pool_sample && !select_sample) sample_col else character()
    )
    cols_need <- unique(c(denom_by, pool_margins))

    DT_all <- data[, ..cols_need]
    DT_tst <- data[in_test, ..cols_need]

    make_share_combo <- function(DTsub, wsub) {
      DTsub <- data.table::as.data.table(DTsub)
      DTsub[, w := wsub]

      by_num <- unique(c(denom_by, pool_margins))
      num <- DTsub[, .(w_sum = sum(w)), by = by_num]

      if (length(denom_by) == 0L) {
        num[, den_w := sum(DTsub$w)]
      } else {
        den <- DTsub[, .(den_w = sum(w)), by = denom_by]
        num <- num[den, on = denom_by]
      }

      num[, share := ifelse(is.finite(den_w) & den_w > 0, w_sum / den_w, NA_real_)]
      num[, c("w_sum", "den_w") := NULL]
      num
    }

    sh     <- make_share_combo(DT_tst, wv[in_test])
    sh_all <- make_share_combo(DT_all, wv)

    data.table::setnames(sh,     "share", "share")
    data.table::setnames(sh_all, "share", "share_all")

    by_merge <- unique(c(denom_by, pool_margins))
    shares <- merge(sh, sh_all, by = by_merge, all = TRUE)
    shares[is.na(share),     share := 0.0]
    shares[is.na(share_all), share_all := 0.0]

    data.table::setcolorder(shares, c(denom_by, pool_margins, "share", "share_all"))
  }

  out <- list(
    results = results_out,
    grid    = grid_out,
    global  = global,
    Xmeans  = Xmeans,
    Xmeans_all = Xmeans_all,
    XSD = XSD
  )
  if (!is.null(shares)) out$shares <- shares
  if (!store_grid) out$grid <- NULL

  out
}


##USED by CART_test and forest_test to estimate global means.
global_means_crv1 <- function(
    data,
    scorev,
    wv,
    clv,
    by_cols = character(),     # character vector of grouping cols (e.g. test_by)
    sample_col = "sample",     # name of sample column (only used to fill NA if absent)
    crv1_mean_fun             # function(score, w, cl) -> list(coef,se,t,G,N)
) {
  stopifnot(data.table::is.data.table(data))
  stopifnot(is.numeric(scorev), is.numeric(wv), length(scorev) == length(wv), length(wv) == length(clv))
  stopifnot(is.character(by_cols))
  stopifnot(is.character(sample_col), length(sample_col) == 1L)
  stopifnot(is.function(crv1_mean_fun))

  # Build small table
  dt_all <- data.table::data.table(score = scorev, w = wv, cl = clv)

  # Attach grouping columns (standard eval; no .. confusion)
  by_cols <- unique(by_cols)
  if (length(by_cols) > 0L) {
    missing <- setdiff(by_cols, names(data))
    if (length(missing)) stop("global_means_crv1: missing grouping columns in data: ",
                              paste(missing, collapse = ", "))
    dt_all[, (by_cols) := data[, ..by_cols]]
  }

  # Compute
  if (length(by_cols) == 0L) {
    o <- crv1_mean_fun(dt_all$score, dt_all$w, dt_all$cl)
    global_dt <- data.table::data.table(
      train = FALSE,
      sample = NA_integer_,
      G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
      tau_cutoff = NA_real_,
      p.raw = stats::pnorm(o$t)
    )
  } else {
    global_dt <- dt_all[, {
      o <- crv1_mean_fun(score, w, cl)
      data.table::data.table(
        train = FALSE,
        G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
        tau_cutoff = NA_real_,
        p.raw = stats::pnorm(o$t)
      )
    }, by = by_cols]

    # If sample_col isn't in the grouping, add it so schema matches your results tables
    if (!(sample_col %chin% names(global_dt))) global_dt[, (sample_col) := NA_integer_]
  }

  global_dt
}

#######ONE SIDED NONCOMPLIANCE HELPERS #######
test_one_sided_noncompliance <- function(data, D, Z, margins = character(0)) {
  dt <- data.table::as.data.table(data)

  stopifnot(is.character(D), length(D) == 1, D %in% names(dt))
  stopifnot(is.character(Z), length(Z) == 1, Z %in% names(dt))

  if (is.null(margins)) margins <- character(0)
  stopifnot(is.character(margins), all(margins %in% names(dt)))

  # binary check (strict, allow NA)
  if (!all(dt[[D]] %in% c(0L, 1L), na.rm = TRUE)) stop("D must be binary 0/1 (allowing NA).")
  if (!all(dt[[Z]] %in% c(0L, 1L), na.rm = TRUE)) stop("Z must be binary 0/1 (allowing NA).")

  core <- function(.SD) {
    z <- .SD[[Z]]
    d <- .SD[[D]]

    n_z0 <- sum(z == 0L)
    n_z1 <- sum(z == 1L)

    p_always <- if (n_z0 > 0L) mean(d[z == 0L] == 1L) else NA_real_  # P(D=1|Z=0)
    p_never  <- if (n_z1 > 0L) mean(d[z == 1L] == 0L) else NA_real_  # P(D=0|Z=1)

    only_always <- is.finite(p_always) && is.finite(p_never) && (p_always > 0) && (p_never == 0)
    only_never  <- is.finite(p_always) && is.finite(p_never) && (p_never  > 0) && (p_always == 0)

    one_sided <- only_always || only_never

    direction <- if (only_never) {
      "Never-takers only"
    } else if (only_always) {
      "Always-takers only"
    } else if (is.finite(p_always) && is.finite(p_never) && p_always > 0 && p_never > 0) {
      "Two-sided noncompliance"
    } else {
      "Perfect compliance"
    }

    # your encoding: eq==1 trivial for only never-takers; eq==0 trivial for only always-takers
    trivial_equation <- if (only_never) 1L else if (only_always) 0L else NA_integer_

    list(
      n_z0 = n_z0,
      n_z1 = n_z1,
      p_always = p_always,
      p_never = p_never,
      one_sided = one_sided,
      direction = direction,
      trivial_equation = trivial_equation
    )
  }

  Dcol <- D
  Zcol <- Z

  dt2 <- dt[!is.na(get(Dcol)) & !is.na(get(Zcol))]

  if (length(margins) == 0L) {
    res <- data.table::as.data.table(core(dt2))
  } else {
    res <- dt2[, core(.SD), by = margins]
  }

  os <- res[one_sided == TRUE]
  if (nrow(os) == 0L) {
    message(if (length(margins) == 0L)
      "No one-sided noncompliance detected overall (no margins)."
      else
        "No one-sided noncompliance detected in any margins cell."
    )
  } else {
    message(if (length(margins) == 0L)
      "One-sided noncompliance detected overall (no margins):"
      else
        "One-sided noncompliance detected in the following margins cells:"
    )
    print(os)
  }

  invisible(res)
}

#####binarize var #######

binarize_var <- function(data,
                         var,
                         ngroups,
                         gridtype = c("equisized", "equidistant"),
                         wvar = NA_character_,
                         newvar = NULL) {

  gridtype <- match.arg(gridtype)

  dt <- data.table::as.data.table(data)
  dt <- data.table::copy(dt)

  if (!var %in% names(dt)) {
    stop("`var` not found in `data`.")
  }

  if (!is.na(wvar) && !wvar %in% names(dt)) {
    stop("`wvar` not found in `data`.")
  }

  # decide output variable
  outvar <- if (is.null(newvar)) var else newvar

  vals_nonmiss <- dt[[var]][!is.na(dt[[var]])]

  if (length(unique(vals_nonmiss)) <= ngroups) {
    if (!identical(outvar, var)) {
      dt[, (outvar) := get(var)]
    }
    return(dt)
  }

  if (gridtype == "equidistant") {

    rng <- range(dt[[var]], na.rm = TRUE)

    breaks <- seq(
      from = rng[1] - 0.001,
      to   = rng[2] + 0.001,
      length.out = ngroups + 1
    )

    bins <- as.integer(cut(dt[[var]], breaks = breaks)) - 1L

  } else {

    w <- if (is.na(wvar)) NULL else dt[[wvar]]

    inner_breaks <- unique(
      weighted_quantile(
        x     = dt[[var]],
        w     = w,
        probs = (1:(ngroups - 1)) / ngroups
      )
    )

    breaks <- c(-Inf, inner_breaks, Inf)

    bins <- as.integer(cut(dt[[var]], breaks = breaks)) - 1L
  }

  dt[, (outvar) := bins]

  dt
}
#============ Liu and Xie (2019, JASA) cauchy combination test p-value ============#
#https://pubmed.ncbi.nlm.nih.gov/33012899/
cct_pvalue <- function(p, w = NULL,eps=1e-15) {
  p <- as.numeric(p)
  stopifnot(all(is.finite(p)), all(p >= 0), all(p <= 1))

  d <- length(p)
  if (is.null(w)) w <- rep(1 / d, d)
  w <- as.numeric(w)
  stopifnot(length(w) == d, all(w >= 0))
  w <- w / sum(w)

  # Stabilize extreme values (tan() blows up at p=0 or p=1)
  p2 <- pmin(pmax(p, eps), 1 - eps)

  Tstat <- sum(w * tan((0.5 - p2) * pi))
  p_cct <- 0.5 - atan(Tstat) / pi   # Eq. (3)
  pmin(pmax(p_cct, 0), 1)
}

