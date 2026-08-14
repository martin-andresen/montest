

##################################################### HELPER FUNCTIONS ######################################

## helper: all subsets of Y support
all_subsets <- function(vals, min_size = 1L, max_size = length(vals)) {
  out <- vector("list", 0L)
  for (k in min_size:max_size) {
    out <- c(out, combn(vals, k, simplify = FALSE))
  }
  out
}

append_to_main_rhs <- function(fml, y_name) {
  fml <- if (inherits(fml, "formula")) fml else stats::as.formula(fml)

  split_top_pipe <- function(expr) {
    if (is.call(expr) && identical(expr[[1L]], as.name("|"))) {
      c(split_top_pipe(expr[[2L]]), list(expr[[3L]]))
    } else {
      list(expr)
    }
  }

  rebuild_pipe <- function(parts) {
    if (length(parts) == 1L) return(parts[[1L]])
    Reduce(function(a, b) call("|", a, b), parts)
  }

  has_lhs <- length(fml) == 3L
  lhs <- if (has_lhs) fml[[2L]] else NULL
  rhs <- if (has_lhs) fml[[3L]] else fml[[2L]]

  parts <- split_top_pipe(rhs)

  parts[[1L]] <- if (identical(parts[[1L]], quote(1))) {
    as.name(y_name)
  } else {
    call("+", parts[[1L]], as.name(y_name))
  }

  new_rhs <- rebuild_pipe(parts)

  if (has_lhs) {
    stats::as.formula(call("~", lhs, new_rhs), env = environment(fml))
  } else {
    stats::as.formula(call("~", new_rhs), env = environment(fml))
  }
}

append_to_x_expr <- function(x_expr, y_name, env = parent.frame()) {
  fml <- if (is.null(x_expr) || identical(x_expr, quote(1))) {
    stats::as.formula(~ 1, env = env)
  } else {
    stats::as.formula(call("~", x_expr), env = env)
  }

  fml2 <- append_to_main_rhs(fml, y_name)
  fml2[[2L]]
}

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
    select  = NULL,
    cp = 0.0,
    maxrankcp = 10L,
    alpha = 0.05,
    prune = TRUE,
    minsize = 50L,
    screen = c("stepdown","none", "minimum", "negative", "nonpositive", "fgk_relevant"),
    store_trees = FALSE,
    verbose = FALSE,
    xval = 10,
    rpart_options = NULL,
    fe_expr = NULL,
    fe_rank_adj = TRUE
) {
  stopifnot(data.table::is.data.table(data))
  screen <- match.arg(screen)

  `%||%` <- function(x, y) if (is.null(x)) y else x

  sample_col <- as.character(sample)
  scores_col <- as.character(scores)
  stopifnot(length(sample_col) == 1L, sample_col %in% names(data))
  stopifnot(length(scores_col) == 1L, scores_col %in% names(data))

  x_names <- as.character(x_names)
  stopifnot(length(x_names) >= 1L, all(x_names %in% names(data)))

  if (is.null(margins)) margins <- character()
  margins <- unique(as.character(margins))
  if (length(margins) > 0L) stopifnot(all(margins %in% names(data)))

  weight_col <- if (is.null(weight)) NULL else as.character(weight)
  if (!is.null(weight_col)) stopifnot(length(weight_col) == 1L, weight_col %in% names(data))

  cluster_col <- if (is.null(cluster)) NULL else as.character(cluster)
  if (!is.null(cluster_col)) stopifnot(length(cluster_col) == 1L, cluster_col %in% names(data))

  stopifnot(is.numeric(minsize), length(minsize) == 1L, is.finite(minsize), minsize >= 1)
  minsize <- as.integer(minsize)

  stopifnot(is.numeric(alpha), length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  stopifnot(is.numeric(cp), length(cp) == 1L, is.finite(cp), cp >= 0)
  stopifnot(is.numeric(maxrankcp), length(maxrankcp) == 1L, is.finite(maxrankcp), maxrankcp >= 1)
  maxrankcp <- as.integer(maxrankcp)

  svec <- as.integer(data[[sample_col]])
  stopifnot(all(svec %in% c(1L, 2L)))

  # ---------------- validate select ----------------

  if (is.null(select)) select <- character()
  select <- unique(as.character(select))

  allowed_select <- c(sample_col, margins)

  bad_select <- setdiff(select, allowed_select)
  if (length(bad_select) > 0L) {
    stop(
      "`select` may only contain the sample column and variables in `margins`. Invalid entries: ",
      paste(bad_select, collapse = ", ")
    )
  }

  select_sample <- sample_col %in% select
  select_margins <- intersect(select, margins)

  # CART_test does not pool across sample halves.
  pool_sample <- FALSE
  pool_margins <- character()

  # Margins are always fit separately. There is no margin pooling in CART_test.
  sep_margins <- margins

  # Margins selected over are not final reporting strata.
  final_keep_margins <- setdiff(margins, select_margins)

  # ---------------- helpers ----------------

  select_leaves_local <- function(dt_train_leaf) {
    L <- nrow(dt_train_leaf)
    if (L == 0L) return(integer())

    if (screen == "none") {
      return(as.integer(dt_train_leaf$leaf))
    }

    if (screen == "minimum") {
      return(as.integer(dt_train_leaf[which.min(t), leaf]))
    }

    if (screen == "stepdown") {
      dtf <- dt_train_leaf[is.finite(t)]
      if (nrow(dtf) == 0L) return(integer())

      data.table::setorder(dtf, t)
      p <- stats::pnorm(dtf$t)

      m_keep <- 1L

      if (length(p) >= 2L) {
        for (m in 2L:length(p)) {
          # Holm on the current prefix of length m
          holm_ok <- all(p[seq_len(m)] <= alpha / rev(seq_len(m)))
          if (!holm_ok) break
          m_keep <- m
        }
      }

      return(as.integer(dtf$leaf[seq_len(m_keep)]))
    }

    if (screen == "negative") {
      thr <- stats::qnorm(alpha / L)
      keep <- dt_train_leaf[t <= thr, leaf]
      if (length(keep) == 0L) keep <- dt_train_leaf[which.min(t), leaf]
      return(as.integer(keep))
    }

    if (screen == "nonpositive") {
      thr <- stats::qnorm(1 - alpha / L)
      keep <- dt_train_leaf[t < thr, leaf]
      if (length(keep) == 0L) keep <- dt_train_leaf[which.min(t), leaf]
      return(as.integer(keep))
    }

    if (screen == "fgk_relevant") {
      dtf <- dt_train_leaf[is.finite(t)]
      if (nrow(dtf) == 0L) return(integer())

      L <- nrow(dtf)
      thr <- stats::qnorm(1 - alpha / L)

      keep <- dtf[t < thr, leaf]

      ## No local fallback. Fallback is handled later at the final
      ## non-selected margin/sample grouping level.
      return(as.integer(keep))
    }

    as.integer(dt_train_leaf$leaf)
  }

  choose_group_selection <- function(dt_grp) {
    dtf <- dt_grp[is.finite(t)]

    empty_out <- data.table::data.table(
      best_fit_id = NA_integer_,
      best_train_s = NA_integer_,
      sel_leaves = list(integer())
    )

    if (nrow(dtf) == 0L) return(empty_out)

    if (screen == "minimum") {
      j <- which.min(dtf$t)
      return(data.table::data.table(
        best_fit_id = as.integer(dtf$fit_id[j]),
        best_train_s = as.integer(dtf$train_s[j]),
        sel_leaves = list(as.integer(dtf$leaf[j]))
      ))
    }

    if (screen == "stepdown") {
      data.table::setorder(dtf, t)
      p <- stats::pnorm(dtf$t)

      m_keep <- 1L

      if (length(p) >= 2L) {
        for (m in 2L:length(p)) {
          holm_ok <- all(p[seq_len(m)] <= alpha / rev(seq_len(m)))
          if (!holm_ok) break
          m_keep <- m
        }
      }

      keep <- dtf[seq_len(m_keep)]

      return(keep[, .(
        sel_leaves = list(as.integer(leaf))
      ), by = .(
        best_fit_id = fit_id,
        best_train_s = train_s
      )])
    }

    if (screen == "none") {
      j <- which.min(dtf$t)
      bf <- as.integer(dtf$fit_id[j])
      bt <- as.integer(dtf$train_s[j])
      return(data.table::data.table(
        best_fit_id = bf,
        best_train_s = bt,
        sel_leaves = list(as.integer(dtf[fit_id == bf & train_s == bt, leaf]))
      ))
    }

    if (screen == "negative") {
      thr <- stats::qnorm(alpha / nrow(dtf))
      keep <- dtf[t <= thr]
      if (nrow(keep) == 0L) {
        j <- which.min(dtf$t)
        return(data.table::data.table(
          best_fit_id = as.integer(dtf$fit_id[j]),
          best_train_s = as.integer(dtf$train_s[j]),
          sel_leaves = list(as.integer(dtf$leaf[j]))
        ))
      }

      j <- which.min(keep$t)
      bf <- as.integer(keep$fit_id[j])
      bt <- as.integer(keep$train_s[j])
      return(data.table::data.table(
        best_fit_id = bf,
        best_train_s = bt,
        sel_leaves = list(as.integer(keep[fit_id == bf & train_s == bt, leaf]))
      ))
    }

    if (screen == "nonpositive") {
      thr <- stats::qnorm(1 - alpha / nrow(dtf))
      keep <- dtf[t < thr]
      if (nrow(keep) == 0L) {
        j <- which.min(dtf$t)
        return(data.table::data.table(
          best_fit_id = as.integer(dtf$fit_id[j]),
          best_train_s = as.integer(dtf$train_s[j]),
          sel_leaves = list(as.integer(dtf$leaf[j]))
        ))
      }

      j <- which.min(keep$t)
      bf <- as.integer(keep$fit_id[j])
      bt <- as.integer(keep$train_s[j])
      return(data.table::data.table(
        best_fit_id = bf,
        best_train_s = bt,
        sel_leaves = list(as.integer(keep[fit_id == bf & train_s == bt, leaf]))
      ))
    }

    if (screen == "fgk_relevant") {
      if ("fgk_keep" %in% names(dtf)) {
        keep <- dtf[fgk_keep == TRUE]

        if (nrow(keep) > 0L) {
          return(keep[, .(
            sel_leaves = list(as.integer(leaf))
          ), by = .(
            best_fit_id = fit_id,
            best_train_s = train_s
          )])
        }
      }

      ## Fallback only within the final choice group.
      ## Therefore if sample/yval/dval are in select, this can become one
      ## global fallback instead of one fallback per local tree.
      j <- which.min(dtf$t)
      data.table::data.table(
        best_fit_id = as.integer(dtf$fit_id[j]),
        best_train_s = as.integer(dtf$train_s[j]),
        sel_leaves = list(as.integer(dtf$leaf[j]))
      )
    }

    j <- which.min(dtf$t)
    data.table::data.table(
      best_fit_id = as.integer(dtf$fit_id[j]),
      best_train_s = as.integer(dtf$train_s[j]),
      sel_leaves = list(as.integer(dtf$leaf[j]))
    )
  }

  add_key_cols <- function(dt, key_dt, all_margin_names) {
    if (length(all_margin_names) == 0L) return(dt)

    if (is.null(key_dt) || ncol(key_dt) == 0L) {
      stop("add_key_cols(): key_dt is empty although margins are required: ",
           paste(all_margin_names, collapse = ", "))
    }

    missing_from_key <- setdiff(all_margin_names, names(key_dt))
    if (length(missing_from_key) > 0L) {
      stop("add_key_cols(): key_dt is missing margin columns: ",
           paste(missing_from_key, collapse = ", "))
    }

    for (cc in all_margin_names) {
      val <- key_dt[[cc]][1]
      if (is.na(val)) {
        print(key_dt)
        stop("add_key_cols(): margin column ", cc, " is NA in key_dt.")
      }
      dt[, (cc) := val]
    }

    dt
  }
  build_tree_design <- function(df) {
    as.data.frame(df[, ..x_names])
  }

  fit_one_tree_group <- function(df_group, train_sample) {
    dtr <- df_group[df_group[[sample_col]] == train_sample]
    if (nrow(dtr) < 2L) return(NULL)

    Xtr <- build_tree_design(dtr)
    df_tr <- data.frame(.y = as.numeric(dtr[[scores_col]]), Xtr, check.names = FALSE)
    w_tr <- if (!is.null(weight_col)) as.numeric(dtr[[weight_col]]) else NULL

    fml <- stats::as.formula(paste(".y ~", paste(colnames(Xtr), collapse = " + ")))

    tree <- rpart::rpart(
      formula = fml,
      data = df_tr,
      method = "anova",
      weights = w_tr,
      control = do.call(
        rpart::rpart.control,
        modifyList(
          rpart_options %||% list(),
          list(cp = cp, minbucket = minsize, xval = xval)
        )
      )
    )

    if (!is.null(tree$cptable) && nrow(tree$cptable) > 0L) {
      maxrankcp2 <- min(maxrankcp, nrow(tree$cptable))
      maxcp <- tree$cptable[maxrankcp2, 1]
      if (isTRUE(prune)) {
        opcpid <- which.min(tree$cptable[, 4])
        opcp <- tree$cptable[opcpid, 1]
        tree <- rpart::prune(tree, cp = opcp)
      } else {
        tree <- rpart::prune(tree, cp = maxcp)
      }
    }

    Xall <- build_tree_design(df_group)
    leaf_all <- as.integer(rpart:::pred.rpart(tree, as.matrix(Xall)))
    list(tree = tree, leaf_all = leaf_all)
  }

  global_means_one_cell <- function(df_cell, key_dt , idx) {
    dtg <- data.table::data.table(
      rowid  = idx,
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
      rank_adj <- if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
        fe_rank(
          data = data,
          idx = rowid,
          fe_expr = fe_expr,
          weight_col = weight_col,
          cluster_vals = if (!is.null(cluster_col)) data[[cluster_col]] else NULL
        )$rank
      } else {
        0L
      }

      o <- crv1_mean(score, w, cl, rank_adj = rank_adj)

      data.table::data.table(
        train = FALSE,
        G = o$G,
        N = o$N,
        fe_rank = rank_adj,
        df = o$df,
        coef = o$coef,
        stderr = o$se,
        t = o$t,
        p.raw = stats::pnorm(o$t)
      )
    }, by = "sample"]

    if (!is.null(key_dt) && ncol(key_dt) > 0L) {
      for (cc in names(key_dt)) out[, (cc) := key_dt[[cc]][1]]
      data.table::setcolorder(out, c(names(key_dt), "train", "sample"))
    }
    out
  }

  # ---------------- global output on fully separated cells ----------------
  n <- nrow(data)
  gid_full <- if (length(margins) == 0L) {
    rep(1L, n)
  } else {
    data.table::frankv(data[, ..margins], ties.method = "dense")
  }
  idx_full_list <- split(seq_len(n), gid_full)
  n_full_cells <- length(idx_full_list)

  global_list <- vector("list", n_full_cells)
  for (g in seq_len(n_full_cells)) {
    idx <- idx_full_list[[g]]
    key_dt <- if (length(margins) > 0L) data[idx[1L], ..margins] else NULL
    df_cell <- data[idx]
    global_list[[g]] <- global_means_one_cell(df_cell, key_dt, idx)
  }
  global_out <- data.table::rbindlist(global_list, use.names = TRUE, fill = TRUE)

  # ---------------- fit-groups for trees ----------------
  gid_fit <- if (length(sep_margins) == 0L) {
    rep(1L, n)
  } else {
    data.table::frankv(data[, ..sep_margins], ties.method = "dense")
  }
  idx_fit_list <- split(seq_len(n), gid_fit)
  n_fit <- length(idx_fit_list)

  fit_objs   <- vector("list", n_fit)
  trees_list <- if (store_trees) vector("list", n_fit) else NULL

  for (g in seq_len(n_fit)) {
    idx <- idx_fit_list[[g]]
    df_group <- data[idx]
    key_dt <- if (length(sep_margins) > 0L) data[idx[1L], ..sep_margins] else NULL

    dirs <- vector("list", 2L)
    if (store_trees) trees_list[[g]] <- vector("list", 2L)

    for (train_s in c(1L, 2L)) {
      est_s <- if (train_s == 1L) 2L else 1L

      fit <- fit_one_tree_group(df_group, train_s)
      if (is.null(fit)) {
        dirs[[train_s]] <- NULL
        next
      }
      if (store_trees) trees_list[[g]][[train_s]] <- fit$tree

      leaf_all <- fit$leaf_all
      s_group <- as.integer(df_group[[sample_col]])
      score_group <- as.numeric(df_group[[scores_col]])
      w_group <- if (!is.null(weight_col)) as.numeric(df_group[[weight_col]]) else rep(1.0, nrow(df_group))
      w_group[!is.finite(w_group)] <- 0
      cl_group <- if (!is.null(cluster_col)) df_group[[cluster_col]] else NULL

      dt0 <- data.table::data.table(
        rowid  = idx,
        leaf   = leaf_all,
        sample = s_group,
        score  = score_group,
        w      = w_group
      )
      if (!is.null(cluster_col)) {
        dt0[, cl := as.integer(factor(cl_group, exclude = NULL))]
      } else {
        dt0[, cl := .I]
      }

      dt_train <- dt0[sample == train_s & !is.na(leaf)]

      train_leaf <- dt_train[, {
        G_here <- data.table::uniqueN(cl)
        N_here <- .N

        if (G_here < minsize) {
          list(
            G = G_here,
            N = N_here,
            fe_rank = NA_integer_,
            df = NA_real_,
            coef = NA_real_,
            stderr = NA_real_,
            t = NA_real_
          )
        } else {
          rank_adj <- 0L

          if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
            rank_adj <- fe_rank(
              data = data,
              idx = rowid,
              fe_expr = fe_expr,
              weight_col = weight_col,
              cluster_vals = if (!is.null(cluster_col)) data[[cluster_col]] else NULL
            )$rank
          }

          o <- crv1_mean(
            score = score,
            w = w,
            cl = cl,
            rank_adj = rank_adj
          )

          list(
            G = o$G,
            N = o$N,
            fe_rank = rank_adj,
            df = o$df,
            coef = o$coef,
            stderr = o$se,
            t = o$t
          )
        }
      }, by = leaf]

      train_leaf[, `:=`(train = TRUE, sample = train_s)]
      train_leaf[, p.raw := stats::pnorm(t)]

      prom <- select_leaves_local(train_leaf[is.finite(t), .(leaf, t)])
      train_leaf <- add_key_cols(train_leaf, key_dt, margins)

      dirs[[train_s]] <- list(
        idx = idx,
        key_dt = key_dt,
        leaf_all = leaf_all,
        train_leaf = train_leaf,
        prom = as.integer(prom),
        fit_id = g,
        train_s = train_s,
        est_s = est_s
      )
    }

    fit_objs[[g]] <- dirs

    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("CART_test: processed %d/%d fit groups", g, n_fit))
    }
  }

  # ---------------- collect locally screened candidate leaves ----------------
  leaf_rows <- list()
  all_leaf_rows <- list()
  lr <- 0L
  alr <- 0L

  for (g in seq_len(n_fit)) {
    for (train_s in c(1L, 2L)) {
      obj <- fit_objs[[g]][[train_s]]
      if (is.null(obj)) next

      tl_all <- data.table::copy(obj$train_leaf[is.finite(t)])
      if (nrow(tl_all) == 0L) next

      tl_all[, `:=`(
        fit_id = g,
        train_s = train_s,
        est_s = obj$est_s,
        fgk_keep = leaf %in% obj$prom
      )]

      all_leaf_rows[[alr <- alr + 1L]] <- tl_all

      if (length(obj$prom) == 0L) next

      tl <- data.table::copy(
        obj$train_leaf[leaf %in% obj$prom & is.finite(t)]
      )
      if (nrow(tl) == 0L) next

      tl[, `:=`(
        fit_id = g,
        train_s = train_s,
        est_s = obj$est_s,
        fgk_keep = TRUE
      )]

      leaf_rows[[lr <- lr + 1L]] <- tl
    }
  }

  if (length(all_leaf_rows) == 0L) {
    stop("Testing subset is empty: no finite training leaf statistics.")
  }

  leaf_tbl_all <- data.table::rbindlist(
    all_leaf_rows,
    use.names = TRUE,
    fill = TRUE
  )

  leaf_tbl_use <- if (screen == "fgk_relevant") {
    leaf_tbl_all
  } else {
    if (length(leaf_rows) == 0L) {
      stop("Testing subset is empty (after local screening).")
    }
    data.table::rbindlist(leaf_rows, use.names = TRUE, fill = TRUE)
  }

  if (nrow(leaf_tbl_use) == 0L) {
    stop("Testing subset is empty (after local screening).")
  }
  # ---------------- choose among selected margins ----------------
  choice_group_cols <- c(final_keep_margins)
  if (!select_sample) choice_group_cols <- c(choice_group_cols, "train_s")
  choice_group_cols <- unique(choice_group_cols)

  if (length(choice_group_cols) == 0L) {

    sel_tbl <- choose_group_selection(
      leaf_tbl_use[, .(fit_id, train_s, leaf, t, fgk_keep)]
    )[, .(
      best_fit_id,
      best_train_s,
      sel_leaves
    )]

  } else {

    sel_tbl <- leaf_tbl_use[
      ,
      choose_group_selection(.SD[, .(fit_id, train_s, leaf, t, fgk_keep)]),
      by = choice_group_cols
    ]

    key_row <- unique(leaf_tbl_use[, ..choice_group_cols])

    sel_tbl <- merge.data.table(
      key_row,
      sel_tbl,
      by = choice_group_cols,
      all.x = TRUE
    )
  }

  # ---------------- mark relevant train rows and build testing jobs ----------------
  train_rows <- list()
  tr_k <- 0L
  test_jobs <- list()
  tj_k <- 0L

  for (g in seq_len(n_fit)) {
    for (train_s in c(1L, 2L)) {
      obj <- fit_objs[[g]][[train_s]]
      if (is.null(obj)) next

      key_row <- data.table::data.table()
      if (length(final_keep_margins) > 0L && !is.null(obj$key_dt) && ncol(obj$key_dt) > 0L) {
        for (cc in final_keep_margins) {
          if (cc %in% names(obj$key_dt)) key_row[, (cc) := obj$key_dt[[cc]][1]]
        }
      }
      if (!select_sample) key_row[, train_s := train_s]

      if (length(choice_group_cols) == 0L) {
        sel_here <- sel_tbl
      } else {
        sel_here <- merge(
          key_row,
          sel_tbl,
          by = choice_group_cols,
          all.x = TRUE
        )
      }

      sel_leaves <- integer()
      chosen <- FALSE

      if (nrow(sel_here) > 0L) {
        sel_match <- sel_here[
          is.finite(best_fit_id) &
            is.finite(best_train_s) &
            best_fit_id == g &
            best_train_s == train_s
        ]

        if (nrow(sel_match) > 0L) {
          sel_leaves <- as.integer(unlist(sel_match$sel_leaves, use.names = FALSE))
          sel_leaves <- unique(sel_leaves)
          chosen <- length(sel_leaves) > 0L
        }
      }

      tl <- data.table::copy(obj$train_leaf)
      tl[, relevant := as.integer(leaf %in% sel_leaves)]
      train_rows[[tr_k <- tr_k + 1L]] <- tl

      if (chosen) {
        for (lf in sel_leaves) {
          job_dt <- data.table::data.table(
            fit_id   = g,
            train_s  = train_s,
            est_s    = obj$est_s,
            leaf     = as.integer(lf)
          )

          # margins kept in reported output
          for (cc in margins) {
            val <- NA
            if (!is.null(obj$key_dt) && cc %in% names(obj$key_dt)) val <- obj$key_dt[[cc]][1]
            job_dt[, (cc) := val]
          }

          # pooling key for final testing
          for (cc in final_keep_margins) {
            val <- NA
            if (!is.null(obj$key_dt) && cc %in% names(obj$key_dt)) {
              val <- obj$key_dt[[cc]][1]
            }
            job_dt[, (paste0(".pool_", cc)) := val]
          }

          # Leaf IDs are only meaningful within a fitted tree and training direction.
          # Never pool numeric leaf IDs across different trees.
          job_dt[, .pool_leaf := paste(g, train_s, as.integer(lf), sep = "_")]

          # CART_test never pools across sample halves.
          job_dt[, .pool_est_s := obj$est_s]

          test_jobs[[tj_k <- tj_k + 1L]] <- job_dt
        }
      }
    }
  }

  train_out <- data.table::rbindlist(train_rows, use.names = TRUE, fill = TRUE)

  jobs_dt <- data.table::rbindlist(test_jobs, use.names = TRUE, fill = TRUE)

  # ---------------- collapse testing jobs by selected tree/leaf/sample direction ----------------
  pool_cols <- grep("^\\.pool_", names(jobs_dt), value = TRUE)

  # ---------------- testing: one test per pooled job-group ----------------
  test_rows <- list()
  tt_k <- 0L

  if (length(pool_cols) == 0L) {
    unique_job_groups <- data.table::data.table(.dummy = 1L)
    jobs_dt[, .dummy := 1L]
    pool_cols <- ".dummy"
  } else {
    unique_job_groups <- unique(jobs_dt[, ..pool_cols])
  }

  for (jj in seq_len(nrow(unique_job_groups))) {
    key <- unique_job_groups[jj]

    sel_expr <- rep(TRUE, nrow(jobs_dt))
    for (cc in pool_cols) {
      v <- key[[cc]][1]
      if (is.na(v)) {
        sel_expr <- sel_expr & is.na(jobs_dt[[cc]])
      } else {
        sel_expr <- sel_expr & (jobs_dt[[cc]] == v)
      }
    }
    jobs_here <- jobs_dt[sel_expr]
    if (nrow(jobs_here) == 0L) next

    idx_keep_all <- integer()

    for (kk in seq_len(nrow(jobs_here))) {
      job <- jobs_here[kk]
      obj <- fit_objs[[job$fit_id]][[job$train_s]]

      idx <- obj$idx
      leaf_all <- obj$leaf_all

      keep_local <- which(
        svec[idx] == job$est_s &
          !is.na(leaf_all) &
          leaf_all == job$leaf
      )
      if (length(keep_local) == 0L) next

      idx_keep_all <- c(idx_keep_all, idx[keep_local])
    }

    idx_keep_all <- sort(unique(idx_keep_all))
    if (length(idx_keep_all) == 0L) next

    y <- as.numeric(data[[scores_col]])[idx_keep_all]
    w <- if (!is.null(weight_col)) as.numeric(data[[weight_col]])[idx_keep_all] else rep(1.0, length(idx_keep_all))
    w[!is.finite(w)] <- 0
    cl <- if (!is.null(cluster_col)) data[[cluster_col]][idx_keep_all] else NULL

    rank_adj <- if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
      fe_rank(
        data = data,
        idx = idx_keep_all,
        fe_expr = fe_expr,
        weight_col = weight_col,
        cluster_vals = if (!is.null(cluster_col)) data[[cluster_col]] else NULL
      )$rank
    } else {
      0L
    }

    o <- crv1_mean(y, w, cl, rank_adj = rank_adj)

    sample_here <- unique(jobs_here$est_s)
    leaf_here   <- unique(jobs_here$leaf)

    rr <- data.table::data.table(
      train = FALSE,
      sample = if (length(sample_here) == 1L) as.integer(sample_here) else NA_integer_,
      leaf = if (length(leaf_here) == 1L) as.integer(leaf_here) else NA_integer_,
      G = o$G,
      N = o$N,
      fe_rank = rank_adj,
      df = o$df,
      coef = o$coef,
      stderr = o$se,
      t = o$t,
      p.raw = stats::pnorm(o$t),
      relevant = NA_integer_
    )

    # carry reported margins forward
    for (cc in margins) {
      vals <- unique(jobs_here[[cc]])
      vals <- vals[!is.na(vals)]
      rr[, (cc) := if (length(vals) == 1L) vals else NA]
    }

    test_rows[[tt_k <- tt_k + 1L]] <- rr
  }

  if (".dummy" %in% names(jobs_dt)) jobs_dt[, .dummy := NULL]

  if (length(test_rows) == 0L) stop("Testing subset is empty (after selection).")
  test_out <- data.table::rbindlist(test_rows, use.names = TRUE, fill = TRUE)

  results_out <- data.table::rbindlist(list(train_out, test_out), use.names = TRUE, fill = TRUE)
  if (!("relevant" %in% names(results_out))) results_out[, relevant := NA_integer_]
  results_out[train == FALSE, relevant := NA_integer_]
  data.table::setorder(results_out, train)
  if ("selected_local" %in% names(results_out)) {
    results_out[, selected_local := NULL]
  }

  front_cols <- c(
    intersect(margins, names(results_out)),
    intersect(sample_col, names(results_out)),
    intersect("leaf", names(results_out)),
    intersect("train", names(results_out))
  )

  other_cols <- setdiff(names(results_out), front_cols)
  data.table::setcolorder(results_out, c(front_cols, other_cols))

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
                             strat_col = NULL,     # NULL = no stratification
                             diag_prefix = "cf_",
                             verbose = TRUE) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.numeric(K), length(K) == 1L, K >= 2)
  stopifnot(is.character(fold_col), length(fold_col) == 1L)

  if (!is.null(by_col)) {
    stopifnot(is.character(by_col), length(by_col) == 1L, by_col %in% names(DT))
  }
  if (!is.null(strat_col)) {
    stopifnot(is.character(strat_col), length(strat_col) == 1L, strat_col %in% names(DT))
  }
  if (!is.null(cluster_name)) {
    stopifnot(is.character(cluster_name), length(cluster_name) == 1L, cluster_name %in% names(DT))
  }

  if (!is.null(cluster_name) && !is.null(strat_col)) {
    stop("make_group_folds(): stratification and clustering cannot currently be used together.")
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

  split_by <- c(by_col, strat_col)

  if (length(split_by) == 0L) {
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
    # ---------- WITHIN-BLOCK SPLIT ----------
    if (has_diag) {
      DT[, c(fold_col, k_col, g_col, n_col) := {
        cl <- if (!is.null(cluster_name)) get(cluster_name) else NULL
        make_one(.N, cl)
      }, by = split_by]
    } else {
      DT[, (fold_col) := {
        cl <- if (!is.null(cluster_name)) get(cluster_name) else NULL
        make_one(.N, cl)
      }, by = split_by]
    }

    if (verbose) {
      if (!is.null(by_col) && !is.null(strat_col)) {
        message(sprintf(
          "Cross-fitting folds: created '%s' with up to %d folds within '%s', stratified by '%s'.",
          fold_col, K, by_col, strat_col
        ))
      } else if (!is.null(by_col)) {
        message(sprintf(
          "Cross-fitting folds: created '%s' with up to %d folds within '%s'.",
          fold_col, K, by_col
        ))
      } else {
        message(sprintf(
          "Cross-fitting folds: created '%s' with up to %d folds stratified by '%s'.",
          fold_col, K, strat_col
        ))
      }
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
                         sample_name = "sample",
                         forest_opts = list(),
                         hat_suffix = ".hat",
                         mode = c("within", "across"),
                         verbose = FALSE) {

  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.character(y_name), length(y_name) == 1L)
  stopifnot(is.character(x_names), length(x_names) >= 1L)
  stopifnot(is.character(sample_name), length(sample_name) == 1L)
  stopifnot(sample_name %in% names(DT))
  stopifnot(all(DT[[sample_name]] %in% c(1L, 2L)))

  mode <- match.arg(mode)

  if (!is.null(folds)) {
    stopifnot(is.character(folds), length(folds) == 1L, folds %in% names(DT))
  }

  if (is.null(i)) i <- DT[, .I]
  stopifnot(is.integer(i) || is.numeric(i))
  i <- as.integer(i)
  n_i <- length(i)
  if (n_i == 0L) return(invisible(DT))

  if (is.null(margins)) margins <- character()
  grp_cols <- as.character(margins)

  missing_cols <- setdiff(
    c(sample_name, y_name, x_names, grp_cols,
      if (!is.null(weight_name)) weight_name,
      if (!is.null(folds)) folds),
    names(DT)
  )
  if (length(missing_cols)) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  hat_col <- paste0(y_name, hat_suffix)
  if (!(hat_col %in% names(DT))) DT[, (hat_col) := NA_real_]

  X_all <- as.matrix(DT[i, ..x_names])

  rid_col <- ".__rid__"
  DT[i, (rid_col) := seq_len(n_i)]
  rid <- DT[[rid_col]]

  sample_all <- DT[[sample_name]]
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
      c(
        list(
          X = X_tr,
          Y = y_tr,
          sample.weights = w_tr,
          compute.oob.predictions = FALSE
        ),
        forest_opts
      )
    )

    X_te <- X_all[rid[idx_te], , drop = FALSE]
    preds[] <- as.numeric(predict(rf, X_te)$predictions)
    preds
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
      c(
        list(
          X = X_all[rid[idx_s], , drop = FALSE],
          Y = y_s,
          sample.weights = w_s,
          compute.oob.predictions = TRUE
        ),
        forest_opts
      )
    )

    p[] <- as.numeric(predict(fit)$predictions)
    p
  }

  within_pred_idx <- function(idx_s) {
    n <- length(idx_s)
    if (n == 0L) return(numeric())

    f <- f_all[idx_s]

    if (data.table::uniqueN(f) < 2L) {
      return(within_oob_idx(idx_s))
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

        mu <- if (length(y_tr) == 0L) {
          NA_real_
        } else {
          if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
        }

        p[te_pos] <- mu
      } else {
        p[te_pos] <- fit_predict_idx(idx_tr, idx_te)
      }
    }

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

leave_cluster_out_mean <- function(DT,
                                   y_name,
                                   cluster_name = NULL,
                                   by = NULL,
                                   out_name = paste0(y_name, ".hat"),
                                   weight_name = NULL) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.character(y_name), length(y_name) == 1L, y_name %in% names(DT))
  stopifnot(is.character(out_name), length(out_name) == 1L)

  by <- unique(as.character(by %||% character()))
  by <- by[by %in% names(DT)]

  if (!is.null(cluster_name)) {
    stopifnot(is.character(cluster_name), length(cluster_name) == 1L)
    stopifnot(cluster_name %in% names(DT))
  }

  if (!is.null(weight_name)) {
    stopifnot(is.character(weight_name), length(weight_name) == 1L)
    stopifnot(weight_name %in% names(DT))
  }

  if (!(out_name %in% names(DT))) {
    DT[, (out_name) := NA_real_]
  }

  cols <- unique(c(y_name, cluster_name, weight_name))

  DT[
    ,
    (out_name) := {
      y <- .SD[[y_name]]
      w <- if (!is.null(weight_name)) .SD[[weight_name]] else rep(1, .N)

      ok <- is.finite(y) & is.finite(w) & w > 0

      if (!any(ok)) {
        rep(NA_real_, .N)

      } else if (is.null(cluster_name)) {
        ## Leave-one-observation-out weighted mean.
        wy <- ifelse(ok, w * y, 0)
        ww <- ifelse(ok, w, 0)

        total_wy <- sum(wy)
        total_w  <- sum(ww)

        denom <- total_w - ww
        out <- (total_wy - wy) / denom
        out[!is.finite(out) | denom <= 0] <- total_wy / total_w
        out

      } else {
        ## Leave-one-cluster-out weighted mean.
        cl <- .SD[[cluster_name]]

        wy <- ifelse(ok, w * y, 0)
        ww <- ifelse(ok, w, 0)

        total_wy <- sum(wy)
        total_w  <- sum(ww)

        cluster_wy <- rowsum(wy, cl, reorder = FALSE)
        cluster_w  <- rowsum(ww, cl, reorder = FALSE)

        cl_key <- as.character(cl)
        denom <- total_w - as.numeric(cluster_w[cl_key, 1L])
        numer <- total_wy - as.numeric(cluster_wy[cl_key, 1L])

        out <- numer / denom
        out[!is.finite(out) | denom <= 0] <- total_wy / total_w
        out
      }
    },
    by = by,
    .SDcols = cols
  ]

  invisible(DT)
}


##VALIDATE IV FORMULA
split_top_pipe <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("|"))) {
    c(split_top_pipe(expr[[2L]]), list(expr[[3L]]))
  } else {
    list(expr)
  }
}

rebuild_pipe <- function(parts) {
  if (length(parts) == 0L) return(NULL)
  if (length(parts) == 1L) return(parts[[1L]])
  Reduce(function(a, b) call("|", a, b), parts)
}

parse_iv_formula <- function(fml) {
  fml <- if (inherits(fml, "formula")) {
    fml
  } else {
    stats::as.formula(fml, env = parent.frame())
  }

  env <- environment(fml)

  ## Handles R's nested parse of:
  ##   Y ~ X | D ~ Z
  ##   Y ~ X | fe | D ~ Z
  ## and one-sided versions:
  ##   ~ X | D ~ Z
  ##   ~ X | fe | D ~ Z
  ##
  ## R parses these roughly as:
  ##   (Y ~ X | D) ~ Z
  ##   (Y ~ X | fe | D) ~ Z
  if (length(fml) == 3L &&
      is.call(fml[[2L]]) &&
      identical(fml[[2L]][[1L]], as.name("~"))) {

    left_fml <- fml[[2L]]

    if (length(left_fml) == 3L) {
      has_lhs <- TRUE
      lhs <- left_fml[[2L]]
      left_rhs <- left_fml[[3L]]
    } else if (length(left_fml) == 2L) {
      has_lhs <- FALSE
      lhs <- NULL
      left_rhs <- left_fml[[2L]]
    } else {
      stop("Malformed nested formula.", call. = FALSE)
    }

    z_expr <- fml[[3L]]
    left_parts <- split_top_pipe(left_rhs)

    ## Last left part is endogenous variable D.
    ## Everything before it is the reduced-form RHS,
    ## including fixed effects if present.
    d_expr <- left_parts[[length(left_parts)]]
    rf_parts <- left_parts[-length(left_parts)]

    iv <- call("~", d_expr, z_expr)

    return(list(
      formula = fml,
      env = env,
      has_lhs = has_lhs,
      lhs = lhs,
      rf_parts = rf_parts,
      iv = iv
    ))
  }

  ## Normal parse, e.g. with parenthesized IV part:
  ##   Y ~ X | fe | (D ~ Z)
  has_lhs <- length(fml) == 3L

  lhs <- if (has_lhs) fml[[2L]] else NULL
  rhs <- if (has_lhs) fml[[3L]] else fml[[2L]]

  parts <- split_top_pipe(rhs)

  iv_flags <- vapply(
    parts,
    function(x) is.call(x) && identical(x[[1L]], as.name("~")),
    logical(1L)
  )

  iv <- if (any(iv_flags)) {
    parts[[tail(which(iv_flags), 1L)]]
  } else {
    NULL
  }

  rf_parts <- if (any(iv_flags)) {
    parts[!iv_flags]
  } else {
    parts
  }

  list(
    formula = fml,
    env = env,
    has_lhs = has_lhs,
    lhs = lhs,
    rf_parts = rf_parts,
    iv = iv
  )
}

validate_iv <- function(fml, data) {
  parsed <- parse_iv_formula(fml)

  errors <- character()

  lhs <- parsed$lhs
  rf_parts <- parsed$rf_parts
  iv <- parsed$iv
  env <- parsed$env
  has_lhs <- parsed$has_lhs

  Y <- if (has_lhs) all.vars(lhs) else NULL

  if (has_lhs && length(Y) != 1L) {
    errors <- c(errors, "If supplied, the outcome side must contain exactly one variable.")
  }

  if (is.null(iv)) {
    errors <- c(errors, "No IV part found. Expected an endogenous variable with instrument, e.g. D ~ Z.")
    D <- NULL
    Z <- NULL
  } else {
    D <- all.vars(iv[[2L]])
    Z <- all.vars(iv[[3L]])

    if (length(D) != 1L) {
      errors <- c(errors, "Exactly one endogenous variable is required.")
    }

    if (!is.symbol(iv[[2L]])) {
      errors <- c(errors, "The endogenous variable must be a bare non-empty variable name, e.g. D ~ Z.")
    }

    if (length(Z) < 1L) {
      errors <- c(errors, "At least one instrument is required.")
    }
  }

  if (length(rf_parts) == 0L) {
    errors <- c(errors, "No reduced-form RHS found before the IV part.")
    covars_expr <- NULL
    fe_expr <- NULL
  } else {
    covars_expr <- rf_parts[[1L]]
    fe_expr <- if (length(rf_parts) >= 2L) rf_parts[[2L]] else NULL
  }

  X <- if (!is.null(covars_expr)) {
    unique(all.vars(covars_expr))
  } else {
    character()
  }

  FE <- if (!is.null(fe_expr)) {
    unique(all.vars(fe_expr))
  } else {
    character()
  }

  variables <- unique(c(Y, X, FE, D, Z))
  missing <- setdiff(variables, names(data))

  if (length(missing)) {
    errors <- c(errors, paste("Missing variables:", paste(missing, collapse = ", ")))
  }

  if (length(errors)) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }

  rf_rhs_expr <- rebuild_pipe(rf_parts)

  reduced_form_rhs <- stats::as.formula(
    call("~", rf_rhs_expr),
    env = env
  )

  reduced_form_formula <- if (has_lhs) {
    stats::as.formula(
      call("~", lhs, rf_rhs_expr),
      env = env
    )
  } else {
    NULL
  }

  # New: formula containing only the main X part before the first pipe.
  # This is useful as the default for fml.Z/fml.Q and for forest features.
  X_formula <- if (!is.null(covars_expr)) {
    stats::as.formula(
      call("~", covars_expr),
      env = env
    )
  } else {
    NULL
  }

  # New: formula containing only the FE part, if present.
  FE_formula <- if (!is.null(fe_expr)) {
    stats::as.formula(
      call("~", fe_expr),
      env = env
    )
  } else {
    NULL
  }

  null_if_empty <- function(x) {
    if (is.null(x) || length(x) == 0L) NULL else x
  }

  list(
    ok = TRUE,
    errors = NULL,
    formula = parsed$formula,

    reduced_form_rhs = reduced_form_rhs,
    reduced_form_formula = reduced_form_formula,

    # New / important for the refactor
    X_expr = covars_expr,
    X_formula = X_formula,
    FE_expr = fe_expr,
    FE_formula = FE_formula,

    variables = null_if_empty(variables),
    X = null_if_empty(X),
    Y = null_if_empty(Y),
    D = null_if_empty(D),
    Z = null_if_empty(Z),
    FE = null_if_empty(FE),

    has_lhs = has_lhs,
    iv = iv
  )
}

semiparametric_hat_stage <- function(DT,
                                     y_tilde,
                                     x_names,
                                     out_sp_hat,
                                     has_fe = FALSE,
                                     by_nuis = NULL,
                                     weight = NULL,
                                     cluster = NULL,
                                     folds = NULL,
                                     crossfit_mode = "within",
                                     forest_opts = list()) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(y_tilde %in% names(DT))

  by_nuis <- unique(as.character(by_nuis %||% character()))
  by_nuis <- by_nuis[by_nuis %in% names(DT)]

  has_x <- !is.null(x_names) && length(x_names) > 0L

  if (has_x) {
    crossfit_hat(
      DT,
      i = i,
      y_name = y_for_rf,
      x_names = x_names,
      folds = foldname,
      sample_name = sample_var,
      margins = by_nuis,
      weight_name = weight,
      mode = ifelse(
        (!is.null(crossfit_label) &&
           crossfit_label %in% crossfit &&
           is.null(foldname)),
        "across",
        "within"
      ),
      forest_opts = forest_opts,
      hat_suffix = ".sp_hat"
    )

    DT[, (out_sp_hat) := get(paste0(y_tilde, ".sp_hat"))]
  } else if (has_fe) {
    DT[, (out_sp_hat) := 0]
  } else {
    oos_mean_hat_by(
      DT = DT,
      y_name = y_tilde,
      out_name = out_sp_hat,
      by = by_nuis,
      weight = weight,
      cluster = cluster,
      folds = folds
    )
  }

  invisible(DT)
}


feols_partial_out <- function(DT,
                              y,
                              rhs_expr = quote(1),
                              fe_expr = NULL,
                              by = NULL,
                              weight = NULL,
                              prefix = NULL,
                              keep = c("resid", "fitted", "both"),
                              fixest_opts = list()) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(all(y %in% names(DT)))

  keep <- match.arg(keep)

  by <- unique(as.character(by %||% character()))
  by <- by[by %in% names(DT)]

  ## Singleton/perfect-fit FE removal only affects reported SEs, not the
  ## fitted/residual values we actually use here ??? but it does silently
  ## shrink the row count `residuals()`/`fitted()` come back with. We only
  ## ever consume residuals/fitted from this fit, so default to keeping
  ## every row unless the caller explicitly asked for something else.
  fixest_opts <- utils::modifyList(list(fixef.rm = "none"), fixest_opts)

  rhs_txt <- if (is.null(rhs_expr) || identical(rhs_expr, quote(1))) {
    "1"
  } else {
    deparse1(rhs_expr)
  }

  has_fe <- !is.null(fe_expr)

  out <- data.table::data.table(
    y = y,
    resid = paste0(prefix %||% y, "_resid"),
    fitted = paste0(prefix %||% y, "_fitted")
  )

  if (keep %in% c("resid", "both")) {
    for (v in out$resid) DT[, (v) := NA_real_]
  }
  if (keep %in% c("fitted", "both")) {
    for (v in out$fitted) DT[, (v) := NA_real_]
  }

  # No FE and no RHS means residual = original, fitted = 0.
  # This is useful for the FWL path where there is nothing to partial out.
  if (!has_fe && rhs_txt == "1") {
    if (keep %in% c("resid", "both")) {
      for (j in seq_along(y)) {
        DT[, (out$resid[j]) := get(y[j])]
      }
    }
    if (keep %in% c("fitted", "both")) {
      for (j in seq_along(y)) {
        DT[, (out$fitted[j]) := 0]
      }
    }

    return(list(
      resid = if (keep %in% c("resid", "both")) out$resid else NULL,
      fitted = if (keep %in% c("fitted", "both")) out$fitted else NULL,
      map = out
    ))
  }

  fml_for <- function(lhs) {
    lhs_q <- paste0("`", gsub("`", "``", lhs), "`")   # backtick-quote so leading "_" etc. can't break the parser
    if (has_fe) {
      stats::as.formula(paste0(lhs_q, " ~ ", rhs_txt, " | ", deparse1(fe_expr)))
    } else {
      stats::as.formula(paste0(lhs_q, " ~ ", rhs_txt))
    }
  }

  idx_list <- if (!length(by)) {
    list(seq_len(nrow(DT)))
  } else {
    DT[, .(idx = list(.I)), by = by]$idx
  }

  for (ii in idx_list) {
    dsub <- DT[ii]

    for (j in seq_along(y)) {
      lhs <- y[j]
      fml <- fml_for(lhs)

      ## feols() drops any row with NA in a variable entering the formula
      ## (or a NA/non-positive weight) and returns residuals()/fitted()
      ## only for the rows it kept. Filter to those rows up front so the
      ## length lines up with what we assign back into DT.
      vars_needed <- intersect(all.vars(fml), names(dsub))
      keep_row <- stats::complete.cases(dsub[, ..vars_needed])
      if (!is.null(weight)) {
        keep_row <- keep_row & !is.na(dsub[[weight]]) & dsub[[weight]] > 0
      }

      dsub_j <- dsub[keep_row]
      ii_j   <- ii[keep_row]

      args <- c(
        list(fml = fml, data = dsub_j, notes = FALSE, warn = FALSE),
        fixest_opts
      )

      if (!is.null(weight)) {
        args$weights <- stats::as.formula(paste0("~", weight))
      }

      fit <- do.call(fixest::feols, args)

      if (keep %in% c("resid", "both")) {
        DT[ii_j, (out$resid[j]) := as.numeric(stats::residuals(fit))]
      }

      if (keep %in% c("fitted", "both")) {
        DT[ii_j, (out$fitted[j]) := as.numeric(stats::fitted(fit))]
      }
    }
  }

  list(
    resid = if (keep %in% c("resid", "both")) out$resid else NULL,
    fitted = if (keep %in% c("fitted", "both")) out$fitted else NULL,
    map = out
  )
}


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


parse_one_sided_rhs <- function(fml, name) {
  if (is.null(fml)) return(NULL)

  if (!inherits(fml, "formula") || length(fml) != 2L) {
    stop(name, " must be a one-sided formula, e.g. ~ X1 + X2.", call. = FALSE)
  }

  fml[[2L]]
}

null_if_empty <- function(x) {
  if (is.null(x) || length(x) == 0L) NULL else x
}

make_safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


estimate_conditional_mean <- function(DT,
                                      y_name,
                                      x_expr,
                                      fe_expr = NULL,
                                      out_hat = paste0(y_name, ".hat"),
                                      out_resid = NULL,
                                      by = NULL,
                                      sample_var = "sample",
                                      weight = NULL,
                                      cluster = NULL,
                                      parametric = FALSE,
                                      foldname = NULL,
                                      crossfit = NULL,
                                      crossfit_label = NULL,
                                      forest_opts = list(),
                                      fixest_opts = list(),
                                      x_names = NULL,
                                      x_prefix = "__x",
                                      keep_x = FALSE,
                                      return_residual = FALSE,
                                      partial_out_y_fe = TRUE,
                                      i = NULL) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(y_name %in% names(DT))

  by <- unique(as.character(by %||% character()))
  by <- by[by %in% names(DT)]

  if (!is.null(sample_var)) {
    stopifnot(is.character(sample_var), length(sample_var) == 1L)
    if (!(sample_var %in% names(DT))) {
      stop("`sample_var` not found in `DT`: ", sample_var, call. = FALSE)
    }
  }

  by_fe <- by
  by_sp <- unique(c(sample_var, by))
  by_sp <- by_sp[by_sp %in% names(DT)]

  has_FE <- !is.null(fe_expr)
  has_X_expr <- !is.null(x_expr) && !identical(x_expr, quote(1))

  if (is.null(i)) {
    i <- seq_len(nrow(DT))
  }

  X_info <- NULL

  ## X residualization is full-sample within by_fe, not by sample half.
  if (!isTRUE(parametric) && has_X_expr && is.null(x_names)) {
    X_info <- make_X_residualized_from_FE(
      DT = DT,
      x_expr = x_expr,
      fe_expr = fe_expr,
      prefix = x_prefix,
      by = by_fe,
      weight = weight,
      fixest_opts = fixest_opts
    )

    x_names <- X_info$x_names
  }

  has_X <- !is.null(x_names) && length(x_names) > 0L

  if (isTRUE(parametric)) {

    if (has_X_expr || has_FE) {
      po <- feols_partial_out(
        DT = DT,
        y = y_name,
        rhs_expr = if (has_X_expr) x_expr else quote(1),
        fe_expr = fe_expr,
        by = by_fe,
        weight = weight,
        prefix = y_name,
        keep = "fitted",
        fixest_opts = fixest_opts
      )

      DT[i, (out_hat) := .SD[[po$fitted]], .SDcols = po$fitted]

    } else {
      leave_cluster_out_mean(
        DT = DT,
        y_name = y_name,
        cluster_name = cluster,
        by = by_sp,
        out_name = out_hat,
        weight_name = weight
      )
    }

  } else {

    if (isTRUE(partial_out_y_fe)) {
      po <- feols_partial_out(
        DT = DT,
        y = y_name,
        rhs_expr = quote(1),
        fe_expr = fe_expr,
        by = by_fe,
        weight = weight,
        prefix = y_name,
        keep = "both",
        fixest_opts = fixest_opts
      )

      y_for_rf <- po$resid
      y_fehat <- po$fitted

    } else {
      y_for_rf <- y_name
      y_fehat <- NULL
    }

    y_sp_hat <- paste0(y_name, ".sp_hat")

    if (has_X) {
      crossfit_hat(
        DT,
        i = i,
        y_name = y_for_rf,
        x_names = x_names,
        folds = foldname,
        margins = by_sp,
        weight_name = weight,
        mode = ifelse(
          (!is.null(crossfit_label) &&
             crossfit_label %in% crossfit &&
             is.null(foldname)),
          "across",
          "within"
        ),
        forest_opts = forest_opts,
        hat_suffix = ".sp_hat"
      )

      DT[i, (y_sp_hat) := .SD[[paste0(y_for_rf, ".sp_hat")]],
         .SDcols = paste0(y_for_rf, ".sp_hat")]

    } else if (has_FE && isTRUE(partial_out_y_fe)) {
      DT[i, (y_sp_hat) := 0]

    } else if (has_FE && !isTRUE(partial_out_y_fe)) {
      po_fe <- feols_partial_out(
        DT = DT,
        y = y_name,
        rhs_expr = quote(1),
        fe_expr = fe_expr,
        by = by_fe,
        weight = weight,
        prefix = paste0(y_name, "_feonly"),
        keep = "fitted",
        fixest_opts = fixest_opts
      )

      DT[i, (y_sp_hat) := .SD[[po_fe$fitted]], .SDcols = po_fe$fitted]

    } else {
      leave_cluster_out_mean(
        DT = DT,
        y_name = y_for_rf,
        cluster_name = cluster,
        by = by_sp,
        out_name = y_sp_hat,
        weight_name = weight
      )
    }

    if (isTRUE(partial_out_y_fe)) {
      DT[i, (out_hat) := .SD[[y_fehat]] + .SD[[y_sp_hat]],
         .SDcols = c(y_fehat, y_sp_hat)]
    } else {
      DT[i, (out_hat) := .SD[[y_sp_hat]], .SDcols = y_sp_hat]
    }
  }

  if (isTRUE(return_residual)) {
    if (is.null(out_resid)) {
      out_resid <- paste0(y_name, ".res")
    }

    DT[i, (out_resid) := .SD[[y_name]] - .SD[[out_hat]],
       .SDcols = c(y_name, out_hat)]
  }

  if (!isTRUE(keep_x) && !is.null(X_info$x_names)) {
    drop_x <- X_info$x_names[X_info$x_names %in% names(DT)]
    if (length(drop_x)) DT[, (drop_x) := NULL]
  }

  list(
    hat = out_hat,
    resid = if (isTRUE(return_residual)) out_resid else NULL,
    x_names = x_names,
    x_info = X_info
  )
}

make_X_residualized_from_FE <- function(DT,
                                        x_expr,
                                        fe_expr = NULL,
                                        prefix = "__x",
                                        by = NULL,
                                        weight = NULL,
                                        fixest_opts = list()) {
  stopifnot(data.table::is.data.table(DT))

  if (is.null(x_expr) || identical(x_expr, quote(1))) {
    return(list(
      x_names = NULL,
      raw_names = NULL,
      resid_names = NULL,
      has_FE = !is.null(fe_expr)
    ))
  }

  by <- unique(as.character(by %||% character()))
  by <- by[by %in% names(DT)]

  ## 1. Build raw model matrix from formula RHS.
  fml_x <- stats::as.formula(
    call("~", x_expr),
    env = parent.frame()
  )

  mm <- stats::model.matrix(fml_x, data = DT)

  ## Drop intercept if present.
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }

  if (ncol(mm) == 0L) {
    return(list(
      x_names = NULL,
      raw_names = NULL,
      resid_names = NULL,
      has_FE = !is.null(fe_expr)
    ))
  }

  raw_names <- paste0(prefix, "_raw_", make.names(colnames(mm), unique = TRUE))
  mm_dt <- data.table::as.data.table(mm)
  data.table::setnames(mm_dt, raw_names)
  DT[, (raw_names) := mm_dt]

  ## 2. If no FE, raw model-matrix columns ARE the forest features.
  if (is.null(fe_expr)) {
    return(list(
      x_names = raw_names,
      raw_names = raw_names,
      resid_names = raw_names,
      has_FE = FALSE
    ))
  }

  ## 3. If FE exist, residualize each model-matrix column from FE.
  resid_names <- paste0(prefix, "_res_", seq_along(raw_names))

  for (jj in seq_along(raw_names)) {
    out_j <- resid_names[jj]

    po <- feols_partial_out(
      DT = DT,
      y = raw_names[jj],
      rhs_expr = quote(1),
      fe_expr = fe_expr,
      by = by,
      weight = weight,
      prefix = out_j,
      keep = "resid",
      fixest_opts = fixest_opts
    )

    ## Depending on your feols_partial_out() return convention:
    ## if po$resid is the created residual column, copy/rename it.
    if (!identical(po$resid, out_j)) {
      DT[, (out_j) := get(po$resid)]
    }
  }

  ## Optional: drop raw model-matrix columns after residualizing.
  DT[, (raw_names) := NULL]

  list(
    x_names = resid_names,
    raw_names = raw_names,
    resid_names = resid_names,
    has_FE = TRUE
  )
}

##Downgrade linear specs
normalize_linear <- function(linear, J, K) {
  linear <- unique(linear)

  out <- character()

  for (ll in linear) {
    if (ll == "none") {
      out <- c(out, "none")

    } else if (ll == "Z") {
      if (K > 2L) {
        out <- c(out, "Z")
      } else {
        out <- c(out, "none")
      }

    } else if (ll == "D") {
      if (J > 2L) {
        out <- c(out, "D")
      } else {
        out <- c(out, "none")
      }

    } else if (ll == "DZ") {
      z_lin <- K > 2L
      d_lin <- J > 2L

      if (z_lin && d_lin) {
        out <- c(out, "DZ")
      } else if (z_lin && !d_lin) {
        out <- c(out, "Z")
      } else if (!z_lin && d_lin) {
        out <- c(out, "D")
      } else {
        out <- c(out, "none")
      }

    } else {
      stop("Unknown linear option: ", ll)
    }
  }

  unique(out)
}

##SCORES computation
make_scores_vec <- function(Y,
                            Z,
                            Y.hat,
                            Z.hat,
                            tau = NULL,
                            target = c("all", "overlap"),
                            target_binary=NULL,
                            target_linear=NULL,
                            z_is_linear = FALSE,
                            Z.var.hat = NULL,
                            weight = NULL,                      # <-- new
                            clip = 1e-3,
                            var_floor = 1e-6) {
  target <- match.arg(target)

  target_binary <- if (is.null(target_binary)) target else match.arg(target_binary, c("all", "overlap"))
  target_linear <- if (is.null(target_linear)) target else match.arg(target_linear, c("all", "overlap"))

  n <- length(Y)

  if (length(z_is_linear) == 1L) {
    z_is_linear <- rep(z_is_linear, n)
  }

  if (is.null(weight)) {                                        # <-- new
    weight <- rep(1, n)                                         # <-- new
  } else if (length(weight) != n) {                              # <-- new
    stop("`weight` must have length equal to `Y`.", call. = FALSE)  # <-- new
  }                                                               # <-- new

  score <- rep(NA_real_, n)

  idx_binary <- which(!z_is_linear)
  idx_linear <- which(z_is_linear)

  tau_needed <- length(idx_binary) > 0L || (target_linear == "all" && length(idx_linear) > 0L)

  if (tau_needed && is.null(tau)) {
    stop("`tau` is required for binary Z rows and for target = 'all' with linear Z rows.",
         call. = FALSE)
  }

  if (!is.null(tau) && length(tau) != n) {
    stop("`tau` must have length equal to `Y`.", call. = FALSE)
  }

  if (length(idx_binary)) {
    ii <- idx_binary

    y <- as.numeric(Y[ii])
    z <- as.numeric(Z[ii])
    m <- as.numeric(Y.hat[ii])
    e <- as.numeric(Z.hat[ii])
    t <- as.numeric(tau[ii])
    wt <- as.numeric(weight[ii])                                 # <-- new

    w <- as.numeric(z > 0.5)

    if (any(e < 0 | e > 1, na.rm = TRUE)) {
      e <- plogis(e)
    }

    if (!is.null(clip)) {
      e <- pmin(pmax(e, clip), 1 - clip)
    }

    m1 <- m + (1 - e) * t
    m0 <- m - e * t

    if (target_binary == "all") {
      score[ii] <-
        t +
        (w / e) * (y - m1) -
        ((1 - w) / (1 - e)) * (y - m0)
    } else {
      h <- e * (1 - e)
      hbar <- stats::weighted.mean(h, w = wt, na.rm = TRUE)       # <-- changed

      score[ii] <-
        (
          h * t +
            (1 - e) * w * (y - m1) -
            e * (1 - w) * (y - m0)
        ) / hbar
    }
  }

  if (length(idx_linear)) {
    ii <- idx_linear

    y <- as.numeric(Y[ii])
    z <- as.numeric(Z[ii])
    m <- as.numeric(Y.hat[ii])
    e <- as.numeric(Z.hat[ii])
    wt <- as.numeric(weight[ii])                                 # <-- new

    zres <- z - e
    yres <- y - m

    if (target_linear== "all") {
      if (is.null(Z.var.hat)) {
        stop("`Z.var.hat` is required for target = 'all' with continuous/linear Z.",
             call. = FALSE)
      }

      t <- as.numeric(tau[ii])
      v <- as.numeric(Z.var.hat[ii])
      v <- pmax(v, var_floor)

      score[ii] <- t + (zres / v) * (yres - t * zres)

    } else {
      zres_c <- zres - stats::weighted.mean(zres, w = wt, na.rm = TRUE)  # <-- changed
      yres_c <- yres - stats::weighted.mean(yres, w = wt, na.rm = TRUE)  # <-- changed

      den <- stats::weighted.mean(zres_c^2, w = wt, na.rm = TRUE)        # <-- changed

      score[ii] <- zres_c * yres_c / den
    }
  }

  score
}

make_scores <- function(DT,
                        y_name,
                        z_name,
                        y_hat_name,
                        z_hat_name,
                        tau_name,
                        target = c("all", "overlap"),
                        zvar_name = NULL,
                        weight_name = NULL,                       # <-- new
                        score_name = "scores",
                        i = NULL,
                        z_is_linear_name = "z_is_linear",
                        clip = 1e-3,
                        var_floor = 1e-6) {
  target <- match.arg(target)

  if (is.null(i)) i <- seq_len(nrow(DT))
  i <- as.integer(i)

  z_is_linear <- if (!is.null(z_is_linear_name) &&
                     z_is_linear_name %in% names(DT)) {
    as.logical(DT[[z_is_linear_name]][i])
  } else {
    rep(FALSE, length(i))
  }

  Z.var.hat <- if (!is.null(zvar_name)) {
    DT[[zvar_name]][i]
  } else {
    NULL
  }

  weight <- if (!is.null(weight_name)) DT[[weight_name]][i] else NULL   # <-- new

  score <- make_scores_vec(
    Y = DT[[y_name]][i],
    Z = DT[[z_name]][i],
    Y.hat = DT[[y_hat_name]][i],
    Z.hat = DT[[z_hat_name]][i],
    tau <- if (!is.null(tau_name)) DT[[tau_name]][i] else NULL,
    target = target,
    z_is_linear = z_is_linear,
    Z.var.hat = Z.var.hat,
    weight = weight,                                                    # <-- new
    clip = clip,
    var_floor = var_floor
  )

  if (!(score_name %in% names(DT))) {
    DT[, (score_name) := NA_real_]
  }

  DT[i, (score_name) := score]

  invisible(DT)
}

# ========= Helper: Estimate causal/regression/instrumental forests =========
fit_models <- function(DT,
                       i = NULL,
                       forest_type = c("causal", "regression"),
                       y_name,
                       x_names,
                       margins = NULL,
                       w_name = NULL,
                       zvar_name = NULL,
                       folds = NULL,
                       weight_name = NULL,
                       cluster_name = NULL,
                       forest_opts = NULL,
                       aipw.clip = 1e-3,
                       shrink = FALSE,
                       verbose = FALSE,
                       target = c("all", "overlap"),
                       target_binary=NULL,
                       target_linear=NULL,
                       z_linear_score_name = "z_use_linear_score") {

  target <- match.arg(target)
  target_binary <- if (is.null(target_binary)) target else match.arg(target_binary, c("all", "overlap"))
  target_linear <- if (is.null(target_linear)) target else match.arg(target_linear, c("all", "overlap"))

  stopifnot(data.table::is.data.table(DT))
  forest_type <- match.arg(forest_type)

  do_scores <- identical(forest_type, "causal")

  stopifnot("sample" %in% names(DT))
  stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))

  if (!is.null(folds)) {
    stopifnot(is.character(folds), length(folds) == 1L, folds %in% names(DT))
  }

  if (!is.null(aipw.clip)) {
    stopifnot(
      is.numeric(aipw.clip),
      length(aipw.clip) == 1L,
      is.finite(aipw.clip),
      aipw.clip >= 0,
      aipw.clip < 1
    )
  }

  stopifnot(is.logical(shrink), length(shrink) == 1L, !is.na(shrink))

  if (is.null(i)) i <- DT[, .I]
  i <- as.integer(i)
  n_i <- length(i)
  if (n_i == 0L) return(invisible(DT))

  if (is.null(margins)) margins <- character()
  grp_cols <- as.character(margins)

  if (is.null(forest_opts)) forest_opts <- list()

  if (!("pred"   %in% names(DT))) DT[, pred   := NA_real_]
  if (!("pred_o" %in% names(DT))) DT[, pred_o := NA_real_]

  if (do_scores && !("scores" %in% names(DT))) {
    DT[, scores := NA_real_]
  }

  if (shrink) {
    if (!("pred_var"   %in% names(DT))) DT[, pred_var   := NA_real_]
    if (!("pred_o_var" %in% names(DT))) DT[, pred_o_var := NA_real_]
  }

  y_hat <- paste0(y_name, ".hat")
  w_hat <- if (!is.null(w_name)) paste0(w_name, ".hat") else NULL

  if (forest_type == "causal") {
    stopifnot(!is.null(w_name), w_name %in% names(DT))
    stopifnot(y_hat %in% names(DT), w_hat %in% names(DT))
  }

  if (is.null(x_names) || length(x_names) == 0L) {
    stop(
      "`fit_models()` requires non-empty `x_names`. ",
      "The main X part of the formula produced no forest features.",
      call. = FALSE
    )
  }

  missing_x <- setdiff(x_names, names(DT))
  if (length(missing_x)) {
    stop("Missing `x_names` columns in `DT`: ",
         paste(missing_x, collapse = ", "),
         call. = FALSE)
  }

  X_all <- as.matrix(DT[i, ..x_names])

  rid_col <- ".__rid__"
  DT[i, (rid_col) := seq_len(n_i)]
  rid <- DT[[rid_col]]

  sample_all <- DT[["sample"]]
  folds_all  <- if (!is.null(folds)) DT[[folds]] else NULL
  wgt_all    <- if (!is.null(weight_name)) as.numeric(DT[[weight_name]]) else NULL
  cl_all     <- if (!is.null(cluster_name)) DT[[cluster_name]] else NULL

  z_is_linear_all <- if (!is.null(z_linear_score_name) &&
                         z_linear_score_name %in% names(DT)) {
    as.logical(DT[[z_linear_score_name]])
  } else if ("z_is_linear" %in% names(DT)) {
    as.logical(DT[["z_is_linear"]])
  } else {
    rep(FALSE, nrow(DT))
  }

  zvar_all <- if (!is.null(zvar_name)) {
    stopifnot(zvar_name %in% names(DT))
    as.numeric(DT[[zvar_name]])
  } else {
    NULL
  }

  if (do_scores && target_linear == "all" && any(z_is_linear_all[i], na.rm = TRUE)) {

    if (is.null(zvar_all)) {
      stop(
        "`zvar_name` is required for continuous/linear-score Z rows w hen target_linear = \"all\".",
        call. = FALSE
      )
    }

    bad_zvar <- z_is_linear_all[i] & !is.finite(zvar_all[i])

    if (any(bad_zvar, na.rm = TRUE)) {
      stop(
        "Non-finite values found in `", zvar_name,
        "` for rows that use the linear-score Z path.",
        call. = FALSE
      )
    }
  }

  y_all <- as.numeric(DT[[y_name]])

  if (forest_type == "causal") {
    w_all    <- as.numeric(DT[[w_name]])
    yhat_all <- as.numeric(DT[[y_hat]])
    what_all <- as.numeric(DT[[w_hat]])

    bad_nuis <- !is.finite(yhat_all[i]) | !is.finite(what_all[i])
    if (any(bad_nuis, na.rm = TRUE)) {
      stop(
        "Non-finite nuisance predictions in `", y_hat,
        "` or `", w_hat, "` for rows passed to `fit_models()`.",
        call. = FALSE
      )
    }
  }

  y_all <- as.numeric(DT[[y_name]])

  if (forest_type == "causal") {
    w_all    <- as.numeric(DT[[w_name]])
    yhat_all <- as.numeric(DT[[y_hat]])
    what_all <- as.numeric(DT[[w_hat]])
  }

  build_forest_idx <- function(type, idx, compute.oob.predictions = FALSE) {
    X <- X_all[rid[idx], , drop = FALSE]
    Y <- y_all[idx]
    sw <- if (is.null(wgt_all)) NULL else wgt_all[idx]
    cl <- if (is.null(cl_all)) NULL else cl_all[idx]

    if (type == "regression") {
      do.call(
        grf::regression_forest,
        c(
          list(
            X = X,
            Y = Y,
            sample.weights = sw,
            clusters = cl,
            compute.oob.predictions = compute.oob.predictions
          ),
          forest_opts
        )
      )
    } else {
      W <- w_all[idx]

      do.call(
        grf::causal_forest,
        c(
          list(
            X = X,
            Y = Y,
            W = W,
            Y.hat = yhat_all[idx],
            W.hat = what_all[idx],
            sample.weights = sw,
            clusters = cl,
            compute.oob.predictions = compute.oob.predictions
          ),
          forest_opts
        )
      )
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
        fit_k <- build_forest_idx(
          forest_type,
          idx_s[tr_pos],
          compute.oob.predictions = FALSE
        )

        out_k <- predict_with_optional_var(fit_k, Xs[te_pos, , drop = FALSE])
        p[te_pos] <- out_k$pred

        if (shrink) vv[te_pos] <- out_k$var
      }
    }

    list(pred = p, var = vv)
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
  score_buf  <- if (do_scores) rep(NA_real_, n_i) else NULL

  if (shrink) {
    pred_var_buf   <- rep(NA_real_, n_i)
    pred_o_var_buf <- rep(NA_real_, n_i)
  }

  for (g in seq_along(group_list)) {
    idx_g <- group_list[[g]]
    if (length(idx_g) == 0L) next

    idx1 <- idx_g[sample_all[idx_g] == 1L]
    idx2 <- idx_g[sample_all[idx_g] == 2L]

    fit1 <- if (length(idx1)) {
      build_forest_idx(
        forest_type,
        idx1,
        compute.oob.predictions = is.null(folds_all)
      )
    } else {
      NULL
    }

    fit2 <- if (length(idx2)) {
      build_forest_idx(
        forest_type,
        idx2,
        compute.oob.predictions = is.null(folds_all)
      )
    } else {
      NULL
    }

    if (is.null(folds_all)) {
      res1 <- if (!is.null(fit1) && length(idx1)) {
        predict_oob_with_optional_var(fit1, length(idx1))
      } else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }

      res2 <- if (!is.null(fit2) && length(idx2)) {
        predict_oob_with_optional_var(fit2, length(idx2))
      } else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }

    } else {
      res1 <- if (length(idx1)) {
        within_pred_idx(idx1)
      } else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }

      res2 <- if (length(idx2)) {
        within_pred_idx(idx2)
      } else {
        list(pred = numeric(), var = if (shrink) numeric() else NULL)
      }
    }

    res1_o <- if (!is.null(fit2) && length(idx1)) {
      predict_with_optional_var(fit2, X_all[rid[idx1], , drop = FALSE])
    } else {
      list(pred = numeric(), var = if (shrink) numeric() else NULL)
    }

    res2_o <- if (!is.null(fit1) && length(idx2)) {
      predict_with_optional_var(fit1, X_all[rid[idx2], , drop = FALSE])
    } else {
      list(pred = numeric(), var = if (shrink) numeric() else NULL)
    }

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

    if (do_scores && length(idx1)) {
      sc1 <- make_scores_vec(
        Y = y_all[idx1],
        Z = w_all[idx1],
        Y.hat = yhat_all[idx1],
        Z.hat = what_all[idx1],
        tau = p1,
        target = target,
        target_binary = target_binary,
        target_linear = target_linear,
        z_is_linear = z_is_linear_all[idx1],
        Z.var.hat = if (!is.null(zvar_all)) zvar_all[idx1] else NULL,
        weight = if (is.null(wgt_all)) NULL else wgt_all[idx1],   # <-- new
        clip = aipw.clip
      )
    }

    if (do_scores && length(idx2)) {
      sc2 <- make_scores_vec(
        Y = y_all[idx2],
        Z = w_all[idx2],
        Y.hat = yhat_all[idx2],
        Z.hat = what_all[idx2],
        tau = p2,
        target = target,
        target_binary = target_binary,
        target_linear = target_linear,
        z_is_linear = z_is_linear_all[idx2],
        Z.var.hat = if (!is.null(zvar_all)) zvar_all[idx2] else NULL,
        weight = if (is.null(wgt_all)) NULL else wgt_all[idx2],   # <-- new
        clip = aipw.clip
      )
    }

    if (length(idx1)) {
      r1 <- rid[idx1]
      pred_buf[r1]   <- p1
      pred_o_buf[r1] <- p1_o

      if (do_scores) {
        score_buf[r1] <- sc1
      }

      if (shrink) {
        pred_var_buf[r1]   <- v1
        pred_o_var_buf[r1] <- v1_o
      }
    }

    if (length(idx2)) {
      r2 <- rid[idx2]
      pred_buf[r2]   <- p2
      pred_o_buf[r2] <- p2_o

      if (do_scores) {
        score_buf[r2] <- sc2
      }

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
    if (do_scores) {
      DT[i, `:=`(
        pred = pred_buf,
        pred_var = pred_var_buf,
        pred_o = pred_o_buf,
        pred_o_var = pred_o_var_buf,
        scores = score_buf
      )]
    } else {
      DT[i, `:=`(
        pred = pred_buf,
        pred_var = pred_var_buf,
        pred_o = pred_o_buf,
        pred_o_var = pred_o_var_buf
      )]
    }
  } else {
    if (do_scores) {
      DT[i, `:=`(
        pred = pred_buf,
        pred_o = pred_o_buf,
        scores = score_buf
      )]
    } else {
      DT[i, `:=`(
        pred = pred_buf,
        pred_o = pred_o_buf
      )]
    }
  }

  DT[, (rid_col) := NULL]

  invisible(DT)
}

############
## ---------------------------------------------------------------------
## Patch: exclude FE terms nested within the clustering variable from the
## FE-rank penalty, in both fe_rank() (static/one-shot) and
## add_running_fe_rank() (training-side running search).
##
## WHY: clustering SEs on the same variable you absorb as an FE (or on a
## coarser variable that contains it, e.g. cluster = state, FE = county)
## is common and valid, but the FE's rank must NOT also be subtracted
## from `G - 1` (G = running/static cluster count) on top of that -- each
## cluster already "costs" exactly one FE parameter, so that cost is
## already implicit in G. Subtracting rank(FE) again double-penalizes df
## and drives it to <= 0 by construction, for ANY data, whenever the FE
## is (fully) nested within the cluster. This mirrors how reghdfe and
## fixest both special-case "FE nested within cluster": they explicitly
## avoid applying that second penalty.
## ---------------------------------------------------------------------

## ===== 1. fe_rank() ====================================================
## Add a `cluster_vals` argument (a vector of raw cluster values aligned
## to `idx`/`data`, NOT a column name -- every call site already has such
## a vector in scope, e.g. `clv[idx]` or a `cl` column, so no new lookups
## are needed at the call sites). When supplied, any FE term whose levels
## are each contained within a single cluster value is excluded from the
## rank penalty (contributes 0), and is still reported in `detail` for
## transparency, labeled as nested/excluded.

fe_rank <- function(data,
                    idx = NULL,
                    fe_expr,
                    weight_col = NULL,
                    verbose = FALSE,
                    ssc = NULL,
                    cluster_vals = NULL) {

  empty_result <- function() {
    list(rank = 0L, detail = data.table::data.table(
      label = character(0), K = integer(0), M = integer(0), contribution = integer(0)
    ))
  }

  if (is.null(fe_expr)) return(empty_result())

  dsub <- if (!is.null(idx)) {
    if (!length(idx)) return(empty_result())
    data[idx]
  } else {
    data
  }
  if (!nrow(dsub)) return(empty_result())

  cl_sub <- if (!is.null(cluster_vals)) {
    if (!is.null(idx)) cluster_vals[idx] else cluster_vals
  } else {
    NULL
  }
  if (!is.null(cl_sub) && length(cl_sub) != nrow(dsub)) {
    stop(
      "fe_rank: `cluster_vals` must align with `idx`/`data` -- ",
      "expected length ", nrow(dsub), ", got ", length(cl_sub), ".",
      call. = FALSE
    )
  }

  .nested_in_cluster <- function(fe_vals) {
    if (is.null(cl_sub) || !length(fe_vals)) return(FALSE)
    chk <- data.table::data.table(a = fe_vals, b = cl_sub)
    data.table::uniqueN(chk) == data.table::uniqueN(fe_vals)
  }

  ## ---- language-tree helpers -------------------------------------------

  .split_top <- function(expr, op_name) {
    op_sym <- as.name(op_name)
    if (is.call(expr) && length(expr) == 3L && identical(expr[[1]], op_sym)) {
      c(.split_top(expr[[2]], op_name), .split_top(expr[[3]], op_name))
    } else {
      list(expr)
    }
  }

  .classify_fe_term <- function(term) {
    if (is.symbol(term)) {
      return(list(kind = "fe", vars = deparse(term)))
    }
    if (is.call(term) && length(term) == 3L) {
      if (identical(term[[1]], as.name("^"))) {
        parts <- .split_top(term, "^")
        vars <- vapply(parts, deparse, character(1))
        return(list(kind = "fe_interaction", vars = vars))
      }
      if (identical(term[[1]], quote(`[`))) {
        if (is.symbol(term[[3]])) {
          return(list(kind = "slope_bundled", group = deparse(term[[2]]), xvar = deparse(term[[3]])))
        }
      }
      if (identical(term[[1]], quote(`[[`))) {
        if (is.symbol(term[[3]])) {
          return(list(kind = "slope_only", group = deparse(term[[2]]), xvar = deparse(term[[3]])))
        }
      }
    }
    list(kind = "fe", vars = paste(all.vars(term), collapse = "."), fallback = TRUE)
  }

  ## ---- union-find (path compression + union by rank) --------------------

  .uf_make <- function(n) list(parent = seq_len(n), rank = integer(n))

  .uf_find <- function(uf, x) {
    parent <- uf$parent
    root <- x
    while (parent[root] != root) root <- parent[root]
    while (parent[x] != root) {
      nxt <- parent[x]
      parent[x] <- root
      x <- nxt
    }
    uf$parent <- parent
    list(uf = uf, root = root)
  }

  .uf_union <- function(uf, x, y) {
    fx <- .uf_find(uf, x); uf <- fx$uf; rx <- fx$root
    fy <- .uf_find(uf, y); uf <- fy$uf; ry <- fy$root
    if (rx == ry) return(uf)
    if (uf$rank[rx] < uf$rank[ry]) { tmp <- rx; rx <- ry; ry <- tmp }
    uf$parent[ry] <- rx
    if (uf$rank[rx] == uf$rank[ry]) uf$rank[rx] <- uf$rank[rx] + 1L
    uf
  }

  .n_connected_components <- function(a, b) {
    fa <- factor(a); fb <- factor(b)
    n_a <- nlevels(fa); n_b <- nlevels(fb)
    ia <- as.integer(fa); ib <- as.integer(fb) + n_a
    pairs <- unique(data.table::data.table(x = ia, y = ib))
    n_nodes <- n_a + n_b
    uf <- .uf_make(n_nodes)
    if (nrow(pairs)) {
      for (k in seq_len(nrow(pairs))) uf <- .uf_union(uf, pairs$x[k], pairs$y[k])
    }
    roots <- vapply(seq_len(n_nodes), function(v) .uf_find(uf, v)$root, integer(1))
    data.table::uniqueN(roots)
  }

  ## ---- parse the FE formula into terms -----------------------------------

  rhs <- if (inherits(fe_expr, "formula")) fe_expr[[length(fe_expr)]] else fe_expr
  parsed <- lapply(.split_top(rhs, "+"), .classify_fe_term)

  cat_terms <- list()    # categorical / interacted-categorical terms, in order
  slope_rows <- list()   # locally-adjusted slope-term contributions
  nested_rows <- list()  # terms excluded because they're nested within cluster_vals

  for (p in parsed) {

    if (p$kind %in% c("fe", "fe_interaction")) {
      vars <- p$vars
      missing_vars <- setdiff(vars, names(dsub))
      if (length(missing_vars)) {
        if (verbose) message("fe_rank: skipping term, missing column(s): ",
                             paste(missing_vars, collapse = ", "))
        next
      }
      if (length(vars) == 1L) {
        combined <- dsub[[vars]]
        label <- vars
      } else {
        combined <- do.call(paste, c(lapply(vars, function(v) dsub[[v]]), sep = "\r"))
        label <- paste(vars, collapse = "^")
      }

      if (.nested_in_cluster(combined)) {
        k_here <- data.table::uniqueN(combined)
        nested_rows[[length(nested_rows) + 1L]] <- data.table::data.table(
          label = paste0(label, " (nested in cluster -- excluded from FE-rank penalty)"),
          K = k_here, M = k_here, contribution = 0L
        )
        next
      }

      cat_terms[[length(cat_terms) + 1L]] <- list(label = label, values = combined)

    } else if (p$kind %in% c("slope_bundled", "slope_only")) {
      grp <- p$group; xv <- p$xvar
      if (!(grp %in% names(dsub)) || !(xv %in% names(dsub))) {
        if (verbose) message("fe_rank: skipping slope term, missing column(s): ", grp, ", ", xv)
        next
      }
      grp_vals <- dsub[[grp]]
      k_raw <- data.table::uniqueN(grp_vals)

      if (.nested_in_cluster(grp_vals)) {
        label <- if (p$kind == "slope_bundled") {
          sprintf("%s[%s] (nested in cluster -- excluded from FE-rank penalty)", grp, xv)
        } else {
          sprintf("%s[[%s]] (nested in cluster -- excluded from FE-rank penalty)", grp, xv)
        }
        nested_rows[[length(nested_rows) + 1L]] <- data.table::data.table(
          label = label, K = k_raw, M = k_raw, contribution = 0L
        )
        next  ## also skip registering the plain-intercept companion term below
      }

      if (p$kind == "slope_bundled") {
        cat_terms[[length(cat_terms) + 1L]] <- list(label = grp, values = grp_vals)
        degenerate <- dsub[, .(const = data.table::uniqueN(get(xv)) <= 1L), by = grp]$const
        n_degenerate <- sum(degenerate)
        label <- sprintf("%s[%s] (slope, bundled intercept)", grp, xv)
      } else {
        degenerate <- dsub[, .(allzero = all(get(xv) == 0)), by = grp]$allzero
        n_degenerate <- sum(degenerate)
        label <- sprintf("%s[[%s]] (slope only)", grp, xv)
      }

      k_eff <- k_raw - n_degenerate
      slope_rows[[length(slope_rows) + 1L]] <- data.table::data.table(
        label = label, K = k_raw, M = n_degenerate, contribution = k_eff
      )
    }
  }

  ## ---- combine categorical terms: exact for the first pair, conservative
  ##      pairwise-min approximation for term 3+ --------------------------

  cat_rows <- list()
  running_total <- 0L

  if (length(cat_terms)) {
    K1 <- data.table::uniqueN(cat_terms[[1]]$values)
    cat_rows[[1]] <- data.table::data.table(label = cat_terms[[1]]$label, K = K1, M = 0L, contribution = K1)
    running_total <- K1

    if (length(cat_terms) >= 2L) {
      for (j in 2:length(cat_terms)) {
        Kj <- data.table::uniqueN(cat_terms[[j]]$values)
        pairwise_C <- vapply(
          seq_len(j - 1L),
          function(i) .n_connected_components(cat_terms[[i]]$values, cat_terms[[j]]$values),
          numeric(1)
        )
        Mj <- min(pairwise_C)
        contribution <- Kj - Mj
        running_total <- running_total + contribution
        cat_rows[[length(cat_rows) + 1L]] <- data.table::data.table(
          label = cat_terms[[j]]$label, K = Kj, M = Mj, contribution = contribution
        )
      }
    }
  }

  detail <- data.table::rbindlist(c(cat_rows, slope_rows, nested_rows), use.names = TRUE, fill = TRUE)
  total <- running_total + sum(vapply(slope_rows, function(r) r$contribution, numeric(1)))

  if (isTRUE(verbose) && nrow(detail)) {
    print(detail)
    message(sprintf("fe_rank: total rank / DoF lost to FEs = %d", as.integer(total)))
  }

  list(rank = as.integer(total), detail = detail)
}


## ===== 2. add_running_fe_rank() =========================================
## Same idea, `cluster_vals` is a vector aligned to `DT`'s rows (pass
## `dt$cl`, which the caller already builds). Any `fe_vars` entry nested
## within it is dropped BEFORE the running-count columns are built, so it
## never enters the sum at all -- contributes exactly 0 to
## `.fe_rank_running`, at every row.

add_running_fe_rank <- function(DT,
                                fe_expr,
                                by = NULL,
                                out = ".fe_rank_running",
                                conservative = TRUE,
                                cluster_vals = NULL) {
  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.character(out), length(out) == 1L)

  if (is.null(fe_expr)) {
    DT[, (out) := 0L]
    return(invisible(DT))
  }

  by <- unique(as.character(by %||% character()))
  by <- by[by %in% names(DT)]

  fe_vars <- unique(all.vars(fe_expr))
  fe_vars <- fe_vars[fe_vars %in% names(DT)]

  if (!length(fe_vars)) {
    DT[, (out) := 0L]
    return(invisible(DT))
  }

  if (!is.null(cluster_vals) && length(cluster_vals) != nrow(DT)) {
    stop("add_running_fe_rank: `cluster_vals` must have the same length as `DT`.", call. = FALSE)
  }

  if (!is.null(cluster_vals)) {
    nested <- vapply(fe_vars, function(ff) {
      chk <- data.table::data.table(a = DT[[ff]], b = cluster_vals)
      data.table::uniqueN(chk) == data.table::uniqueN(DT[[ff]])
    }, logical(1))
    fe_vars <- fe_vars[!nested]
  }

  if (!length(fe_vars)) {
    DT[, (out) := 0L]
    return(invisible(DT))
  }

  ## Running count of unique observed FE values within each by-cell.
  ## Assumes DT is already sorted in the desired cutoff order.
  rank_cols <- character(0)

  for (ff in fe_vars) {
    cc <- paste0(".__rank_seen__", make.names(ff), "__")

    if (length(by)) {
      DT[
        ,
        (cc) := as.integer(!duplicated(.SD[[ff]])),
        by = by,
        .SDcols = ff
      ]
    } else {
      DT[, (cc) := as.integer(!duplicated(.SD[[ff]])), .SDcols = ff]
    }

    if (length(by)) {
      DT[, (cc) := cumsum(get(cc)), by = by]
    } else {
      DT[, (cc) := cumsum(get(cc))]
    }

    rank_cols <- c(rank_cols, cc)
  }

  if (conservative) {
    DT[, (out) := rowSums(.SD), .SDcols = rank_cols]
  } else {
    DT[
      ,
      (out) := pmax(0L, rowSums(.SD) - length(rank_cols)),
      .SDcols = rank_cols
    ]
  }

  DT[, (rank_cols) := NULL]

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
    verbose = FALSE,
    screen = c("stepdown", "none", "minimum", "negative", "nonpositive", "fgk_relevant"),
    alpha = 0.05,
    fe_expr = NULL,
    fe_rank_adj = !is.null(fe_expr),
    fe_rank_conservative = TRUE
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

  if (sample_col %in% adaptive_dims) {
    stop(
      "`sample` cannot appear in both `pool` and `select`.\n",
      "This would adaptively choose between sample-pooled and sample-separated ",
      "testing using both sample halves. Put `sample` in either `pool` or `select`, ",
      "but not both."
    )
  }

  design_score <- function(fit, pool_eff, select_eff, design_row) {
    tr <- data.table::as.data.table(fit$results)[train == TRUE & relevant == 1L]

    if (nrow(tr) == 0L) return(Inf)
    if (!("p.raw" %in% names(tr))) return(Inf)

    pp <- tr$p.raw
    pp <- pp[is.finite(pp)]

    if (!length(pp)) return(Inf)

    min(stats::p.adjust(pp, method = "holm"), na.rm = TRUE)
  }

  run_core <- function(pool_eff, select_eff) {
    forest_test_core(
      data = data,
      cluster = cluster,
      weight = weight,
      sample = sample,
      pred = pred,
      pred_o = pred_o,
      scores = scores,
      x_names = x_names,
      minsize = minsize,
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
      fe_rank_conservative = fe_rank_conservative
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

crv1_mean <- function(score,
                      w = NULL,
                      cl = NULL,
                      rank_adj = 0) {
  score <- as.numeric(score)
  n0 <- length(score)

  if (n0 == 0L) {
    return(list(
      coef = NA_real_, se = NA_real_, t = NA_real_,
      N = 0L, G = 0L, rank_adj = rank_adj, df = NA_real_
    ))
  }

  if (is.null(w)) {
    w <- rep(1.0, n0)
  } else {
    w <- as.numeric(w)
  }

  if (length(w) != n0) {
    stop("`w` must be NULL or have the same length as `score`.", call. = FALSE)
  }

  if (is.null(cl)) {
    cl <- seq_len(n0)
  } else {
    if (length(cl) != n0) {
      stop("`cl` must be NULL or have the same length as `score`.", call. = FALSE)
    }
  }

  ok <- is.finite(score) & is.finite(w) & w >= 0 & !is.na(cl)

  if (!any(ok)) {
    return(list(
      coef = NA_real_, se = NA_real_, t = NA_real_,
      N = 0L, G = 0L, rank_adj = rank_adj, df = NA_real_
    ))
  }

  score <- score[ok]
  w <- w[ok]
  cl <- cl[ok]

  cl <- as.integer(factor(cl, exclude = NULL))

  rank_adj <- as.numeric(rank_adj)
  if (length(rank_adj) != 1L || !is.finite(rank_adj) || rank_adj < 0) {
    stop("`rank_adj` must be a single nonnegative finite number.", call. = FALSE)
  }

  dt0 <- data.table::data.table(score = score, w = w, cl = cl)

  gb <- dt0[, .(
    U = sum(w * score),
    W = sum(w)
  ), by = cl]

  G <- nrow(gb)
  U <- sum(gb$U)
  W <- sum(gb$W)

  if (!is.finite(W) || W <= 0) {
    return(list(
      coef = NA_real_, se = NA_real_, t = NA_real_,
      N = length(score), G = G, rank_adj = rank_adj, df = NA_real_
    ))
  }

  theta <- U / W
  ug <- gb$U - theta * gb$W

  df <- G - 1 - rank_adj

  if (G < 2L || df <= 0 || !is.finite(df)) {
    se <- NA_real_
  } else {
    se <- sqrt((G / df) * sum(ug^2)) / abs(W)
    if (!is.finite(se)) se <- NA_real_
  }

  list(
    coef = theta,
    se = se,
    t = theta / se,
    N = length(score),
    G = G,
    rank_adj = rank_adj,
    df = df
  )
}


forest_test_core <- function(
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
    verbose = FALSE,
    screen = c("stepdown", "none", "minimum", "negative", "nonpositive", "fgk_relevant"),
    alpha = 0.05,
    fe_expr = NULL,
    fe_rank_adj = !is.null(fe_expr),
    fe_rank_conservative = TRUE
) {

  select_keep_ids <- function(cand, choose_by, id_by, alpha) {
    choose_by <- intersect(choose_by, names(cand))
    id_by <- intersect(id_by, names(cand))

    if (length(id_by) == 0L) {
      stop("No valid id columns found in cand.")
    }

    if (length(choose_by) == 0L) {
      out <- select_groups_screen(cand, alpha = alpha)
      return(unique(out[, .SD, .SDcols = id_by]))
    }

    ## Columns in choose_by are excluded from .SD by data.table.
    ## Therefore only ask .SD for id columns that are not also by-columns.
    id_nonby <- setdiff(id_by, choose_by)

    out <- cand[
      ,
      {
        sel <- select_groups_screen(.SD, alpha = alpha)

        if (length(id_nonby) > 0L) {
          sel[, .SD, .SDcols = id_nonby]
        } else {
          ## Needed when all id_by columns are in choose_by.
          ## Return one row per selected candidate; data.table will add choose_by.
          data.table::data.table(.selected_row__ = seq_len(nrow(sel)))
        }
      },
      by = choose_by
    ]

    if (".selected_row__" %in% names(out)) {
      out[, .selected_row__ := NULL]
    }

    unique(out[, .SD, .SDcols = id_by])
  }

  stopifnot(data.table::is.data.table(data))
  screen <- match.arg(screen)

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

  pool_sample   <- sample_col %in% pool
  select_sample <- sample_col %in% select

  if (pool_sample && select_sample) {
    stop(
      "forest_test_core() does not accept sample in both `pool` and `select`. ",
      "Use forest_test(), which evaluates sample-pooled and sample-separated ",
      "as separate adaptive designs."
    )
  }

  pool_margins   <- intersect(pool, margins)
  select_margins <- intersect(select, margins)

  overlap_margins <- intersect(pool_margins, select_margins)

  if (length(overlap_margins) > 0L) {
    stop(
      "forest_test_core() does not allow any margin or sample in both `pool` and `select`: ",
      paste(overlap_margins, collapse = ", ")
    )
  }

  ## Special valid version of pool="sample" with margin selection:
  ## select margins separately by training sample, test separately by holdout sample,
  ## then raw-pool only margin keys selected by both sample directions.
  sample_pool_after_margin_select <-
    pool_sample &&
    length(select_margins) > 0L &&
    !select_sample

  pool_before_margins <- pool_margins
  separate_select_margins <- select_margins
  cell_cols <- setdiff(margins, pool_before_margins)
  adjust_cols <- setdiff(cell_cols, select_margins)

  test_by <- unique(c(
    adjust_cols,
    separate_select_margins,
    if (!pool_sample || sample_pool_after_margin_select) sample_col else character()
  ))

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

  fe_vars_needed <- if (!is.null(fe_expr)) unique(all.vars(fe_expr)) else character(0)
  fe_vars_needed <- fe_vars_needed[fe_vars_needed %in% names(data)]
  fe_dt_full <- if (length(fe_vars_needed)) data[, ..fe_vars_needed] else NULL


  select_groups_screen <- function(dt_grp, t_col = "sel_t", alpha = 0.05) {
    tt <- dt_grp[[t_col]]
    ok <- is.finite(tt)
    if (!any(ok)) return(dt_grp[0])

    d <- data.table::copy(dt_grp[ok])
    data.table::setorderv(d, t_col)

    if (screen == "none") return(d)
    if (screen == "minimum") return(d[1L])

    if (screen == "negative") {
      thr <- stats::qnorm(alpha / nrow(d))
      keep <- d[get(t_col) <= thr]
      if (nrow(keep) == 0L) keep <- d[1L]
      return(keep)
    }

    if (screen == "nonpositive") {
      thr <- stats::qnorm(1 - alpha / nrow(d))
      keep <- d[get(t_col) < thr]
      if (nrow(keep) == 0L) keep <- d[1L]
      return(keep)
    }

    if (screen == "stepdown") {
      keep_n <- 1L
      if (nrow(d) >= 2L) {
        p <- stats::pnorm(d[[t_col]])
        for (r in 2:nrow(d)) {
          p_top <- p[seq_len(r)]
          p_holm <- stats::p.adjust(p_top, method = "holm")
          if (p_holm[r] <= alpha) keep_n <- r else break
        }
      }
      return(d[seq_len(keep_n)])
    }

    if (screen == "fgk_relevant") {
      thr <- stats::qnorm(1 - alpha / nrow(d))
      keep <- d[get(t_col) < thr]

      ## Fallback at the level at which select_groups_screen()
      ## is called. Therefore if select includes yval/dval/sample,
      ## this fallback is global over those dimensions.
      if (nrow(keep) == 0L) keep <- d[1L]

      return(keep)
    }

    d[1L]
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
    gid <- data.table::frankv(data[, .SD, .SDcols = cell_cols], ties.method = "dense")
    idx_list <- split(seq_len(n), gid)
    keys_list <- lapply(idx_list, function(idx) data[idx[1L], .SD, .SDcols = cell_cols])
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

    if (!is.null(fe_dt_full)) {
      dt[, (fe_vars_needed) := fe_dt_full[idx]]
    }

    data.table::setorder(dt, sample, pred)

    if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
      add_running_fe_rank(
        DT = dt,
        fe_expr = fe_expr,
        by = "sample",
        out = ".fe_rank_running",
        conservative = fe_rank_conservative,
        cluster_vals = dt$cl
      )
    } else {
      dt[, .fe_rank_running := 0L]
    }

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

    dt[, `:=`(
      TA2 = cumsum(dTA2),
      TB2 = cumsum(dTB2),
      TAB = cumsum(dTAB)
    ), by = sample]

    dt[, G := cumsum(!duplicated(cl)), by = sample]

    dt[, sumS2 := (TA2 - 2 * m * TAB + (m^2) * TB2) / (SW^2)]

    tol_var <- 1e-12

    if (any(dt$sumS2 < -tol_var & dt$G >= minsize, na.rm = TRUE)) {
      warning(
        "Negative cluster variance encountered in an eligible-size prefix. ",
        "This may indicate numerical instability or problematic weights."
      )
    }

    dt[
      is.finite(sumS2) & sumS2 < 0 & sumS2 >= -tol_var,
      sumS2 := 0
    ]

    dt[, df := G - 1 - .fe_rank_running]

    dt[, se := NA_real_]

    dt[
      G >= 2L &
        df > 0 &
        SW > 0 &
        is.finite(sumS2) &
        sumS2 >= 0,
      se := sqrt((G / df) * sumS2)
    ]

    dt[
      !is.finite(se) | se <= 0,
      se := NA_real_
    ]

    dt[, t_stat := NA_real_]
    dt[
      is.finite(se),
      t_stat := m / se
    ]

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

    elig <- (dt$pred >= dt$tau_constraint) &
      (dt$G >= minsize) &
      dt$cand_search &
      is.finite(dt$t_stat)

    res <- dt[elig, .SD[which.min(t_stat)], by = sample,
              .SDcols = c("G", "N", "m", "se", "t_stat", "pred",
                          ".fe_rank_running", "df")]
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

        data.table::data.table(
          tau = pred[idx_all],
          t = t_stat[idx_all]
        )
      }, by = "sample"]

      if (!is.null(key_dt) && ncol(key_dt) > 0L) {
        for (cc in names(key_dt)) grid_out[, (cc) := key_dt[[cc]][1]]
        data.table::setcolorder(grid_out, c(names(key_dt), "sample", "t", "tau"))
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
      fe_rank = res$.fe_rank_running,
      df = res$df,
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

    list(
      train = train_out,
      grid = grid_out,
      idx_test = idx_test_by_sample
    )
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
  train_all[, train_row_id := .I]



  grid_out <- if (store_grid) {
    data.table::rbindlist(grid_list, use.names = TRUE, fill = TRUE)
  } else {
    NULL
  }


  selected_train <- train_all
  has_selection <- (length(select_margins) > 0L) || select_sample

  if (has_selection) {

     if (sample_pool_after_margin_select) {

      ## Do NOT select margins using sample-pooled training statistics.
      ## Select separately by training sample.
      cand_by <- unique(c(
        adjust_cols,
        sample_col,
        select_margins,
        "cell_id"
      ))

      cand <- train_all[
        ,
        .(sel_t = mean(t, na.rm = TRUE)),
        by = cand_by
      ]

      choose_by <- unique(c(
        adjust_cols,
        sample_col
      ))

      id_by <- cand_by

      keep_ids <- select_keep_ids(
        cand = cand,
        choose_by = choose_by,
        id_by = id_by,
        alpha = alpha
      )

      selected_train <- merge(
        train_all,
        unique(keep_ids),
        by = id_by,
        all = FALSE,
        sort = FALSE
      )

    } else {

      agg_by <- unique(c(
        adjust_cols,
        sample_col,
        select_margins,
        "cell_id"
      ))

      cand <- train_all[
        ,
        .(sel_t = mean(t, na.rm = TRUE)),
        by = agg_by
      ]

      choose_by <- c(
        adjust_cols,
        if (!select_sample && !pool_sample) sample_col else character()
      )

      id_by <- unique(c(
        adjust_cols,
        sample_col,
        select_margins,
        "cell_id"
      ))

      if (length(choose_by) == 0L) {
        keep_ids <- select_groups_screen(
          cand,
          alpha = alpha
        )[, .SD, .SDcols = id_by]
      } else {
        keep_ids <- cand[
          ,
          select_groups_screen(.SD, alpha = alpha),
          by = choose_by
        ][, .SD, .SDcols = id_by]
      }

      selected_train <- merge(
        train_all,
        unique(keep_ids),
        by = id_by,
        all = FALSE,
        sort = FALSE
      )
    }
  }

  sample_design <- if (sample_pool_after_margin_select) {
    "pooled_matching_after_sample_specific_margin_selection"
  } else if (pool_sample) {
    "pooled"
  } else if (select_sample) {
    "selected_then_separate"
  } else {
    "separate"
  }

  train_out <- data.table::copy(train_all)
  train_out[train_row_id %in% selected_train$train_row_id, relevant := 1L]

  train_out[, c("cell_id", "train_row_id") := NULL]

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

  selected_test_keys <- NULL

  if (length(separate_select_margins) > 0L) {
    key_cols <- intersect(
      unique(c(
        adjust_cols,
        if (!pool_sample || sample_pool_after_margin_select) sample_col else character(),
        separate_select_margins
      )),
      names(selected_train)
    )

    selected_test_keys <- unique(selected_train[, .SD, .SDcols = key_cols])

    ## selected_train$sample is the training sample.
    ## The final test row lives in the opposite holdout sample.
    if (sample_col %in% names(selected_test_keys)) {
      selected_test_keys[
        ,
        (sample_col) := data.table::fifelse(
          get(sample_col) == 1L,
          2L,
          1L
        )
      ]
    }

  } else {
    selected_test_keys <- NULL
  }

  dt_test <- data.table::data.table(
    rowid = seq_len(n),
    score = scorev,
    w = wv,
    cl = clv
  )

  if (length(test_by) > 0L) {
    dt_test[, (test_by) := data[, .SD, .SDcols = test_by]]
  }

  dt_test <- dt_test[in_test]

  ## First build the valid final tests.
  ## In sample_pool_after_margin_select mode this is intentionally sample-specific.
  if (length(test_by) == 0L) {
    rank_adj <- if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
      fe_rank(
        data = data,
        idx = dt_test$rowid,
        fe_expr = fe_expr,
        weight_col = weight_col,
        cluster_vals = clv
      )$rank
    } else {
      0L
    }

    o <- crv1_mean(
      score = dt_test$score,
      w = dt_test$w,
      cl = dt_test$cl,
      rank_adj = rank_adj
    )

    test_out <- data.table::data.table(
      train = FALSE,
      relevant = 1L,
      sample = NA_integer_,
      G = o$G,
      N = o$N,
      fe_rank = rank_adj,
      df = o$df,
      coef = o$coef,
      stderr = o$se,
      t = o$t,
      tau_cutoff = NA_real_,
      p.raw = stats::pnorm(o$t)
    )

    if (!is.null(selected_test_keys) && nrow(selected_test_keys) == 1L) {
      for (cc in separate_select_margins) {
        test_out[, (cc) := selected_test_keys[[cc]][1L]]
      }
    } else if (length(separate_select_margins) > 0L) {
      for (cc in separate_select_margins) test_out[, (cc) := NA]
    }

  } else {
    test_out <- dt_test[, {
      rank_adj <- if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
        fe_rank(
          data = data,
          idx = rowid,
          fe_expr = fe_expr,
          weight_col = weight_col,
          cluster_vals = clv
        )$rank
      } else {
        0L
      }

      o <- crv1_mean(
        score = score,
        w = w,
        cl = cl,
        rank_adj = rank_adj
      )

      data.table::data.table(
        train = FALSE,
        relevant = 1L,
        G = o$G,
        N = o$N,
        fe_rank = rank_adj,
        df = o$df,
        coef = o$coef,
        stderr = o$se,
        t = o$t,
        tau_cutoff = NA_real_,
        p.raw = stats::pnorm(o$t)
      )
    }, by = test_by]

    if (!(sample_col %in% names(test_out))) {
      test_out[, (sample_col) := NA_integer_]
    }

    if (!is.null(selected_test_keys)) {
      merge_by <- intersect(
        unique(c(
          adjust_cols,
          if (!pool_sample || sample_pool_after_margin_select) sample_col else character(),
          separate_select_margins
        )),
        names(test_out)
      )
      merge_by <- intersect(merge_by, names(selected_test_keys))

      if (length(merge_by) == 0L) {
        if (nrow(selected_test_keys) == 1L) {
          for (cc in separate_select_margins) {
            test_out[, (cc) := selected_test_keys[[cc]][1L]]
          }
        } else {
          for (cc in separate_select_margins) test_out[, (cc) := NA]
        }
      } else {
        selected_test_keys_u <- unique(selected_test_keys)

        test_out <- merge(
          test_out,
          selected_test_keys_u,
          by = merge_by,
          all.x = TRUE,
          sort = FALSE
        )
      }
    }
  }

  ## Now implement option 2:
  ## if both holdout sample directions selected the same non-sample key,
  ## raw-pool those observations; otherwise keep sample-specific rows.
  if (sample_pool_after_margin_select) {
    pool_key <- setdiff(test_by, sample_col)
    pool_key <- intersect(pool_key, names(test_out))

    if (length(pool_key) > 0L) {
      common_keys <- test_out[
        ,
        .(n_sample = data.table::uniqueN(get(sample_col))),
        by = pool_key
      ][
        n_sample == 2L,
        .SD,
        .SDcols = pool_key
      ]

      if (nrow(common_keys) > 0L) {
        dt_common <- merge(
          dt_test,
          common_keys,
          by = pool_key,
          all = FALSE,
          sort = FALSE
        )

        test_out_pool <- dt_common[, {
          rank_adj <- if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
            fe_rank(
              data = data,
              idx = rowid,
              fe_expr = fe_expr,
              weight_col = weight_col,
              cluster_vals=clv
            )$rank
          } else {
            0L
          }

          o <- crv1_mean(
            score = score,
            w = w,
            cl = cl,
            rank_adj = rank_adj
          )

          data.table::data.table(
            train = FALSE,
            relevant = 1L,
            G = o$G,
            N = o$N,
            fe_rank = rank_adj,
            df = o$df,
            coef = o$coef,
            stderr = o$se,
            t = o$t,
            tau_cutoff = NA_real_,
            p.raw = stats::pnorm(o$t)
          )
        }, by = pool_key]

        test_out_pool[, (sample_col) := NA_integer_]


        test_out_single <- test_out[
          !common_keys,
          on = pool_key
        ]

        test_out <- data.table::rbindlist(
          list(test_out_pool, test_out_single),
          use.names = TRUE,
          fill = TRUE
        )

      }

    } else {
      ## Degenerate case: no non-sample key.
      ## There is nothing meaningful to selectively match on.
      test_out[, sample_pool_status := "separate_unmatched"]
    }
  }

  results_out <- data.table::rbindlist(
    list(train_out, test_out),
    use.names = TRUE,
    fill = TRUE
  )

  global <- global_means_crv1(
    data = data,
    scorev = scorev,
    wv = wv,
    clv = clv,
    by_cols = test_by,
    sample_col = sample_col,
    crv1_mean_fun = crv1_mean,
    fe_expr = fe_expr,
    fe_rank_adj = fe_rank_adj,
    weight_col = weight_col
  )

  Xmeans <- Xmeans_all <- XSD <- NULL

  if (!is.null(x_names)) {
    cols_need <- unique(c(test_by, x_names))
    df_all <- data[, .SD, .SDcols = cols_need]
    df_tst <- data[in_test, .SD, .SDcols = cols_need]

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
      separate_select_margins,
      if (!pool_sample || sample_pool_after_margin_select) sample_col else character()
    )

    denom_by <- intersect(denom_by, names(data))
    cols_need <- intersect(unique(c(denom_by, pool_margins)), names(data))

    DT_all <- data[, .SD, .SDcols = cols_need]
    DT_tst <- data[in_test, .SD, .SDcols = cols_need]

    make_share_combo <- function(DTsub, wsub) {
      DTsub <- data.table::as.data.table(DTsub)
      DTsub[, w := wsub]

      by_num <- intersect(unique(c(denom_by, pool_margins)), names(DTsub))

      if (length(by_num) == 0L) {
        num <- DTsub[, .(w_sum = sum(w))]
      } else {
        num <- DTsub[, .(w_sum = sum(w)), by = by_num]
      }

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

    sh <- make_share_combo(DT_tst, wv[in_test])
    sh_all <- make_share_combo(DT_all, wv)

    data.table::setnames(sh, "share", "share")
    data.table::setnames(sh_all, "share", "share_all")

    by_merge <- intersect(
      unique(c(denom_by, pool_margins)),
      intersect(names(sh), names(sh_all))
    )

    if (length(by_merge) == 0L) {
      shares <- cbind(sh, share_all = sh_all$share_all[1L])
    } else {
      shares <- merge(sh, sh_all, by = by_merge, all = TRUE)
    }

    shares[is.na(share), share := 0.0]
    shares[is.na(share_all), share_all := 0.0]

    data.table::setcolorder(shares, c(denom_by, pool_margins, "share", "share_all"))
  }

  out <- list(
    results = results_out,
    grid = grid_out,
    global = global,
    Xmeans = Xmeans,
    Xmeans_all = Xmeans_all,
    XSD = XSD,
    sample_design = sample_design,
    sample_pool_after_margin_select = sample_pool_after_margin_select
  )

  if (!is.null(shares)) out$shares <- shares
  if (!store_grid) out$grid <- NULL

  out
}

##USED by overlap and forest_test to estimate global means.
global_means_crv1 <- function(
    data,
    scorev,
    wv,
    clv,
    by_cols = character(),
    sample_col = "sample",
    crv1_mean_fun,
    fe_expr = NULL,
    fe_rank_adj = !is.null(fe_expr),
    weight_col = NULL
) {
  stopifnot(data.table::is.data.table(data))
  stopifnot(
    is.numeric(scorev),
    is.numeric(wv),
    length(scorev) == length(wv),
    length(wv) == length(clv),
    length(scorev) == nrow(data)
  )
  stopifnot(is.character(by_cols))
  stopifnot(is.character(sample_col), length(sample_col) == 1L)
  stopifnot(is.function(crv1_mean_fun))

  by_cols <- unique(by_cols)

  if (length(by_cols) > 0L) {
    missing <- setdiff(by_cols, names(data))
    if (length(missing)) {
      stop(
        "global_means_crv1: missing grouping columns in data: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.null(weight_col)) {
    stopifnot(is.character(weight_col), length(weight_col) == 1L)
    stopifnot(weight_col %in% names(data))
  }

  dt_all <- data.table::data.table(
    rowid = seq_len(nrow(data)),
    score = scorev,
    w = wv,
    cl = clv
  )

  if (length(by_cols) > 0L) {
    dt_all[, (by_cols) := data[, .SD, .SDcols = by_cols]]
  }

    rank_for_rows <- function(rowid,cl_vals) {
      if (isTRUE(fe_rank_adj) && !is.null(fe_expr)) {
        fe_rank(data = data, idx = rowid, fe_expr = fe_expr, weight_col = weight_col,cluster_vals=cl_vals)$rank
      } else {
        0L
      }
    }

  if (length(by_cols) == 0L) {
    rank_adj <- rank_for_rows(dt_all$rowid,clv)

    o <- crv1_mean_fun(
      score = dt_all$score,
      w = dt_all$w,
      cl = dt_all$cl,
      rank_adj = rank_adj
    )

    global_dt <- data.table::data.table(
      train = FALSE,
      sample = NA_integer_,
      G = o$G,
      N = o$N,
      fe_rank = rank_adj,
      df = o$df,
      coef = o$coef,
      stderr = o$se,
      t = o$t,
      tau_cutoff = NA_real_,
      p.raw = stats::pnorm(o$t)
    )

  } else {
    global_dt <- dt_all[, {
      rank_adj <- rank_for_rows(rowid,clv)

      o <- crv1_mean_fun(
        score = score,
        w = w,
        cl = cl,
        rank_adj = rank_adj
      )

      data.table::data.table(
        train = FALSE,
        G = o$G,
        N = o$N,
        fe_rank = rank_adj,
        df = o$df,
        coef = o$coef,
        stderr = o$se,
        t = o$t,
        tau_cutoff = NA_real_,
        p.raw = stats::pnorm(o$t)
      )
    }, by = by_cols]

    if (!(sample_col %in% names(global_dt))) {
      global_dt[, (sample_col) := NA_integer_]
    }
  }

  global_dt
}



#######ONE SIDED NONCOMPLIANCE HELPERS #######
test_one_sided_noncompliance <- function(
    data,
    D,
    Z,
    zmargin_var = NULL
) {
  dt <- data.table::as.data.table(data)

  stopifnot(is.character(D), length(D) == 1, D %in% names(dt))
  stopifnot(is.character(Z), length(Z) == 1, Z %in% names(dt))

  has_zmargin <- !is.null(zmargin_var)
  if (has_zmargin) {
    stopifnot(is.character(zmargin_var), length(zmargin_var) == 1, zmargin_var %in% names(dt))
  }

  Dcol <- D
  Zcol <- Z
  Zmcol <- zmargin_var

  if (!is.numeric(dt[[Dcol]]) && !is.integer(dt[[Dcol]])) {
    stop("D must be numeric/integer ordered treatment values.")
  }

  zvals <- sort(unique(dt[[Zcol]][!is.na(dt[[Zcol]])]))
  if (!all(zvals %in% c(0, 1))) {
    stop("Z must already be binary 0/1 within each zmargin slice.")
  }

  d_support <- sort(unique(dt[[Dcol]][!is.na(dt[[Dcol]])]))
  if (length(d_support) < 2L) stop("D must have at least two support values.")
  dmargin_vals <- d_support[-1L]

  keep_vars <- unique(c(Dcol, Zcol, if (has_zmargin) Zmcol else character(0)))
  dt2 <- dt[stats::complete.cases(dt[, ..keep_vars])]

  if (!has_zmargin) {
    z_support <- sort(unique(dt2[[Zcol]]))
    if (!all(z_support %in% c(0, 1)) || length(z_support) != 2L) {
      stop("When zmargin_var is NULL, Z must be binary 0/1.")
    }
  }

  d_index <- data.table::data.table(dmargin = dmargin_vals, dummy__ = 1L)
  dt2[, dummy__ := 1L]
  tmp <- d_index[dt2, on = "dummy__", allow.cartesian = TRUE][, dummy__ := NULL]
  dt2[, dummy__ := NULL]

  tmp[, Z_thr__ := as.integer(get(Zcol))]
  tmp[, D_thr__ := as.integer(get(Dcol) >= dmargin)]

  by_threshold <- if (has_zmargin) {
    tmp[, zmargin := get(Zmcol)]
    c("zmargin", "dmargin")
  } else {
    "dmargin"
  }

  core <- function(.SD) {
    z <- .SD[["Z_thr__"]]
    d <- .SD[["D_thr__"]]

    n_z0 <- sum(z == 0L)
    n_z1 <- sum(z == 1L)

    p_always <- if (n_z0 > 0L) mean(d[z == 0L] == 1L) else NA_real_
    p_never  <- if (n_z1 > 0L) mean(d[z == 1L] == 0L) else NA_real_

    only_always <- is.finite(p_always) && is.finite(p_never) && (p_always > 0) && (p_never == 0)
    only_never  <- is.finite(p_always) && is.finite(p_never) && (p_never > 0) && (p_always == 0)

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

    list(
      n_z0 = n_z0,
      n_z1 = n_z1,
      p_always = p_always,
      p_never = p_never,
      one_sided = one_sided,
      direction = direction
    )
  }

  threshold_res <- tmp[, core(.SD), by = by_threshold]

  if (has_zmargin) {
    margin_cells <- unique(tmp[, .(zmargin)])
    margin_cells[, dummy__ := 1L]
    exact_index <- data.table::data.table(dval = d_support, dummy__ = 1L)
    exact_res <- exact_index[margin_cells, on = "dummy__", allow.cartesian = TRUE][, dummy__ := NULL]

    d_map <- data.table::data.table(
      dval = d_support,
      next_d = c(d_support[-1L], NA)
    )
    exact_res <- d_map[exact_res, on = "dval"]

    below <- threshold_res[
      direction == "Never-takers only",
      .(zmargin, dval = dmargin, trivial_from_below = TRUE)
    ]

    above <- threshold_res[
      direction == "Always-takers only",
      .(zmargin, next_d = dmargin, trivial_from_above = TRUE)
    ]

    exact_res <- merge(exact_res, below, by = c("zmargin", "dval"), all.x = TRUE)
    exact_res <- merge(exact_res, above, by = c("zmargin", "next_d"), all.x = TRUE)
  } else {
    exact_res <- data.table::data.table(dval = d_support)

    d_map <- data.table::data.table(
      dval = d_support,
      next_d = c(d_support[-1L], NA)
    )
    exact_res <- d_map[exact_res, on = "dval"]

    below <- threshold_res[
      direction == "Never-takers only",
      .(dval = dmargin, trivial_from_below = TRUE)
    ]

    above <- threshold_res[
      direction == "Always-takers only",
      .(next_d = dmargin, trivial_from_above = TRUE)
    ]

    exact_res <- merge(exact_res, below, by = "dval", all.x = TRUE)
    exact_res <- merge(exact_res, above, by = "next_d", all.x = TRUE)
  }

  exact_res[is.na(trivial_from_below), trivial_from_below := FALSE]
  exact_res[is.na(trivial_from_above), trivial_from_above := FALSE]
  exact_res[, trivial_exact := trivial_from_below | trivial_from_above]
  exact_res[, next_d := NULL]

  if (length(d_support) == 2L && "dmargin" %in% names(threshold_res)) {
    threshold_res[, dmargin := NULL]
  }

  os <- threshold_res[one_sided == TRUE]
  if (nrow(os) != 0L) {
    message("One-sided noncompliance detected in the following threshold-margin cells:")
    print(os)
  }

  invisible(list(
    threshold = threshold_res,
    exact = exact_res
  ))
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

