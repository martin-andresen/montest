

##################################################### HELPER FUNCTIONS ######################################

# THIS FUNCTION ESTIMATES A TREE USING CART, TEST IN EACH NODE, DETERMINE relevance and perform the test in the other sample
esttree=function(data,testsample,cp,maxrankcp,alpha,prune,minsize,preselect){
  tree=rpart(scores~.,data=as.data.frame(data[sample==testsample,!c("id_","sample")]),method="anova",cp=cp,minbucket=minsize,weights=weight) # run tree with transformed outcome
  
  maxrankcp=min(maxrankcp, length(tree$cp[, 4]))
  maxcp=tree$cp[maxrankcp, 1]
  
  if (prune==TRUE){ #prune the tree based on cross-validation
    opcpid=which.min(tree$cp[, 4])
    opcp=tree$cp[opcpid, 1]
    tree=prune(tree, cp = max(maxcp,opcp))
  }      else { # prune the tree based on maximum complexity parameter only (no cross-validation)
    tree=prune(tree, cp = maxcp)
  }
  data[,leaf:=rpart.predict.leaves(tree,newdata=data,type = "where")] # create dummies for leaves (based on predictions))
  
  ##Drop leaves with no variation in scores
  if (sum(data[sample==1,sd(scores)==0,by=leaf][,V1])>0) data[sample==1&(leaf %in% data[sample==1,sd(scores)==0,by=leaf][V1==TRUE,leaf]),leaf:=NA]
  if (sum(data[sample==2,sd(scores)==0,by=leaf][,V1])>0) data[sample==2&(leaf %in% data[sample==2,sd(scores)==0,by=leaf][V1==TRUE,leaf]),leaf:=NA]
  
  data[,sample:=ifelse(sample==testsample,"train","est")]
  if (is.null(cluster)==FALSE) {
    data[,G:=cumsum(!duplicated(get(cluster))),by=.(leaf,sample)]
  }else {
    data[,G:=.N,by=.(leaf,sample)]
  }
  
  res=data[is.na(leaf)==FALSE,data.table(cbind(G=max(G),N=.N,feols(scores~1,data=.SD,vcov=clust,weights=weight)$coeftable)),by=.(leaf,sample)]
  colnames(res)=c("leaf","sample","G","N","coef","stderr","t","p")
  res=res[,p:=NULL]
  res=dcast(res,leaf~sample,value.var=c("G","N","coef","stderr","t"))
  if (preselect!="none"&nrow(res)>1) {
    if (preselect=="nonpositive") {res[,relevant:=ifelse(t_train>=qnorm(1-alpha/nrow(res))&(t_train>min(t_train)),0,1)]}
    else if (preselect=="negative") {res[,relevant:=ifelse(t_train>=qnorm(alpha/nrow(res))&(t_train>min(t_train)),0,1)]}
    else if (preselect=="minimum") {res[,relevant:=ifelse(t_train>min(t_train),0,1)]}
  } else {res[,relevant:=1]}
  res[relevant==0,c("coef_est","stderr_est","t_est"):=NA]
  setcolorder(res, c("leaf","G_train","N_train","coef_train","stderr_train","t_train","relevant","G_est","N_est","coef_est","stderr_est","t_est"))
  return(list(tree=tree,res=res))
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
                         folds = NULL,              # name of fold-id column (within-sample crossfit)
                         forest_opts = list(),
                         hat_suffix = ".hat",
                         verbose = FALSE) {
  
  stopifnot(data.table::is.data.table(DT))
  stopifnot(is.character(y_name), length(y_name) == 1L)
  stopifnot(is.character(x_names), length(x_names) >= 1L)
  stopifnot("sample" %chin% names(DT))
  stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))
  
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
  
  # ---- Precompute X for rows i ONCE ----
  X_all <- as.matrix(DT[i, ..x_names])
  
  # rowid within i subset for fast X indexing and preds_buf indexing
  rid_col <- ".__rid__"
  DT[i, (rid_col) := seq_len(n_i)]
  rid <- DT[[rid_col]]  # local alias for speed
  
  # pull vectors once (and coerce once)
  sample_all <- DT[["sample"]]
  y_all <- as.numeric(DT[[y_name]])
  w_all <- if (!is.null(weight_name)) as.numeric(DT[[weight_name]]) else NULL
  f_all <- if (!is.null(folds)) DT[[folds]] else NULL
  
  # output buffer aligned with i via rid
  preds_buf <- rep(NA_real_, n_i)
  
  # Fit on idx_tr, predict on idx_te
  fit_predict_idx <- function(idx_tr, idx_te) {
    n_te <- length(idx_te)
    preds <- rep(NA_real_, n_te)
    if (n_te == 0L) return(preds)
    
    n_tr <- length(idx_tr)
    if (n_tr == 0L) return(preds)
    
    y_tr <- y_all[idx_tr]
    w_tr <- if (is.null(w_all)) NULL else w_all[idx_tr]
    
    # fallback: too few training rows
    if (n_tr < 2L) {
      mu <- if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
      preds[] <- mu
      return(preds)
    }
    
    # fallback: constant outcome (ignore non-finite) without unique()
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
        sample.weights = w_tr
      ), forest_opts)
    )
    
    X_te <- X_all[rid[idx_te], , drop = FALSE]
    preds[] <- as.numeric(predict(rf, X_te)$predictions)
    preds
  }
  
  # Within-sample K-fold predictions for one (sample, margins-cell) subset
  within_pred_idx <- function(idx_s) {
    n <- length(idx_s)
    if (n == 0L) return(numeric())
    
    # if no folds, fit once and predict in-sample
    if (is.null(f_all)) {
      fit <- do.call(
        grf::regression_forest,
        c(list(
          X = X_all[rid[idx_s], , drop = FALSE],
          Y = y_all[idx_s],
          sample.weights = if (is.null(w_all)) NULL else w_all[idx_s]
        ), forest_opts)
      )
      return(as.numeric(predict(fit, X_all[rid[idx_s], , drop = FALSE])$predictions))
    }
    
    f <- f_all[idx_s]
    if (data.table::uniqueN(f) < 2L) {
      fit <- do.call(
        grf::regression_forest,
        c(list(
          X = X_all[rid[idx_s], , drop = FALSE],
          Y = y_all[idx_s],
          sample.weights = if (is.null(w_all)) NULL else w_all[idx_s]
        ), forest_opts)
      )
      return(as.numeric(predict(fit, X_all[rid[idx_s], , drop = FALSE])$predictions))
    }
    
    # precompute positions per fold once; avoid repeated which(f==k)
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
      
      # fallback if training too small
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
  
  # ---- Build group index lists once (kept as-is; straightforward) ----
  idx_all <- i
  if (length(grp_cols) == 0L) {
    group_list <- list(idx_all)
  } else {
    group_dt <- DT[i, ..grp_cols]
    gid <- data.table::frankv(group_dt, ties.method = "dense")
    group_list <- split(idx_all, gid)
  }
  
  # ---- Main loop over groups ----
  for (g in seq_along(group_list)) {
    idx_g <- group_list[[g]]
    if (length(idx_g) == 0L) next
    
    idx1 <- idx_g[sample_all[idx_g] == 1L]
    idx2 <- idx_g[sample_all[idx_g] == 2L]
    
    if (is.null(folds)) {
      # ACROSS-SAMPLE mode: train 1 -> predict 2 and train 2 -> predict 1
      if (length(idx1) > 1L && length(idx2) > 0L) {
        preds_buf[rid[idx2]] <- fit_predict_idx(idx1, idx2)
      } else if (length(idx1) == 1L && length(idx2) > 0L) {
        # 1 training obs: mean fallback
        mu <- y_all[idx1]
        preds_buf[rid[idx2]] <- mu
      }
      
      if (length(idx2) > 1L && length(idx1) > 0L) {
        preds_buf[rid[idx1]] <- fit_predict_idx(idx2, idx1)
      } else if (length(idx2) == 1L && length(idx1) > 0L) {
        mu <- y_all[idx2]
        preds_buf[rid[idx1]] <- mu
      }
      
    } else {
      # WITHIN-SAMPLE crossfit mode: out-of-fold predictions within each sample part
      if (length(idx1) > 0L) preds_buf[rid[idx1]] <- within_pred_idx(idx1)
      if (length(idx2) > 0L) preds_buf[rid[idx2]] <- within_pred_idx(idx2)
    }
    
    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("crossfit_hat: processed %d/%d groups", g, length(group_list)))
    }
  }
  
  # writeback + cleanup
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
                       folds = NULL,              # name of fold-id column (within-sample crossfit for pred only)
                       weight_name = NULL,
                       cluster_name = NULL,
                       forest_opts = NULL,
                       aipw_clip = 1e-3,
                       verbose = FALSE) {
  
  stopifnot(data.table::is.data.table(DT))
  forest_type <- match.arg(forest_type)
  
  stopifnot("sample" %chin% names(DT))
  stopifnot(all(DT[["sample"]] %in% c(1L, 2L)))
  
  if (!is.null(folds)) stopifnot(is.character(folds), length(folds) == 1L, folds %chin% names(DT))
  if (!is.null(aipw_clip)) {
    stopifnot(is.numeric(aipw_clip), length(aipw_clip) == 1L, is.finite(aipw_clip), aipw_clip > 0, aipw_clip < 1)
  }
  
  if (is.null(i)) i <- DT[, .I]
  i <- as.integer(i)
  n_i <- length(i)
  if (n_i == 0L) return(invisible(DT))
  
  if (is.null(margins)) margins <- character()
  grp_cols <- as.character(margins)
  
  if (is.null(forest_opts)) forest_opts <- list()
  
  # Output cols
  if (!("pred"   %chin% names(DT))) DT[, pred   := NA_real_]
  if (!("pred_o" %chin% names(DT))) DT[, pred_o := NA_real_]
  if (!("scores" %chin% names(DT))) DT[, scores := NA_real_]
  
  # Nuisance names (must exist for causal/IV)
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
  
  # ---- Precompute X for rows i ONCE ----
  X_all <- as.matrix(DT[i, ..x_names])
  
  # rowid within the i-subset (used for both X indexing + output buffers)
  rid_col <- ".__rid__"
  DT[i, (rid_col) := seq_len(n_i)]
  rid <- DT[[rid_col]]  # local alias
  
  # Pull commonly used vectors once (coerce once where useful)
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
  
  # ---- Forest builder taking DT-row indices ----
  build_forest_idx <- function(type, idx) {
    X <- X_all[rid[idx], , drop = FALSE]
    Y <- y_all[idx]
    sw <- if (is.null(wgt_all)) NULL else wgt_all[idx]
    cl <- if (is.null(cl_all)) NULL else cl_all[idx]
    
    if (type == "regression") {
      do.call(grf::regression_forest,
              c(list(X = X, Y = Y, sample.weights = sw, clusters = cl), forest_opts))
    } else if (type == "causal") {
      W <- w_all[idx]
      do.call(grf::causal_forest,
              c(list(X = X, Y = Y, W = W,
                     Y.hat = yhat_all[idx],
                     W.hat = what_all[idx],
                     sample.weights = sw,
                     clusters = cl), forest_opts))
    } else { # instrumental
      W <- w_all[idx]
      Z <- z_all[idx]
      do.call(grf::instrumental_forest,
              c(list(X = X, Y = Y, W = W, Z = Z,
                     Y.hat = yhat_all[idx],
                     W.hat = what_all[idx],
                     Z.hat = zhat_all[idx],
                     sample.weights = sw,
                     clusters = cl), forest_opts))
    }
  }
  
  # ---- Within-sample predictions (pred) ----
  within_pred_idx <- function(idx_s) {
    n <- length(idx_s)
    if (n == 0L) return(numeric())
    
    Xs <- X_all[rid[idx_s], , drop = FALSE]
    
    # no folds => single fit, in-sample preds
    if (is.null(folds_all)) {
      fit <- build_forest_idx(forest_type, idx_s)
      return(as.numeric(predict(fit, Xs)$predictions))
    }
    
    f <- folds_all[idx_s]
    if (data.table::uniqueN(f) < 2L) {
      fit <- build_forest_idx(forest_type, idx_s)
      return(as.numeric(predict(fit, Xs)$predictions))
    }
    
    # precompute positions per fold once; avoid repeated which(f==k)
    pos_by_fold <- split(seq_len(n), f)
    
    p <- rep(NA_real_, n)
    tr_mask <- rep.int(TRUE, n)
    
    for (k in names(pos_by_fold)) {
      te_pos <- pos_by_fold[[k]]
      tr_mask[te_pos] <- FALSE
      tr_pos <- which(tr_mask)
      tr_mask[te_pos] <- TRUE
      
      # if training set too small, fallback to mean Y (weighted if available)
      if (length(tr_pos) < 2L) {
        idx_tr <- idx_s[tr_pos]
        y_tr <- y_all[idx_tr]
        w_tr <- if (is.null(wgt_all)) NULL else wgt_all[idx_tr]
        mu <- if (length(y_tr) == 0L) NA_real_ else {
          if (is.null(w_tr)) mean(y_tr) else stats::weighted.mean(y_tr, w_tr)
        }
        p[te_pos] <- mu
      } else {
        fit_k <- build_forest_idx(forest_type, idx_s[tr_pos])
        p[te_pos] <- as.numeric(predict(fit_k, Xs[te_pos, , drop = FALSE])$predictions)
      }
    }
    p
  }
  
  # ---- AIPW scores using OUT-OF-SAMPLE tau (pred_o) ----
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
    
    # if e not in [0,1], assume logit index and convert
    if (any(e < 0 | e > 1, na.rm = TRUE)) e <- plogis(e)
    if (!is.null(aipw_clip)) e <- pmin(pmax(e, aipw_clip), 1 - aipw_clip)
    
    tau <- as.numeric(tau_vec)
    stopifnot(length(tau) == n)
    
    m1 <- m + (1 - e) * tau
    m0 <- m - e * tau
    
    tau + (W / e) * (Y - m1) - ((1 - W) / (1 - e)) * (Y - m0)
  }
  
  # ---- Build group index lists once (kept for simplicity) ----
  idx_all <- i
  if (length(grp_cols) == 0L) {
    group_list <- list(idx_all)
  } else {
    group_dt <- DT[i, ..grp_cols]
    gid <- data.table::frankv(group_dt, ties.method = "dense")
    group_list <- split(idx_all, gid)
  }
  
  # ---- Output buffers aligned with rid ----
  pred_buf   <- rep(NA_real_, n_i)
  pred_o_buf <- rep(NA_real_, n_i)
  score_buf  <- rep(NA_real_, n_i)
  
  # ---- Main loop ----
  for (g in seq_along(group_list)) {
    idx_g <- group_list[[g]]
    if (length(idx_g) == 0L) next
    
    idx1 <- idx_g[sample_all[idx_g] == 1L]
    idx2 <- idx_g[sample_all[idx_g] == 2L]
    
    # within-sample (possibly crossfit) predictions
    p1 <- if (length(idx1)) within_pred_idx(idx1) else numeric()
    p2 <- if (length(idx2)) within_pred_idx(idx2) else numeric()
    
    # fit full forests in each sample for out-of-sample prediction (pred_o)
    fit1 <- if (length(idx1)) build_forest_idx(forest_type, idx1) else NULL
    fit2 <- if (length(idx2)) build_forest_idx(forest_type, idx2) else NULL
    
    p1_o <- if (!is.null(fit2) && length(idx1)) {
      as.numeric(predict(fit2, X_all[rid[idx1], , drop = FALSE])$predictions)
    } else numeric()
    
    p2_o <- if (!is.null(fit1) && length(idx2)) {
      as.numeric(predict(fit1, X_all[rid[idx2], , drop = FALSE])$predictions)
    } else numeric()
    
    sc1 <- if (length(idx1)) compute_aipw_vec(idx1, p1_o) else numeric()
    sc2 <- if (length(idx2)) compute_aipw_vec(idx2, p2_o) else numeric()
    
    # write into buffers using rid (no match())
    if (length(idx1)) {
      r1 <- rid[idx1]
      pred_buf[r1]   <- p1
      pred_o_buf[r1] <- p1_o
      score_buf[r1]  <- sc1
    }
    if (length(idx2)) {
      r2 <- rid[idx2]
      pred_buf[r2]   <- p2
      pred_o_buf[r2] <- p2_o
      score_buf[r2]  <- sc2
    }
    
    if (verbose && (g %% 25L == 0L)) {
      message(sprintf("fit_models: processed %d/%d groups", g, length(group_list)))
    }
  }
  
  # ---- Single writeback ----
  DT[i, `:=`(pred = pred_buf, pred_o = pred_o_buf, scores = score_buf)]
  
  # cleanup
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
  
  allowed_pool <- unique(c(margins, sample_col))
  if (is.null(pool)) pool <- allowed_pool
  pool <- unique(as.character(pool))
  bad_pool <- setdiff(pool, allowed_pool)
  if (length(bad_pool) > 0L) stop("pool contains invalid names: ", paste(bad_pool, collapse = ", "))
  
  pool_sample <- sample_col %in% pool
  pool_margins <- intersect(pool, margins)
  
  cell_cols <- setdiff(margins, pool_margins)
  test_by <- c(cell_cols, if (!pool_sample) sample_col else character())
  
  share_denom_by <- c(cell_cols, if (!pool_sample) sample_col else character())
  do_shares <- (length(pool_margins) > 0L)
  
  search_grid_on <- !is.null(gridpoints)
  if (search_grid_on) {
    stopifnot(is.numeric(gridpoints), length(gridpoints) == 1L, is.finite(gridpoints), gridpoints >= 2)
    gridpoints <- as.integer(gridpoints)
  }
  
  # ---- base vectors once (big speed win) ----
  n <- nrow(data)
  samp <- as.integer(data[[sample_col]])
  stopifnot(all(samp %in% c(1L, 2L)))
  
  predv   <- as.numeric(data[[pred_col]])
  predov  <- as.numeric(data[[pred_o_col]])
  scorev  <- as.numeric(data[[scores_col]])
  
  wv <- if (!is.null(weight_col)) as.numeric(data[[weight_col]]) else rep(1.0, n)
  wv[!is.finite(wv)] <- 0
  
  clv <- if (!is.null(cluster_col)) {
    # precompute integer cluster ids ONCE
    as.integer(factor(data[[cluster_col]], exclude = NULL))
  } else {
    # iid fallback
    seq_len(n)
  }
  
  # helpers ------------------------------------------------------------
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
  
  # weighted mean (robust to NA / non-finite)
  w_mean <- function(x, w) {
    ok <- is.finite(x) & is.finite(w)
    if (!any(ok)) return(NA_real_)
    sum(w[ok] * x[ok]) / sum(w[ok])
  }
  
  # weighted SD (population version, matches previous code)
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
  
  # ---- build cell index lists once (no copying of data[idx]) ----
  if (length(cell_cols) == 0L) {
    idx_list <- list(seq_len(n))
    keys_list <- list(NULL)
  } else {
    gid <- data.table::frankv(data[, ..cell_cols], ties.method = "dense")
    idx_list <- split(seq_len(n), gid)
    keys_list <- lapply(idx_list, function(idx) data[idx[1L], ..cell_cols])
  }
  
  # ---- per-cell runner without copying df_cell ----
  run_one_cell_idx <- function(idx, key_dt = NULL) {
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
      rid    = seq_along(idx)   # mapping back to idx
    )
    
    data.table::setorder(dt, sample, pred)
    
    # recursion
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
    
    # dual constraint
    tau_tr <- dt[G >= minsize, .(tau_tr = min(pred, na.rm = TRUE)), by = sample]
    
    dt_est <- dt[, .(sample, pred_o, cl, w)]
    data.table::setorder(dt_est, sample, pred_o)
    dt_est[, G_est := cumsum(!duplicated(cl)), by = sample]
    tau_est <- dt_est[G_est >= minsize, .(tau_est = min(pred_o, na.rm = TRUE)), by = sample]
    
    if (nrow(tau_tr) < 2L || nrow(tau_est) < 2L) {
      msg <- "minsize too large: cannot reach minsize clusters in train or est order for both samples."
      if (!is.null(key_dt) && ncol(key_dt) > 0L) {
        msg <- paste0(msg, " Cell: ",
                      paste(paste(names(key_dt), as.character(key_dt[1, ]), sep="="), collapse=", "))
      }
      stop(msg)
    }
    
    tau_vec <- c(
      max(as.numeric(c(tau_tr[sample == 1L, tau_tr], tau_est[sample == 2L, tau_est]))),
      max(as.numeric(c(tau_tr[sample == 2L, tau_tr], tau_est[sample == 1L, tau_est])))
    )
    dt[, tau_constraint := data.table::fifelse(sample == 1L, tau_vec[1], tau_vec[2])]
    
    # search candidates
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
              .SDcols = c("G","N","m","se","t_stat","pred")]
    
    res[, train := TRUE]
    
    cut1 <- res[sample == 1L, pred]
    cut2 <- res[sample == 2L, pred]
    if (length(cut1) == 0L || !is.finite(cut1)) cut1 <- tau_vec[1]
    if (length(cut2) == 0L || !is.finite(cut2)) cut2 <- tau_vec[2]
    
    # stored grid (100 + chosen)
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
    
    # apply swapped cutoffs and map indices correctly using rid (BUG FIX)
    dt[, tau_test := data.table::fifelse(sample == 1L, cut2, cut1)]
    in_test <- (dt$pred_o <= dt$tau_test)
    idx_test <- idx[ dt$rid[in_test] ]
    
    # training output
    train_out <- data.table::data.table(
      train = TRUE,
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
      data.table::setcolorder(train_out, c(names(key_dt), "train", "sample"))
    }
    
    list(train = train_out, grid = grid_out, idx_test = idx_test)
  }
  
  # ---- run all cells ----
  n_cells <- length(idx_list)
  train_list <- vector("list", n_cells)
  grid_list  <- if (store_grid) vector("list", n_cells) else NULL
  idx_test_list <- vector("list", n_cells)
  
  for (g in seq_len(n_cells)) {
    ans <- run_one_cell_idx(idx_list[[g]], keys_list[[g]])
    train_list[[g]] <- ans$train
    idx_test_list[[g]] <- ans$idx_test
    if (store_grid) grid_list[[g]] <- ans$grid
    if (verbose && (g %% 25L == 0L)) message(sprintf("forest_test: processed %d/%d cells", g, n_cells))
  }
  
  train_out <- data.table::rbindlist(train_list, use.names = TRUE, fill = TRUE)
  grid_out  <- if (store_grid) data.table::rbindlist(grid_list, use.names = TRUE, fill = TRUE) else NULL
  
  idx_test_all <- unlist(idx_test_list, use.names = FALSE)
  in_test <- rep(FALSE, n)
  if (length(idx_test_all) > 0L) in_test[idx_test_all] <- TRUE
  if (!any(in_test)) stop("Testing subset is empty.")
  
  # ---- testing data (built once) ----
  dt_test <- data.table::data.table(
    score = scorev,
    w     = wv,
    cl    = clv
  )
  if (length(test_by) > 0L) dt_test[, (test_by) := data[, ..test_by]]
  dt_test <- dt_test[in_test]
  
  # testing output
  if (length(test_by) == 0L) {
    o <- crv1_mean(dt_test$score, dt_test$w, dt_test$cl)
    test_out <- data.table::data.table(
      train = FALSE,
      sample = NA_integer_,
      G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
      tau_cutoff = NA_real_,
      p.raw = stats::pnorm(o$t)
    )
  } else {
    test_out <- dt_test[, {
      o <- crv1_mean(score, w, cl)
      data.table::data.table(
        train = FALSE,
        G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
        tau_cutoff = NA_real_,
        p.raw = stats::pnorm(o$t)
      )
    }, by = test_by]
    if (!(sample_col %in% names(test_out))) test_out[, (sample_col) := NA_integer_]
  }
  
  results_out <- data.table::rbindlist(list(train_out, test_out), use.names = TRUE, fill = TRUE)
  
  # ---- global means by same strata as testing ----
  dt_all <- data.table::data.table(score = scorev, w = wv, cl = clv)
  if (length(test_by) > 0L) dt_all[, (test_by) := data[, ..test_by]]
  
  if (length(test_by) == 0L) {
    o <- crv1_mean(dt_all$score, dt_all$w, dt_all$cl)
    global_dt <- data.table::data.table(
      train = FALSE,
      sample = NA_integer_,
      G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
      tau_cutoff = NA_real_,
      p.raw = stats::pnorm(o$t)
    )
  } else {
    global_dt <- dt_all[, {
      o <- crv1_mean(score, w, cl)
      data.table::data.table(
        train = FALSE,
        G = o$G, N = o$N, coef = o$coef, stderr = o$se, t = o$t,
        tau_cutoff = NA_real_,
        p.raw = stats::pnorm(o$t)
      )
    }, by = test_by]
    if (!(sample_col %in% names(global_dt))) global_dt[, (sample_col) := NA_integer_]
  }
  
  # ---- X summaries ----
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
  
  # ---- shares: long format by pooled margins within denom groups (cell_cols +/- sample) ----
  shares <- NULL
  if (do_shares) {
    denom_by <- c(cell_cols, if (!pool_sample) sample_col else character())
    cols_need <- unique(c(denom_by, pool_margins))
    DT_all  <- data[, ..cols_need]
    DT_tst  <- data[in_test, ..cols_need]
    
    make_share_long <- function(DTsub, wsub, vars) {
      DTsub <- data.table::as.data.table(DTsub)
      DTsub[, w := wsub]
      out_list <- vector("list", length(vars))
      for (j in seq_along(vars)) {
        v <- vars[j]
        by_num <- unique(c(denom_by, v))
        num <- DTsub[, .(w_sum = sum(w)), by = by_num]
        if (length(denom_by) == 0L) {
          num[, den_w := sum(DTsub$w)]
        } else {
          den <- DTsub[, .(den_w = sum(w)), by = denom_by]
          num <- num[den, on = denom_by]
        }
        num[, share := ifelse(is.finite(den_w) & den_w > 0, w_sum / den_w, NA_real_)]
        num[, c("w_sum", "den_w") := NULL]
        num[, var := v]
        data.table::setnames(num, v, "level")
        out_list[[j]] <- num
      }
      data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
    }
    
    sh     <- make_share_long(DT_tst, wv[in_test], pool_margins)
    sh_all <- make_share_long(DT_all, wv,          pool_margins)
    
    data.table::setnames(sh, "share", "share")
    data.table::setnames(sh_all, "share", "share_all")
    
    by_merge <- c(denom_by, "var", "level")
    shares <- merge(sh, sh_all, by = by_merge, all = TRUE)
    data.table::setcolorder(shares, c(denom_by, "var", "level", "share", "share_all"))
  }
  
  out <- list(
    results = results_out,
    grid    = grid_out,
    global  = global_dt,
    Xmeans  = Xmeans,
    Xmeans_all = Xmeans_all,
    XSD = XSD
  )
  if (!is.null(shares)) out$shares <- shares
  if (!store_grid) out$grid <- NULL
  
  out
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

