## montest() and seqtest() over many seeds, combined into one p-value per correction method.
## Keep this file free of `:=` and of `[` with two index slots: data-raw/reserved_names.R scans every function of
## the package, and new hits would change which user columns montest() aliases.

#' Run montest() over many seeds and combine the p-values
#'
#' The results of \code{\link{montest}} depend on \code{seed} through the sample split and the forests.
#' \code{montest_seeds()} runs the same call once for every seed in a set fixed in advance, and combines the
#' overall minimum p-values (\code{$minp}) with a rule that is also fixed in advance.
#'
#' \strong{Combined p-value.} For each correction method, \code{p_agg = min(1, 2 * median(p))}, where \code{p}
#' holds that method's \code{minp} for all seeds. This is valid under arbitrary dependence between the seeds, for
#' any number of seeds (Rüger 1978; Meinshausen, Meier and Bühlmann 2009), provided each single-seed p-value is
#' valid. Chernozhukov, Demirer, Duflo and Fernández-Val (2025) likewise take medians of p-values over many data
#' splits and adjust the level. \code{p.holm} is the headline. \code{p.raw} is not corrected for multiple testing when a run has more
#' than one test, and each of the other methods is only as valid as its single-seed p-value is. A seed that fails,
#' or has a missing p-value, counts as p = 1. \code{reject_share} is the share of seeds with a p-value of at most
#' \code{alpha}; it describes the variation between seeds and is not a test. Cells are not matched across seeds,
#' because each seed selects its own subgroups, so nothing is combined cell by cell.
#'
#' \strong{Seeds.} Choose the set of seeds before seeing any result; adding or removing seeds afterwards voids the
#' guarantee. A \code{seed} inside one of the \code{*parameters} lists holds that forest's randomness fixed across
#' runs. The state of the caller's random number generator is restored on exit.
#'
#' \strong{Storage and resume.} With \code{file}, the runs are saved to that rds file after every seed, starting
#' with the first successful one. Calling again with the same arguments skips the successful seeds, retries the
#' failed ones and returns the full object. The call stops, without touching the file, if the formula or the
#' arguments in \code{...} (apart from \code{progress} and \code{memcheck}), the data (all columns, in row order),
#' the seeds, the code of the montest package, the version of R, grf, fixest, data.table, rpart or mvtnorm, the
#' random number generator or a \code{montest.*} option differ from those that created the file. Variables or
#' functions that a formula takes from its environment are not checked, nor is the linear algebra library. Only one
#' process may write to a file at a time.
#'
#' \strong{Computation.} The seeds are run one at a time. The forest predictions use all cores whatever
#' \code{num.threads} is set to.
#'
#' @inheritParams montest
#' @param ... Further arguments to \code{\link{montest}}, passed unchanged. All must be named in full, and
#'   \code{seed} may not be among them.
#' @param seeds Integer vector of at least two distinct seeds. Default \code{1:51}, which does not include
#'   \code{montest()}'s default seed.
#' @param file Optional path of an rds file in which to save the runs, see "Storage and resume".
#' @param keep_fits Logical: also keep the full output of every seed (default \code{FALSE}). So that the data are
#'   not saved along with it, the environments of its formulas are replaced by the global environment and the calls
#'   stored in it are replaced or removed. With \code{file}, the whole file is rewritten after every seed.
#'
#' @return An object of class \code{"montest_seeds"}: a list with
#'   \item{p_agg}{Combined p-value of every method.}
#'   \item{reject_share}{Share of seeds rejecting at \code{alpha}, by method.}
#'   \item{summary}{One row per seed: \code{seed}, \code{ok}, the six minimum p-values, the number of tests,
#'     \code{N}, \code{G}, the number of warnings, the time in seconds and the error message.}
#'   \item{results}{The test-sample rows of \code{$results} of every successful seed, with a column \code{seed}.}
#'   \item{log}{Warnings and errors of all runs: \code{seed}, \code{type}, \code{text}.}
#'   \item{fits}{With \code{keep_fits}, the \code{montest()} output by seed, otherwise \code{NULL}.}
#'   \item{seeds, alpha, n_failed, call}{The seeds, the level, the number of failed seeds and the call.}
#'
#' @references Rüger, B. (1978). Das maximale Signifikanzniveau des Tests: "Lehne H0 ab, wenn k unter n gegebenen
#'   Tests zur Ablehnung führen". \emph{Metrika} 25, 171-178.
#'
#'   Meinshausen, N., Meier, L. and Bühlmann, P. (2009). p-Values for high-dimensional regression. \emph{Journal of
#'   the American Statistical Association} 104, 1671-1681.
#'
#'   Chernozhukov, V., Demirer, M., Duflo, E. and Fernández-Val, I. (2025). Generic machine learning inference on
#'   heterogeneous treatment effects in randomized experiments, with an application to immunization in India.
#'   \emph{Econometrica} 93(4), 1121-1164.
#'
#' @seealso montest, seqtest_seeds
#' @examples
#' \dontrun{
#' r <- montest_seeds(Y ~ X | D ~ Z, data = df, condition = "KR", file = "KR_seeds.rds")
#' r$p_agg["p.holm"]
#' }
#' @export
montest_seeds <- function(fml, data, ..., seeds = 1:51, file = NULL, keep_fits = FALSE) {
  call <- mts_clean_call(match.call(), quote(montest_seeds))
  ## Called with `...` so that montest() sees exactly the arguments the user supplied (it uses missing()).
  mts_run("montest", function(s) montest(fml, data, ..., seed = s), fml, data, list(...),
          allowed = names(formals(montest)), seeds = seeds, file = file, keep_fits = keep_fits, call = call)
}

#' Run seqtest() over many seeds and combine the p-values
#'
#' \code{seqtest_seeds()} is to \code{\link{seqtest}} what \code{\link{montest_seeds}} is to
#' \code{\link{montest}}: it runs the same call once for every seed in a set fixed in advance and combines the
#' minimum p-values. The combination rule, the choice of seeds, failed seeds, saving and resuming, and computation
#' are as described in \code{\link{montest_seeds}}; a file is resumed only by the function that created it.
#'
#' \strong{Which p-values.} \code{p_agg} combines \code{seqtest()}'s top-level \code{$minp}, which corrects all
#' cells of all requested conditions as one family; \code{p_agg["p.holm"]} is the headline.
#' \code{p_agg_condition} combines the minimum p-values of each condition by itself. Each of those is valid for a
#' condition chosen in advance, and none is corrected for the other conditions: do not report the smallest one.
#'
#' @inheritParams seqtest
#' @inheritParams montest_seeds
#' @inherit montest_seeds references
#' @param ... Further arguments to \code{\link{seqtest}} (such as \code{condition}) or, through it, to
#'   \code{\link{montest}}, passed unchanged. All must be named in full, and \code{seed} may not be among them.
#'
#' @return An object of class \code{"montest_seeds"} with the elements described in \code{\link{montest_seeds}},
#'   where \code{$results} and \code{$log} have a column \code{seqtest_condition} (\code{condition} is the
#'   underlying \code{montest()} condition; it is \code{NA} in \code{$log} for what cannot be attributed to one
#'   condition), \code{n_tests} counts the tests of all conditions, and \code{N} and \code{G} are \code{NA} if
#'   they differ between conditions. In addition:
#'   \item{p_agg_condition}{Matrix of combined p-values, conditions by correction methods.}
#'   \item{summary_condition}{The minimum p-values of every successful seed and condition.}
#'
#' @seealso seqtest, montest_seeds
#' @examples
#' \dontrun{
#' r <- seqtest_seeds(Y ~ X | D1 + D2 ~ Z, data = df, condition = c("KRD", "FSD"), file = "seq_seeds.rds")
#' r$p_agg["p.holm"]
#' r$p_agg_condition
#' }
#' @export
seqtest_seeds <- function(fml, data, ..., seeds = 1:51, file = NULL, keep_fits = FALSE) {
  call <- mts_clean_call(match.call(), quote(seqtest_seeds))
  mts_run("seqtest", function(s) seqtest(fml, data, ..., seed = s), fml, data, list(...),
          allowed = c(names(formals(seqtest)), names(formals(montest))), seeds = seeds, file = file,
          keep_fits = keep_fits, call = call)
}

## The work of both wrappers. `fun` is "montest" or "seqtest" and `run_seed` runs it for one seed.
mts_run <- function(fun, run_seed, fml, data, dots, allowed, seeds, file, keep_fits, call) {
  wrapper <- paste0(fun, "_seeds")
  if (length(dots) && (is.null(names(dots)) || !all(nzchar(names(dots))))) {
    stop("All arguments in `...` must be named.", call. = FALSE)
  }
  if ("seed" %in% names(dots)) stop("`seed` cannot be passed to ", wrapper, "(): use `seeds`.", call. = FALSE)
  unknown <- setdiff(names(dots), setdiff(allowed, "..."))
  if (length(unknown)) {
    stop("Not an argument of ", if (fun == "seqtest") "seqtest() or ", "montest(): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  seeds <- mts_check_seeds(seeds)
  if (!isTRUE(keep_fits) && !isFALSE(keep_fits)) stop("`keep_fits` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file) || dir.exists(file)) {
      stop("`file` must be NULL or the path of one file.", call. = FALSE)
    }
    if (!dir.exists(dirname(file))) stop("The directory of `file` does not exist: ", dirname(file), call. = FALSE)
  }
  alpha <- if (is.null(dots[["alpha"]])) eval(formals(montest)$alpha) else dots[["alpha"]]
  show <- isTRUE(if ("progress" %in% names(dots)) dots[["progress"]] else interactive())

  ## The caller's random number generator is left as it was.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
  }

  runs <- list()
  if (!is.null(file)) {
    fingerprint <- mts_fingerprint(fun, fml, data, dots, seeds)
    if (file.exists(file)) runs <- mts_load(file, fingerprint)
  }
  done <- vapply(as.character(seeds), function(k) {
    isTRUE(runs[[k]]$summary$ok) && (!keep_fits || !is.null(runs[[k]]$fit))
  }, logical(1L))
  todo <- seeds[!done]
  if (show && any(done)) message(wrapper, ": ", sum(done), " of ", length(seeds), " seeds already in `file`.")

  for (i in seq_along(todo)) {
    s <- todo[i]
    t0 <- proc.time()[[3L]]
    ## Warnings are not signalled further; they are kept in the log of each run.
    warned <- character()
    fit <- tryCatch(
      withCallingHandlers(run_seed(s), warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) e)
    run <- mts_reduce(fun, fit, s, proc.time()[[3L]] - t0, keep_fits, call, warned)
    rm(fit)
    k <- as.character(s)
    ## A seed that is re-run only to obtain its fit keeps its earlier result if the re-run fails.
    if (run$summary$ok || !isTRUE(runs[[k]]$summary$ok)) runs[[k]] <- run
    if (show) {
      message(wrapper, ": seed ", s, " (", i, " of ", length(todo), ") ",
              if (run$summary$ok) "done" else "failed", " after ", round(run$summary$sec, 1), " s")
    }
    ## Nothing is written before a seed has succeeded, so a mistaken call leaves no file behind.
    if (!is.null(file) && any(vapply(runs, function(r) r$summary$ok, logical(1L)))) {
      mts_save(list(format = 1L, fingerprint = fingerprint, runs = runs), file)
    }
  }

  runs <- runs[as.character(seeds)]
  per_seed <- do.call(rbind, lapply(runs, function(r) r$summary))
  rownames(per_seed) <- NULL
  n_failed <- sum(!per_seed$ok)
  if (n_failed) {
    first <- per_seed$error[which(!per_seed$ok)[1L]]
    if (n_failed == length(seeds)) stop(fun, "() failed for every seed. First error: ", first, call. = FALSE)
    warning(n_failed, " of ", length(seeds), " seeds failed and count as p = 1. First error: ", first, call. = FALSE)
  }
  run_log <- do.call(rbind, lapply(runs, function(r) r$log))
  rownames(run_log) <- NULL
  fits <- NULL
  if (keep_fits) fits <- Filter(Negate(is.null), lapply(runs, function(r) r$fit))
  agg <- mts_aggregate(per_seed, alpha)

  out <- list(
    p_agg = agg$p_agg,
    reject_share = agg$reject_share,
    summary = per_seed,
    results = data.table::rbindlist(lapply(runs, function(r) r$results), use.names = TRUE, fill = TRUE),
    log = run_log,
    fits = fits,
    seeds = seeds,
    alpha = alpha,
    n_failed = n_failed,
    call = call
  )
  if (fun == "seqtest") {
    by_condition <- do.call(rbind, lapply(runs, function(r) r$by_condition))
    rownames(by_condition) <- NULL
    out$p_agg_condition <- mts_aggregate_condition(by_condition, seeds)
    out$summary_condition <- by_condition
  }
  structure(out, class = "montest_seeds")
}

#' @export
print.montest_seeds <- function(x, ...) {
  s <- x$summary
  cat(if (is.null(x$p_agg_condition)) "montest" else "seqtest", " over ", length(x$seeds), " seeds (", min(x$seeds), " to ", max(x$seeds), "): ",
      x$n_failed, " failed (counted as p = 1), ", sum(s$n_warnings > 0L), " with warnings\n\n", sep = "")
  p <- mts_pmatrix(s)
  ord <- c("p.holm", setdiff(mts_methods(), "p.holm"))
  labels <- sub("^p[.]", "", ord)
  labels[ord == "p.holm"] <- "holm (headline)"
  labels[ord == "p.raw"] <- "raw (uncorrected)"
  tab <- data.frame(p_agg = signif(unname(x$p_agg[ord]), 4L),
                    median_p = signif(unname(apply(p, 2L, stats::median)[ord]), 4L),
                    reject_share = round(unname(x$reject_share[ord]), 3L),
                    row.names = labels)
  print(tab)
  cat("\np_agg = min(1, 2 * median p over seeds). reject_share: share of seeds with p <= ", x$alpha, ".\n", sep = "")
  if (!is.null(x$p_agg_condition) && nrow(x$p_agg_condition) > 1L) {
    by_condition <- signif(as.data.frame(x$p_agg_condition)[ord], 4L)
    names(by_condition) <- sub("^p[.]", "", ord)
    cat("\np_agg by condition (not corrected across conditions):\n")
    print(by_condition)
    cat("\n")
  }
  n_tests <- range(s$n_tests, na.rm = TRUE)
  if (n_tests[1L] != n_tests[2L]) cat("Tests per seed: ", n_tests[1L], " to ", n_tests[2L], ".\n", sep = "")
  cat("Per seed: $summary. Test rows: $results. Warnings and errors: $log.\n")
  invisible(x)
}

## The correction methods of montest()'s `minp`.
mts_methods <- function() paste0("p.", c("raw", "holm", "hochberg", "BH", "BY", "CCT"))

mts_check_seeds <- function(seeds) {
  ok <- is.numeric(seeds) && length(seeds) >= 2L && all(is.finite(seeds)) &&
    all(seeds == round(seeds)) && all(abs(seeds) <= .Machine$integer.max)
  if (!ok) stop("`seeds` must be a vector of at least two integers without missing values.", call. = FALSE)
  if (anyDuplicated(seeds)) stop("`seeds` must not contain duplicates.", call. = FALSE)
  as.integer(seeds)
}

## Replaces the environment of every formula in `x` (a formula or a list), which would otherwise be saved along
## with it, data and all.
mts_strip_env <- function(x) {
  if (inherits(x, "formula")) {
    environment(x) <- globalenv()
  } else if (is.list(x) && !is.object(x)) {
    x[] <- lapply(x, mts_strip_env)
  }
  x
}

## The wrapper's call without objects that were inlined into it, e.g. by do.call().
mts_clean_call <- function(cl, name) {
  cl[[1L]] <- name
  if (!is.null(cl$data) && !is.language(cl$data)) cl$data <- quote(data)
  for (i in seq_along(cl)[-1L]) {
    if (inherits(cl[[i]], "formula") || is.list(cl[[i]])) cl[[i]] <- mts_strip_env(cl[[i]])
  }
  cl
}

## A montest() or seqtest() object without what would save the data along with it: formula environments and the
## calls, into which seqtest() inlines the data.
mts_slim <- function(x) {
  if (!is.null(x$options)) x$options <- mts_strip_env(x$options)
  x$call <- NULL
  x$montest_call <- NULL
  x
}

mts_log_rows <- function(seed, type, text, condition = NULL) {
  out <- data.frame(seed = rep(seed, length(type)))
  if (!is.null(condition)) out$seqtest_condition <- rep(condition, length(type))
  out$type <- type
  out$text <- text
  out
}

## What is stored of one run. `fit` is the object returned or the error condition of a failed call; `warned` holds
## the messages of the warnings signalled during the call.
mts_reduce <- function(fun, fit, seed, sec, keep_fits, call, warned) {
  failed <- inherits(fit, "error")
  seq_run <- fun == "seqtest"
  ## The montest() objects of the run, named by seqtest() condition.
  parts <- if (failed) {
    list()
  } else if (!seq_run) {
    list(fit)
  } else if (length(fit$condition) == 1L) {
    stats::setNames(list(fit), fit$condition)
  } else {
    unclass(fit)[fit$condition]
  }
  unknown_condition <- if (seq_run) NA_character_

  ## Warnings and errors, from the logs of the montest() objects or of the failed montest() call.
  rows <- function(full, condition) {
    keep <- full$type %in% c("warning", "error")
    mts_log_rows(seed, full$type[keep], full$text[keep], condition)
  }
  run_log <- if (!failed) {
    do.call(rbind, lapply(seq_along(parts), function(i) rows(parts[[i]]$log, if (seq_run) names(parts)[i])))
  } else if (!is.null(fit$log)) {
    rows(fit$log, unknown_condition)
  } else {
    mts_log_rows(seed, "error", conditionMessage(fit), unknown_condition)
  }
  ## Warnings that no montest() log holds: seqtest()'s own, and those before a failure.
  warned <- sub("\n+$", "", warned)
  for (tx in run_log$text[run_log$type == "warning"]) {
    j <- match(tx, warned)
    if (!is.na(j)) warned <- warned[-j]
  }
  run_log <- rbind(run_log, mts_log_rows(seed, rep("warning", length(warned)), warned, unknown_condition))
  rownames(run_log) <- NULL

  tests <- lapply(seq_along(parts), function(i) {
    res <- parts[[i]]$results
    tb <- res[which(res$train == FALSE)]
    data.table::set(tb, j = "seed", value = seed)
    if (seq_run) data.table::set(tb, j = "seqtest_condition", value = names(parts)[i])
    data.table::setcolorder(tb, c("seed", if (seq_run) "seqtest_condition"))
    tb
  })
  ## N and G of the run, if the same for all its montest() objects.
  common <- function(name) {
    v <- unique(vapply(parts, function(x) as.numeric(x$obs[[name]]), numeric(1L)))
    if (length(v) == 1L) v else NA_real_
  }
  p <- if (failed) rep(NA_real_, length(mts_methods())) else unname(fit$minp[mts_methods()])
  row <- data.frame(seed = seed, ok = !failed)
  row[mts_methods()] <- as.list(p)
  row$n_tests <- if (failed) NA_integer_ else sum(vapply(tests, nrow, integer(1L)))
  row$N <- common("N")
  row$G <- common("G")
  row$n_warnings <- sum(run_log$type == "warning")
  row$sec <- unname(sec)
  row$error <- if (failed) conditionMessage(fit) else NA_character_

  by_condition <- NULL
  if (seq_run && !failed) {
    by_condition <- do.call(rbind, lapply(names(parts), function(cn) {
      one <- data.frame(seed = seed, seqtest_condition = cn)
      one[mts_methods()] <- as.list(unname(parts[[cn]]$minp[mts_methods()]))
      one
    }))
  }

  kept <- NULL
  if (keep_fits && !failed) {
    call[[1L]] <- as.name(fun)
    call$seeds <- call$file <- call$keep_fits <- NULL
    call$seed <- seed
    if (seq_run && length(parts) > 1L) {
      kept <- fit
      for (cn in names(parts)) kept[[cn]] <- mts_slim(fit[[cn]])
    } else {
      kept <- mts_slim(fit)
    }
    kept$call <- call
  }
  list(summary = row, results = if (!failed) data.table::rbindlist(tests, use.names = TRUE, fill = TRUE),
       log = run_log, by_condition = by_condition, fit = kept)
}

mts_pmatrix <- function(per_seed) {
  p <- as.matrix(per_seed[mts_methods()])
  p[is.na(p)] <- 1
  p
}

## min(1, 2 * median) and the share of rejections, with missing p-values (failed seeds) counted as 1.
mts_aggregate <- function(per_seed, alpha) {
  p <- mts_pmatrix(per_seed)
  list(p_agg = pmin(2 * apply(p, 2L, stats::median), 1),
       reject_share = colMeans(p <= alpha))
}

## p_agg of every seqtest() condition; seeds without a result count as 1.
mts_aggregate_condition <- function(by_condition, seeds) {
  conditions <- unique(by_condition$seqtest_condition)
  agg <- vapply(conditions, function(cn) {
    hit <- which(by_condition$seqtest_condition == cn)
    vapply(mts_methods(), function(m) {
      p <- rep(1, length(seeds))
      p[match(by_condition$seed[hit], seeds)] <- by_condition[[m]][hit]
      p[is.na(p)] <- 1
      min(2 * stats::median(p), 1)
    }, numeric(1L))
  }, numeric(length(mts_methods())))
  t(agg)
}

mts_md5 <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  saveRDS(x, f, version = 2L, compress = FALSE)
  unname(tools::md5sum(f))
}

## The functions and constants of the package, apart from those of the seed wrappers.
mts_code <- function() {
  ns <- asNamespace("montest")
  nm <- sort(ls(ns, all.names = TRUE), method = "radix")
  nm <- nm[!grepl("^(montest_seeds$|seqtest_seeds$|print[.]montest_seeds$|mts_|[.]__|[.]packageName$)", nm)]
  objects <- mget(nm, envir = ns)
  objects <- objects[!vapply(objects, is.environment, logical(1L))]
  lapply(objects, function(f) {
    if (!is.function(f)) return(f)
    f <- utils::removeSource(f)
    list(deparse(formals(f)), deparse(body(f)))
  })
}

## One checksum for each thing that must be unchanged for saved seeds to be combined with new ones.
mts_fingerprint <- function(fun, fml, data, dots, seeds) {
  args <- c(list(fml = fml), dots[setdiff(names(dots), c("progress", "memcheck"))])
  args <- mts_strip_env(args[sort(names(args), method = "radix")])
  opts <- options()
  packages <- c("grf", "fixest", "data.table", "rpart", "mvtnorm")
  env <- list(R.version.string, R.version$platform,
              vapply(packages, function(pk) as.character(utils::packageVersion(pk)), character(1L)),
              RNGkind(), opts[sort(grep("^montest[.]", names(opts), value = TRUE), method = "radix")])
  list(wrapper = mts_md5(fun),
       args = mts_md5(args),
       data = mts_md5(list(names(data), lapply(seq_along(data), function(j) data[[j]]))),
       seeds = mts_md5(seeds),
       code = mts_md5(mts_code()),
       env = mts_md5(env))
}

## The runs saved in `file`, if it was created by the same call, data, seeds and code.
mts_load <- function(file, fingerprint) {
  saved <- tryCatch(suppressWarnings(readRDS(file)), error = function(e) NULL)
  if (!is.list(saved) || !identical(saved$format, 1L) || !is.list(saved$fingerprint) || !is.list(saved$runs)) {
    stop("`file` exists but was not written by montest_seeds() or seqtest_seeds(), and is not overwritten: ", file,
         call. = FALSE)
  }
  same <- vapply(names(fingerprint), function(n) identical(fingerprint[[n]], saved$fingerprint[[n]]), logical(1L))
  if (!all(same)) {
    what <- c(wrapper = "function (montest_seeds() or seqtest_seeds())", args = "formula or arguments", data = "data", seeds = "seeds", code = "montest code",
              env = "R, package versions, random number generator or montest options")
    stop("The runs in `file` cannot be continued: changed since it was created: ",
         paste(what[names(same)[!same]], collapse = "; "), ". Use a new `file`: ", file, call. = FALSE)
  }
  saved$runs
}

## Written to a temporary file first, so that an interrupted write cannot damage `file`.
mts_save <- function(x, file) {
  tmp <- paste0(file, ".tmp")
  saveRDS(x, tmp)
  if (!file.rename(tmp, file)) stop("Could not write `file`: ", file, call. = FALSE)
}
