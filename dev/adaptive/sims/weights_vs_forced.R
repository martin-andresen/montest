## New weighted-Bonferroni adaptive method vs forced pool, forced select (old stepdown screen)
## and the old adaptive default. Usage: Rscript sims/weights_vs_forced.R pkg reps n threads scenario out
args <- commandArgs(TRUE)
pkg <- args[1]; reps <- as.integer(args[2]); n <- as.integer(args[3])
threads <- min(as.integer(args[4]), 16L); sc <- args[5]; out <- args[6]
suppressMessages(pkgload::load_all(pkg, quiet = TRUE))
data.table::setDTthreads(2L)

dgp <- function(sc) {
  x <- rnorm(n); z <- sample(0:4, n, TRUE); u <- runif(n)
  p <- 0.2 + 0.1 * z; bad <- x > 0
  if (sc == "diffuse")      p[bad] <- 0.6 - 0.05 * z[bad]
  if (sc == "concentrated") p[bad] <- c(0.3, 0.4, 0.5, 0.35, 0.45)[z[bad] + 1]
  d <- as.integer(u < p)
  y <- rnorm(n) + 0.5 * d + if (sc == "null") 0 else bad * d
  data.frame(y, x, d, z)
}
fit <- function(dat, method = "legacy", ...) {
  options(montest.adaptive_method = method)
  tryCatch(montest(y ~ x | d ~ z, data = dat, condition = "MW", progress = FALSE,
                   Rparameters = list(num.threads = threads), ...), error = function(e) { message(conditionMessage(e)); NULL })
}
pv <- function(o, nm) if (is.null(o)) NA_real_ else o$minp[[nm]]
res <- list()
for (r in seq_len(reps)) {
  t_r <- proc.time()[3]
  set.seed(1000 + r); dat <- dgp(sc)
  op <- fit(dat, pool = c("sample", "zmargin", "equation"), select = "none")
  os <- fit(dat, pool = "sample", select = c("zmargin", "equation"))
  oa <- fit(dat)
  ow <- fit(dat, method = "weights")
  wsum <- if (is.null(ow)) NA_character_ else {
    aw <- ow$adaptive_weights; aw <- aw[, .(w = sum(w)), by = design_id]
    paste(sprintf("d%d=%.2f", aw$design_id, aw$w), collapse = " ") }
  res[[r]] <- data.frame(sc, r, p_pool = pv(op, "p.holm"), p_select = pv(os, "p.holm"),
                         p_adaptive_old = pv(oa, "p.holm"), p_weights = pv(ow, "p.wbonf"), p_weights_hi = pv(ow, "p.wbonf_hi"), p_equal = pv(ow, "p.equal"), wsum = wsum)
  message(sprintf("rep %d: %.0fs", r, proc.time()[3] - t_r))
  if (r %% 5 == 0) saveRDS(do.call(rbind, res), paste0(out, ".rds"))
}
res <- do.call(rbind, res); saveRDS(res, paste0(out, ".rds"))
print(sapply(res[, 3:8], function(p) c(rate = mean(p < 0.05, na.rm = TRUE), n_na = sum(is.na(p)))))
