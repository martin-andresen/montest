## MW vs KR, trivalued D (0,1,2), binary or trivalued Y. Small and quick.
## Best case for KR: with a binary / 3-point Y the KR sets are few and need no binning.
## Both methods are run with the same pool ("none"), so the pooling question is held out.
##
## Usage: Rscript sims/mw_vs_kr.R [NREP] [N] [scenario substring, e.g. "interior"]
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[1]) else 30L
N     <- if (length(args) >= 2) as.integer(args[2]) else 1500L
ONLY  <- if (length(args) >= 3 && args[3] != "all") args[3] else NULL
MODE  <- if (length(args) >= 4) args[4] else "none"      # "none" | "pool" | "select": pool/select over ALL margins
stopifnot(MODE %in% c("none", "pool", "select"))
CORES <- 8L                                   # <= 16 (project rule)
TREES <- 200L

library(parallel)
cl <- makeCluster(CORES)
clusterEvalQ(cl, suppressMessages(devtools::load_all("C:/Users/martiea/Dropbox/Prosjekter/working/montest", quiet = TRUE)))

## DGP: Z randomised, D = 1{V+Z>-.5}+1{V+Z>.7} (monotone), Y confounded through V.
## `cell` = D value in which Z shifts the latent outcome directly (exclusion violation); NA = valid IV.
## ky = number of Y categories (2 or 3); fs = first-stage strength (a weak first stage shrinks the
## complier "budget" that hides interior violations).
gen <- function(n, ky, cell = NA, gamma = 0, fs = 1) {
  x <- rnorm(n); Z <- rbinom(n, 1, 0.5); V <- rnorm(n)
  D <- as.integer(V + fs * Z > -0.5) + as.integer(V + fs * Z > 0.7)
  lat <- 0.5 * D + 0.5 * x + 0.5 * V + rlogis(n)
  if (!is.na(cell)) lat <- lat + gamma * Z * (D == cell)
  Y <- if (ky == 2) as.integer(lat > 0.5) else as.integer(lat > -0.3) + as.integer(lat > 1.3)
  data.frame(Y, D, Z, x)
}

scen <- data.frame(scenario = c("null", "interior (D=1)", "endpoint (D=2)", "interior, weak 1st stage"),
                   cell = c(NA, 1, 2, 1), gamma = c(0, 3, 3, 8), fs = c(1, 1, 1, 0.5), stringsAsFactors = FALSE)
if (!is.null(ONLY)) scen <- scen[grepl(ONLY, scen$scenario), ]
grid <- expand.grid(rep = seq_len(NREP), s = seq_len(nrow(scen)), ky = c(2, 3))
grid$id <- seq_len(nrow(grid))
clusterExport(cl, c("gen", "scen", "grid", "N", "TREES", "MODE"))

run_one <- function(i) {
  g <- grid[i, ]
  set.seed(5000 + g$id)
  d <- gen(N, g$ky, scen$cell[g$s], scen$gamma[g$s], scen$fs[g$s])
  GRF <- list(num.trees = TREES, num.threads = 1)
  ps <- switch(MODE, none = list(pool = "none"),
                     pool = list(pool = "all"),
                     select = list(pool = "none", select = "all"))
  common <- c(list(data = d, progress = FALSE, Zparameters = GRF, Qparameters = GRF, Cparameters = GRF), ps)
  one <- function(cond, extra = list()) {
    t0 <- proc.time()[["elapsed"]]
    p <- tryCatch(do.call(montest, c(list(fml = Y ~ x | D ~ Z, condition = cond), common, extra))$minp,
                  error = function(e) c(p.raw = NA, p.holm = NA))
    data.frame(rep = g$rep, scenario = scen$scenario[g$s], ky = g$ky, method = cond,
               p.raw = unname(p["p.raw"]), p.holm = unname(p["p.holm"]), secs = proc.time()[["elapsed"]] - t0)
  }
  rbind(one("MW"), one("KR", list(Ysubsets = g$ky, Dsubsets = 3L)))
}
clusterExport(cl, "run_one")
res <- parLapply(cl, seq_len(nrow(grid)), run_one)
stopCluster(cl)

r <- do.call(rbind, res)
out <- aggregate(cbind(rej_holm05 = p.holm < .05, rej_raw05 = p.raw < .05, secs = secs, failed = is.na(p.holm)) ~ ky + scenario + method,
                 data = r, FUN = mean, na.action = na.pass)
out <- out[order(out$ky, out$scenario, out$method), ]
cat("\nN =", N, " reps =", NREP, " trees =", TREES, " mode =", MODE,
    c(none = "(pool = 'none')", pool = "(pool = 'all')", select = "(pool = 'none', select = 'all')")[[MODE]],
    "; secs = mean seconds per test, 1 thread\n")
print(out, row.names = FALSE, digits = 3)
saveRDS(r, file.path(tempdir(), "mw_vs_kr.rds"))
