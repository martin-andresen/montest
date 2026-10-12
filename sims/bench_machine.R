## Machine benchmark: small parallel montest Monte Carlo (MW + KR, null DGP, 1-thread forests per worker).
## Usage (from the package root): Rscript sims/bench_machine.R [CORES] [NREP] [PKGDIR]
##   CORES default 8 (capped at 16), NREP default 32 (x2 methods), PKGDIR default "." (needs devtools or installed montest)
## Compare the printed wall time across machines (use the same CORES/NREP); also reports serial time of one job.
args  <- commandArgs(trailingOnly = TRUE)
CORES <- min(if (length(args) >= 1) as.integer(args[1]) else 8L, 16L)
NREP  <- if (length(args) >= 2) as.integer(args[2]) else 32L
PKG   <- normalizePath(if (length(args) >= 3) args[3] else ".")
N <- 1500L; TREES <- 200L
library(parallel)
cat(R.version.string, "|", Sys.info()[["sysname"]], R.version$arch, "| cores:", detectCores(), "| using:", CORES, "| reps:", NREP, "\n")

gen <- function(n) {
  x <- rnorm(n); Z <- rbinom(n, 1, 0.5); V <- rnorm(n)
  D <- as.integer(V + Z > -0.5) + as.integer(V + Z > 0.7)
  lat <- 0.5 * D + 0.5 * x + 0.5 * V + rlogis(n)
  data.frame(Y = as.integer(lat > 0.5), D, Z, x)
}
run_one <- function(i) {
  set.seed(5000 + i); d <- gen(N)
  GRF <- list(num.trees = TREES, num.threads = 1)
  common <- list(fml = Y ~ x | D ~ Z, data = d, pool = "none", progress = FALSE,
                 Zparameters = GRF, Qparameters = GRF, Cparameters = GRF)
  p1 <- montest(fml = common$fml, data = d, condition = "MW", pool = "none", progress = FALSE, Zparameters = GRF, Qparameters = GRF, Cparameters = GRF)$minp
  p2 <- do.call(montest, c(common, list(condition = "KR", Ysubsets = 2L, Dsubsets = 3L)))$minp
  c(MW = unname(p1["p.raw"]), KR = unname(p2["p.raw"]))
}
load_pkg <- function(pkg) suppressMessages(if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(pkg, "DESCRIPTION")))
  devtools::load_all(pkg, quiet = TRUE) else library(montest))

load_pkg(PKG)
t1 <- system.time(run_one(0))[["elapsed"]]
cat(sprintf("Serial: 1 rep (MW+KR)        %6.1f s\n", t1))

cl <- makeCluster(CORES)
clusterExport(cl, c("gen", "run_one", "N", "TREES", "PKG", "load_pkg"))
invisible(clusterEvalQ(cl, load_pkg(PKG)))
clusterSetRNGStream(cl, 1)
t2 <- system.time(res <- parLapply(cl, seq_len(NREP), run_one))[["elapsed"]]
stopCluster(cl)
cat(sprintf("Parallel: %d reps on %d cores  %6.1f s   (speedup vs serial est: %.1fx)\n", NREP, CORES, t2, t1 * NREP / t2))
p <- do.call(rbind, res); cat("Sanity, mean raw p (null, ~0.5ish):", round(colMeans(p, na.rm = TRUE), 2), "\n")
