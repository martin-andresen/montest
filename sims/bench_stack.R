## Timing benchmark: stack = TRUE vs stack = FALSE (same pool/select, so results should agree up to RNG).
## Run sequentially (clean timings); forests use THREADS threads (<= 16, project rule).
##
## Usage: Rscript sims/bench_stack.R [REPS] [Ns, comma separated] [scenario substring]
args  <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 2L
NS    <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(1000L, 4000L)
ONLY  <- if (length(args) >= 3 && args[3] != "all") args[3] else NULL
THREADS <- 8L
TREES   <- 200L

suppressMessages(devtools::load_all("C:/Users/martiea/Dropbox/Prosjekter/working/montest", quiet = TRUE))
library(data.table)

## D in 0..(kd-1) monotone in Z, Y with ky categories, valid IV (null) -- timing only.
gen <- function(n, kd, ky, kz = 2) {
  x <- rnorm(n); Z <- sample(0:(kz - 1), n, TRUE); V <- rnorm(n)
  cuts <- seq(-0.5, 0.7, length.out = kd - 1)
  D <- rowSums(outer(V + 0.8 * Z / (kz - 1), cuts, ">"))
  lat <- 0.5 * D + 0.5 * x + 0.5 * V + rlogis(n)
  Y <- if (ky == 0) lat else findInterval(lat, quantile(lat, seq_len(ky - 1) / ky))
  data.frame(Y, D, Z, x)
}

## label, condition, kd, ky, kz, extra montest args
scen <- list(
  list("MW kd=2",               "MW", 2, 0, 2, list()),
  list("MW kd=3",               "MW", 3, 0, 2, list()),
  list("MW kd=5",               "MW", 5, 0, 2, list()),
  list("KR kd=2 ky=2",          "KR", 2, 2, 2, list(Ysubsets = 2L, Dsubsets = 2L)),
  list("KR kd=3 ky=3",          "KR", 3, 3, 2, list(Ysubsets = 3L, Dsubsets = 3L)),
  list("KR kd=3 ky=6 joint=F",  "KR", 3, 6, 2, list(Ysubsets = 6L, Dsubsets = 3L, joint = FALSE)),
  list("KR kd=5 ky=6 joint=F",  "KR", 5, 6, 2, list(Ysubsets = 6L, Dsubsets = 5L, joint = FALSE)),
  list("KR kd=3 ky=3 kz=3",     "KR", 3, 3, 3, list(Ysubsets = 3L, Dsubsets = 3L)),
  list("simple kd=2",           "simple", 2, 0, 2, list())
)
if (!is.null(ONLY)) scen <- Filter(function(s) grepl(ONLY, s[[1]]), scen)

GRF <- list(num.trees = TREES, num.threads = THREADS)
out <- list()
for (n in NS) for (s in scen) for (r in seq_len(REPS)) {
  set.seed(1000 + r)
  d <- gen(n, s[[3]], s[[4]], s[[5]])
  fml <- if (s[[2]] == "simple") x ~ x | D ~ Z else Y ~ x | D ~ Z
  if (s[[2]] == "simple") fml <- Y ~ x | D ~ Z
  for (st in c(TRUE, FALSE)) {
    base <- list(fml = fml, data = d, condition = s[[2]], pool = "sample",
                 select = c("yval", "dval", "equation", "condition"),
                 progress = FALSE, stack = st, seed = 1,
                 Zparameters = GRF, Qparameters = GRF, Cparameters = GRF)
    t0 <- proc.time()[["elapsed"]]
    m <- tryCatch(do.call(montest, c(base, s[[6]])), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    ok <- !inherits(m, "error")
    out[[length(out) + 1L]] <- data.table(
      scenario = s[[1]], n = n, rep = r, stack = st, secs = el,
      cells = if (ok) nrow(m$margins) else NA_integer_,
      nhyp = if (ok) nrow(m$results[train == FALSE]) else NA_integer_,
      p.raw = if (ok) unname(m$minp["p.raw"]) else NA_real_,
      err = if (ok) NA_character_ else conditionMessage(m))
    cat(sprintf("%-22s n=%5d rep=%d stack=%-5s %7.1fs  cells=%s p=%s %s\n", s[[1]], n, r, st, el,
                out[[length(out)]]$cells, format(out[[length(out)]]$p.raw, digits = 3),
                if (ok) "" else paste("ERR:", conditionMessage(m))))
    flush.console()
  }
}
res <- rbindlist(out)
saveRDS(res, "sims/bench_stack.rds")
print(dcast(res[, .(secs = mean(secs), cells = cells[1], p = mean(p.raw)), by = .(scenario, n, stack)],
            scenario + n + cells ~ stack, value.var = c("secs", "p")))
