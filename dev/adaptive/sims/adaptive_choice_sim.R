## Adaptive pool/select: rejection rates AND which design wins.
## Usage: Rscript sims/adaptive_choice_sim.R [pkg_path=.] [reps=200] [n=3000] [threads=2] [out=adaptive_choice] [scenarios=null,diffuse,concentrated]
## Per replication: default (adaptive) fit, forced all-pooled fit, forced all-selected fit.
## Logs p.holm of each, the winning design string, and the scores of all designs.
## Adaptive margins under condition = "MW" are zmargin and equation; sample is always pooled.
args <- commandArgs(TRUE)
pkg     <- if (length(args) >= 1) args[1] else "."
reps    <- if (length(args) >= 2) as.integer(args[2]) else 200L
n       <- if (length(args) >= 3) as.integer(args[3]) else 3000L
threads <- if (length(args) >= 4) min(as.integer(args[4]), 16L) else 2L
out     <- if (length(args) >= 5) args[5] else "adaptive_choice"
scenarios <- if (length(args) >= 6) strsplit(args[6], ",")[[1]] else c("null", "diffuse", "concentrated")

suppressMessages(pkgload::load_all(pkg, quiet = TRUE))

dgp <- function(sc) {
  x <- rnorm(n); z <- sample(0:4, n, TRUE); u <- runif(n)
  p <- 0.2 + 0.1 * z
  bad <- x > 0
  if (sc == "diffuse")       p[bad] <- 0.6 - 0.05 * z[bad]                        # defiers at every z step
  if (sc == "concentrated")  p[bad] <- c(0.3, 0.4, 0.5, 0.35, 0.45)[z[bad] + 1]   # defiers at one z step
  d <- as.integer(u < p)
  y <- rnorm(n) + 0.5 * d + if (sc == "null") 0 else bad * d
  data.frame(y, x, d, z)
}

fit <- function(dat, ...) tryCatch(
  montest(y ~ x | d ~ z, data = dat, condition = "MW", progress = FALSE,
          Rparameters = list(num.threads = threads), ...),
  error = function(e) NULL)

p_of <- function(o) if (is.null(o)) NA_real_ else o$minp[["p.holm"]]

res <- list()
file_rds <- paste0(out, ".rds")

for (sc in scenarios) for (r in seq_len(reps)) {
  set.seed(1000 + r); dat <- dgp(sc)
  oa <- fit(dat)
  op <- fit(dat, pool = c("sample", "zmargin", "equation"), select = "none")
  os <- fit(dat, pool = "sample", select = c("zmargin", "equation"))
  sc_tab <- if (!is.null(oa) && !is.null(oa$adaptive_scores)) oa$adaptive_scores else NULL
  res[[length(res) + 1L]] <- data.frame(
    sc = sc, r = r,
    p_adaptive = p_of(oa), p_pooled = p_of(op), p_selected = p_of(os),
    design = if (is.null(oa) || is.null(oa$adaptive_design)) NA_character_ else oa$adaptive_design,
    scores = if (is.null(sc_tab)) NA_character_ else
      paste(sprintf("%s=%.3f", sc_tab$design, sc_tab$score), collapse = " | "),
    stringsAsFactors = FALSE)
  if (r %% 10L == 0L) saveRDS(do.call(rbind, res), file_rds)
}
res <- do.call(rbind, res)
saveRDS(res, file_rds)

res$kind <- ifelse(is.na(res$design), NA,
            ifelse(!grepl("selected_then_separate", res$design), "all_pooled",
            ifelse(!grepl("pooled_before_search", res$design), "all_selected", "mixed")))

cat("\nRejection rate (p.holm < 0.05):\n")
print(aggregate(cbind(p_adaptive, p_pooled, p_selected) ~ sc, res,
                function(p) mean(p < 0.05), na.action = na.pass))
cat("\nShare of replications by chosen design type:\n")
print(round(prop.table(table(res$sc, res$kind), 1), 3))
cat("\nRejection rate of adaptive by chosen design type:\n")
print(aggregate(p_adaptive ~ sc + kind, res, function(p) mean(p < 0.05)))
cat("\nMonte Carlo s.e. at rate 0.3 with", reps, "reps:", round(sqrt(.3 * .7 / reps), 3), "\n")
