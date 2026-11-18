## Platform-independent (reproducible) summation.
##
## R's sum()/mean()/cumsum()/colSums() accumulate in `long double`, which is
## 80 bit on x86 and 64 bit on Apple Silicon, so results differ in the last
## bits across platforms; the forest/cutoff search then amplifies those
## differences into different cutoffs and p-values. These helpers make the sum
## independent of the accumulator: every element is split into a few
## components lying on fixed power-of-two grids q_k chosen so that ALL partial
## sums of the components of level k are integer multiples of q_k below 2^53 q_k.
## Such sums are exact in both double and long double, hence identical on every
## platform, whatever the summation order. The components' exact sums are then
## added in a fixed order using plain double arithmetic (IEEE, deterministic).
## Accuracy is about 2^-100 relative to n * max|x| (more than needed).
## Non-finite input falls back to base R (results are NA/Inf anyway).

.kdecomp <- function(x) {
  n <- length(x)
  bound <- max(abs(x)) * n
  total <- bound
  out <- list()
  r <- x
  for (k in 1:3) {
    if (!is.finite(bound) || bound <= 0) break
    q <- 2^(ceiling(log2(bound)) - 52)
    if (q < 2^-1000) break
    h <- round(r / q) * q
    out[[k]] <- h
    r <- r - h
    bound <- n * q / 2
    if (bound < total * 2^-75) break
  }
  out
}

ksum <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(0)
  if (!all(is.finite(x))) return(sum(x))
  if (n == 1L) return(as.numeric(x))
  parts <- .kdecomp(as.numeric(x))
  if (length(parts) == 0L) return(0)
  s <- 0
  for (p in parts) s <- s + sum(p)
  s
}

kcumsum <- function(x) {
  n <- length(x)
  if (n <= 1L) return(as.numeric(x))
  if (!all(is.finite(x))) return(cumsum(x))
  parts <- .kdecomp(as.numeric(x))
  if (length(parts) == 0L) return(numeric(n))
  s <- cumsum(parts[[1L]])
  for (p in parts[-1L]) s <- s + cumsum(p)
  s
}

kmean <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0L) return(NaN)
  ksum(x) / length(x)
}

kwmean <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    ok <- !is.na(x) & !is.na(w)
    x <- x[ok]; w <- w[ok]
  }
  if (length(x) == 0L) return(NaN)
  ksum(w * x) / ksum(w)
}

## Per-group sum (like tapply(x, g, sum)), deterministic: rowsum() adds in
## double in a fixed order. Returns a named numeric vector sorted by group.
ksum_by <- function(x, g) {
  r <- rowsum(as.numeric(x), g, reorder = TRUE)
  stats::setNames(as.numeric(r), rownames(r))
}

## Median that does not rely on long double mean() for the even-n average.
kmedian <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(NA_real_)
  s <- sort(x)
  h <- (n + 1L) %/% 2L
  if (n %% 2L == 1L) s[h] else (s[h] + s[h + 1L]) / 2
}

## Per-cluster sums in order of first appearance (same order as data.table
## `by=`), accumulated by rowsum() in plain double in a fixed order (GForce
## gsum uses long double).
gsum <- function(x, g) as.numeric(rowsum(as.numeric(x), g, reorder = FALSE))
