# montest

Tools for detecting violations of monotonicity and LATE assumptions using machine learning.

## Installation

You can install the development version from GitHub:

```r
install.packages("devtools")
pak::pak("martin-andresen/montest")
```

**Version 1.0.0.** The version number does not identify the exact code; the GitHub commit ID does. Find it in `montest_version()$sha`, in `$version$sha` of the results of `montest()`, `seqtest()` and the `*_seeds()` functions, or in the `RemoteSha` field of `packageDescription("montest")`. The commit ID is recorded when installing from GitHub as above, not when installing from a local folder. Look up a commit at `https://github.com/martin-andresen/montest/commit/<ID>`.

## Overview

`montest` searches for violations of monotonicity and LATE assumptions in subsets of the data and over margins using sample splitting and machine learning methods.

`montestplot` provides visualization tools for the detected violations.

## Example

```r
library(montest)
data <- fct_datasim(setup = "A", n = 3000, J = 1, K = 1)

# fml is a fixest-style formula: Y ~ X | FE | D ~ Z
# (the FE part is optional and can be omitted)
fml <- Y ~ Xvar1 + Xvar2 + Xvar3 | D ~ Z

# Simple monotonicity test
out <- montest(
   fml = fml,
   data = data,
   condition = "simple")

# Test multiple conditions, pooling evidence
out2 <- montest(
   fml = fml,
   data = data,
   condition = c("simple", "KR", "MW"))
```

## Main functions

- `montest()` — detects violations of monotonicity / LATE assumptions
- `seqtest()` — wrapper for montest(), to test violations of instrument validity with multiple treatment parts with a sequential structure
- `montestplot()` — visualizes results

## Status

This is a research package under active development. Beware that with large data,
estimation may be computationally heavy. Reduce computation time by reducing Dsubsets,
Zsubsets, Ysubsets, using conditions that are simpler to test (e.g. condition="simple" or 
condition="MW" rather than condition="KR"), or reduce the number of inner folds.

The package is developed as part of the paper "Testing the Monotonicity Assumption in Instrumental Variable Models", joint with Tymon Słoczyński and Martin Huber

## Author
Martin Eckhoff Andresen, Department of Economics, University of Oslo

Please get in touch with comments or feedback m.e.andresen@econ.uio.no
