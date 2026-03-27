# montest

Tools for detecting violations of monotonicity and LATE assumptions using machine learning.

## Installation

You can install the development version from GitHub:

```r
install.packages("devtools")"
devtools::install_github("martin-andresen/montest")
```

## Overview

`montest` searches for violations of monotonicity and LATE assumptions in subsets of the data and over margins using sample splitting and machine learning methods.

`montestplot` provides visualization tools for the detected violations.

## Example

```r
library(montest)
data=fct_datasim(setup="A",dgp=2,n=3000)

# Simple monotonicity test
out <- montest(
   data = data,
   D = "D",
   Z = "Z",
   X = c("Xvar1", "Xvar2", "Xvar3"),
   test = "simple",
   testtype = "forest")
   
# Test multiple conditions, pooling evidence
 out2 <- montest(
   data = data,
   D = "D",
   Z = "Z",
   Y="Y",
   X = c("Xvar1", "Xvar2", "Xvar3"),
   test = c("simple","BP","MW"))
 }

montestplot(out)
```

## Main functions

- `montest()` ??? detects violations of monotonicity / LATE assumptions
- `montestplot()` ??? visualizes results

## Status

This is a research package under active development. Beware that with large data,
estimation may be computationally heavy. Reduce computation time by reducing Dsubsets,
Zsubsets, Ysubsets, using conditions that are simpler to test (e.g. test="simple" or 
test="MW" rather than test="BP"), or reduce the number of inner folds.

The package is developed as part of the paper "Testing the Monotonicity Assumption in Instrumental Variable Models", joint with Tymon S??oczy??ski and Martn Huber

## Author
Martin Eckhoff Andresen, Department of Economics, University of Oslo

Please get in touch with comments or feedback m.e.andresen@econ.uio.no
