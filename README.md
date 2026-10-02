# fastcramer

Fast, exact computation of the Cramér multivariate two-sample statistic
(Baringhaus & Franz 2004), the statistic behind `cramer::cramer.test()`, plus
a permutation test with optional strata and a bootstrap-versus-blocked-null
design for samples too large to test in full.

`cramer.test()` builds an (m+n)² lookup matrix with R-level loops. The same
statistic is three mean pairwise distances, which `fastcramer` evaluates with
BLAS matrix products and one fused C++ pass. On an Apple M-series laptop:

| per sample | cramer.test | fastcramer | speedup |
|---|---|---|---|
| 500 | 3.7 s | 0.005 s | 700× |
| 2,000 | 64 s | 0.04 s | 1,500× |
| 5,000 | ~400 s (extrapolated) | 0.3 s | ~1,300× |

Results agree with `cramer.test(..., just.statistic = TRUE)` to floating-point
rounding for all five kernels (`cramer`, `bahr`, `log`, `fraca`, `fracb`).

```r
library(fastcramer)
x <- matrix(rnorm(5000 * 768), 5000)
y <- matrix(rnorm(5000 * 768, mean = 0.01), 5000)
cramer_statistic(x, y)                       # the statistic, 0.3 s
cramer_test(x[1:300, ], y[1:300, ], 999)     # permutation test
cramer_bootstrap(x, y, replicates = 100, block = sample(1:20, 10000, TRUE))
```

`cramer_bootstrap()` implements the design used in Furnas, Gao, Yin & Wang,
*Partisan disparities in the production and uptake of science*: observed
statistics from repeated draws of `draw` observations per group, compared with
draws taken after permuting group labels within strata (`block`), summarised
by the exceedance probability and the percentile interval of the pairwise
differences.

## Installation

```r
# install.packages("devtools")
devtools::install_local("path/to/fastcramer")
```

Requires a C++ compiler (Xcode Command Line Tools on macOS, Rtools on Windows).
Linking against a fast BLAS (Apple vecLib, OpenBLAS, MKL) helps the pure-R
reference but matters little for the fused kernel, which is memory-bound.
