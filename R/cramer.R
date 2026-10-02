#' Cramér two-sample statistic, computed fast
#'
#' The multivariate two-sample statistic of Baringhaus and Franz (2004), as
#' implemented in [cramer::cramer.test()], computed from three mean pairwise
#' distances with BLAS matrix products and a fused C++ pass. It returns exactly
#' the same number as `cramer.test(x, y, just.statistic = TRUE, kernel = ...)`
#' (to floating-point rounding) at a small fraction of the cost, and scales to
#' thousands of observations per sample.
#'
#' @param x,y Numeric matrices with observations in rows and the same number
#'   of columns (dimensions). Vectors are treated as one-dimensional samples.
#' @param kernel One of `"cramer"` (phi(z) = sqrt(z)/2, the default of
#'   `cramer.test`), `"bahr"` (1 - exp(-z/2)), `"log"` (log(1 + z)),
#'   `"fraca"` (1 - 1/(1 + z)) or `"fracb"` (1 - 1/(1 + z)^2), applied to the
#'   squared Euclidean distance z between observations.
#' @param chunk Number of observations of `x` processed per block; controls
#'   peak memory (`chunk * nrow(y)` doubles), not the result.
#' @return A single number, the statistic.
#' @references Baringhaus, L., & Franz, C. (2004). On a new multivariate
#'   two-sample test. *Journal of Multivariate Analysis*, 88(1), 190-206.
#' @examples
#' x <- matrix(rnorm(200 * 10), 200)
#' y <- matrix(rnorm(150 * 10, mean = 0.2), 150)
#' cramer_statistic(x, y)
#' @export
cramer_statistic <- function(x, y, kernel = c("cramer", "bahr", "log", "fraca", "fracb"),
                             chunk = 2000L) {
  kernel <- match.arg(kernel)
  x <- as_obs_matrix(x); y <- as_obs_matrix(y)
  if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns (dimensions)")
  cramer_statistic_cols(t(x), t(y), kernel = kernel_code(kernel), chunk = as.integer(chunk))
}

#' Permutation test of equality of two multivariate distributions
#'
#' Pools the two samples, reassigns the group labels at random `replicates`
#' times (within the levels of `block`, if given), recomputes
#' [cramer_statistic()] each time, and reports the Monte Carlo p-value
#' `(1 + #\{T_perm >= T_obs\}) / (replicates + 1)`.
#'
#' This is a permutation test rather than the bootstrap-based critical value
#' of [cramer::cramer.test()]; with `block = NULL` the two agree in spirit
#' (exchangeability under the null) but not numerically.
#'
#' @inheritParams cramer_statistic
#' @param replicates Number of label permutations.
#' @param block Optional vector of length `nrow(x) + nrow(y)` (x first) giving
#'   strata; labels are permuted within strata, which preserves each stratum's
#'   group composition.
#' @param seed Optional integer seed.
#' @return An object of class `"cramer_test"`: a list with `statistic`, the
#'   vector of `null` statistics, `p.value`, `replicates`, `kernel`, and
#'   sample sizes `m` and `n`.
#' @examples
#' x <- matrix(rnorm(100 * 5), 100)
#' y <- matrix(rnorm(100 * 5, mean = 0.3), 100)
#' cramer_test(x, y, replicates = 199)
#' @export
cramer_test <- function(x, y, replicates = 999L, kernel = "cramer", block = NULL,
                        chunk = 2000L, seed = NULL) {
  kernel <- match.arg(kernel, c("cramer", "bahr", "log", "fraca", "fracb"))
  x <- as_obs_matrix(x); y <- as_obs_matrix(y)
  if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns (dimensions)")
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(x); n <- nrow(y)
  pooled <- t(rbind(x, y))                       # d x (m + n), observations in columns
  group  <- c(rep(TRUE, m), rep(FALSE, n))
  if (!is.null(block) && length(block) != m + n) stop("block must have length nrow(x) + nrow(y)")
  kc  <- kernel_code(kernel)
  obs <- cramer_statistic_cols(pooled[, group, drop = FALSE], pooled[, !group, drop = FALSE], kc, chunk)
  shuffler <- make_shuffler(block, m + n)
  null <- vapply(seq_len(replicates), function(i) {
    g <- shuffler(group)
    cramer_statistic_cols(pooled[, g, drop = FALSE], pooled[, !g, drop = FALSE], kc, chunk)
  }, numeric(1))
  structure(list(statistic = obs, null = null,
                 p.value = (1 + sum(null >= obs)) / (replicates + 1),
                 replicates = replicates, kernel = kernel, m = m, n = n),
            class = "cramer_test")
}

#' Bootstrap observed statistics against a blocked permutation null
#'
#' The design used for large samples where the full statistic is too costly:
#' each replicate draws `draw` observations with replacement from each group
#' and computes [cramer_statistic()] on the draws ("observed"), or first
#' permutes the group labels within `block` strata and then draws
#' ("permuted"). The two sets of replicates are compared directly: the
#' exceedance probability is the share of observed/permuted pairs in which the
#' observed statistic is larger, and the percentile interval summarises the
#' pairwise differences.
#'
#' @inheritParams cramer_test
#' @param draw Observations drawn per group per replicate. The default,
#'   `min(5000, nrow(x), nrow(y))`, uses the smaller group's size when it is
#'   below 5,000 so that small groups are not resampled far beyond their size.
#' @param replicates Number of observed and of permuted replicates.
#' @return An object of class `"cramer_bootstrap"`: a list with a data frame
#'   `replicates` (columns `type`, `statistic`), `draw`, `kernel`, and a
#'   `summary` list with the observed and permuted means, their difference,
#'   the 2.5th and 97.5th percentiles of the pairwise differences, and the
#'   exceedance probability.
#' @examples
#' x <- matrix(rnorm(600 * 8), 600)
#' y <- matrix(rnorm(400 * 8, mean = 0.2), 400)
#' yr <- sample(2000:2004, 1000, replace = TRUE)       # a stratum per observation
#' cramer_bootstrap(x, y, replicates = 50, draw = 200, block = yr)
#' @export
cramer_bootstrap <- function(x, y, replicates = 100L, draw = NULL, kernel = "cramer",
                             block = NULL, chunk = 2000L, seed = NULL) {
  kernel <- match.arg(kernel, c("cramer", "bahr", "log", "fraca", "fracb"))
  x <- as_obs_matrix(x); y <- as_obs_matrix(y)
  if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns (dimensions)")
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(x); n <- nrow(y)
  if (is.null(draw)) draw <- min(5000L, m, n)
  pooled <- t(rbind(x, y))
  group  <- c(rep(TRUE, m), rep(FALSE, n))
  if (!is.null(block) && length(block) != m + n) stop("block must have length nrow(x) + nrow(y)")
  kc <- kernel_code(kernel)
  draw_statistic <- function(g) {
    xi <- sample(which(g), draw, replace = TRUE)
    yi <- sample(which(!g), draw, replace = TRUE)
    cramer_statistic_cols(pooled[, xi, drop = FALSE], pooled[, yi, drop = FALSE], kc, chunk)
  }
  shuffler <- make_shuffler(block, m + n)
  observed <- vapply(seq_len(replicates), function(i) draw_statistic(group), numeric(1))
  permuted <- vapply(seq_len(replicates), function(i) draw_statistic(shuffler(group)), numeric(1))
  differences <- as.vector(outer(observed, permuted, "-"))
  structure(list(
    replicates = data.frame(type = rep(c("observed", "permuted"), each = replicates),
                            statistic = c(observed, permuted)),
    draw = draw, kernel = kernel, m = m, n = n,
    summary = list(observed_mean = mean(observed), permuted_mean = mean(permuted),
                   difference = mean(observed) - mean(permuted),
                   difference_interval = stats::quantile(differences, c(0.025, 0.975), names = FALSE),
                   exceedance = mean(differences > 0))),
    class = "cramer_bootstrap")
}

#' @export
print.cramer_test <- function(x, ...) {
  cat(sprintf("Cramer two-sample permutation test (kernel = %s)\n", x$kernel))
  cat(sprintf("  m = %d, n = %d, statistic = %.4f\n", x$m, x$n, x$statistic))
  cat(sprintf("  %d permutations, p-value = %.4f\n", x$replicates, x$p.value))
  invisible(x)
}

#' @export
print.cramer_bootstrap <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Cramer bootstrap vs blocked permutation null (kernel = %s, %d per group per draw)\n", x$kernel, x$draw))
  cat(sprintf("  observed mean %.4f, permuted mean %.4f, difference %.4f [%.4f, %.4f]\n",
              s$observed_mean, s$permuted_mean, s$difference, s$difference_interval[1], s$difference_interval[2]))
  cat(sprintf("  P(observed > permuted) = %.3f over %d x %d pairs\n", s$exceedance,
              nrow(x$replicates) / 2, nrow(x$replicates) / 2))
  invisible(x)
}

# ---- internals ---------------------------------------------------------------
as_obs_matrix <- function(a) {
  if (is.null(dim(a))) a <- matrix(a, ncol = 1)
  storage.mode(a) <- "double"
  if (anyNA(a)) stop("missing values are not allowed")
  a
}

kernel_code <- function(kernel) {
  match(kernel, c("cramer", "bahr", "log", "fraca", "fracb")) - 1L
}

# Returns a function that shuffles a logical label vector, within strata when
# `block` is given. Rows ordered by (stratum, random) receive the labels of rows
# ordered by (stratum), so each stratum keeps its label composition.
make_shuffler <- function(block, n) {
  if (is.null(block)) return(function(g) g[sample.int(n)])
  cell <- as.integer(factor(block))
  by_cell <- order(cell)
  function(g) { out <- g; out[order(cell, stats::runif(n))] <- g[by_cell]; out }
}
