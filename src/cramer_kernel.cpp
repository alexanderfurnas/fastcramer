// Fused Cramér two-sample statistic (Baringhaus & Franz 2004).
//
//   T = m n / (m + n) * [ 2/(m n) S(X, Y) - 1/m^2 S(X, X) - 1/n^2 S(Y, Y) ],
//   S(A, B) = sum_{i, j} phi( |a_i - b_j|^2 ),
//
// with phi one of the kernels of the 'cramer' package:
//   cramer: sqrt(z) / 2        bahr: 1 - exp(-z / 2)      log: log(1 + z)
//   fraca:  1 - 1 / (1 + z)    fracb: 1 - 1 / (1 + z)^2
//
// Each cross-product block is computed by R's BLAS (dgemm); the squared
// distance, kernel and sum are accumulated in one pass over the block, so no
// large temporaries are materialised. Observations are passed as columns
// (d x n, column-major) so a subset of observations is a set of contiguous
// columns.
#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

enum Kernel { CRAMER = 0, BAHR = 1, LOG = 2, FRACA = 3, FRACB = 4 };

static inline double phi(double z, int k) {
  switch (k) {
    case CRAMER: return std::sqrt(z) / 2.0;
    case BAHR:   return 1.0 - std::exp(-z / 2.0);
    case LOG:    return std::log1p(z);
    case FRACA:  return 1.0 - 1.0 / (1.0 + z);
    default:     { const double q = 1.0 / (1.0 + z); return 1.0 - q * q; }
  }
}

static std::vector<double> sq_norms(const double* A, int d, int n) {
  std::vector<double> out(n);
  for (int j = 0; j < n; ++j) {
    const double* col = A + (size_t)j * d;
    double s = 0.0;
    for (int i = 0; i < d; ++i) s += col[i] * col[i];
    out[j] = s;
  }
  return out;
}

static double phi_sum(const double* A, int m, const double* B, int n, int d,
                      const std::vector<double>& a2, const std::vector<double>& b2,
                      int chunk, int kernel) {
  std::vector<double> G((size_t)chunk * n);
  const char tA = 'T', tB = 'N';
  const double one = 1.0, zero = 0.0;
  double total = 0.0;
  for (int s = 0; s < m; s += chunk) {
    const int rows = std::min(chunk, m - s);
    F77_CALL(dgemm)(&tA, &tB, &rows, &n, &d, &one, A + (size_t)s * d, &d, B, &d,
                    &zero, G.data(), &rows FCONE FCONE);
    double acc = 0.0;
    for (int j = 0; j < n; ++j) {
      const double* g = G.data() + (size_t)j * rows;
      const double bj = b2[j];
      for (int i = 0; i < rows; ++i) {
        double z = a2[s + i] + bj - 2.0 * g[i];
        if (z < 0.0) z = 0.0;                     // rounding below zero
        acc += phi(z, kernel);
      }
    }
    total += acc;
  }
  return total;
}

//' @rdname cramer_statistic
//' @param xt,yt Numeric matrices with observations in COLUMNS (d x m and d x n).
//' @param kernel Integer kernel code: 0 cramer, 1 bahr, 2 log, 3 fraca, 4 fracb.
//' @param chunk Rows of the cross-product block held in memory at once.
//' @keywords internal
// [[Rcpp::export]]
double cramer_statistic_cols(NumericMatrix xt, NumericMatrix yt, int kernel = 0, int chunk = 2000) {
  const int d = xt.nrow();
  if (yt.nrow() != d) stop("x and y must have the same number of dimensions");
  if (kernel < 0 || kernel > 4) stop("unknown kernel code");
  if (chunk < 1) stop("chunk must be positive");
  const int m = xt.ncol(), n = yt.ncol();
  if (m < 1 || n < 1) stop("each sample needs at least one observation");
  const double* X = xt.begin();
  const double* Y = yt.begin();
  const std::vector<double> x2 = sq_norms(X, d, m), y2 = sq_norms(Y, d, n);
  const double sxy = phi_sum(X, m, Y, n, d, x2, y2, chunk, kernel);
  const double sxx = phi_sum(X, m, X, m, d, x2, x2, chunk, kernel);
  const double syy = phi_sum(Y, n, Y, n, d, y2, y2, chunk, kernel);
  const double mm = m, nn = n;
  return mm * nn / (mm + nn) * (2.0 * sxy / (mm * nn) - sxx / (mm * mm) - syy / (nn * nn));
}
