#include <Rcpp.h>
using namespace Rcpp;

// Indexing helper for z[i, t, k] in R's column-major layout
inline int idx(int i, int t, int k, int N, int T) {
  return i + t * N + k * N * T;
}

// ------------------------------
// [[Rcpp::export]]
NumericMatrix pseudo_dist_unit(NumericVector z) {
  IntegerVector dim = z.attr("dim");
  int N = dim[0], T = dim[1], K1 = dim[2];
  double norm = T * K1;

  NumericMatrix dist_matrix(N, N);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i == j) {
        dist_matrix(i, j) = 0.0;
      } else {
        double max_val = -1e10;
        for (int r = 0; r < N; ++r) {
          if (r == i || r == j) continue;
          double sum = 0.0;
          for (int t = 0; t < T; ++t)
            for (int k = 0; k < K1; ++k)
              sum += (z[idx(i, t, k, N, T)] - z[idx(j, t, k, N, T)]) * z[idx(r, t, k, N, T)];
          double abs_val = std::abs(sum / norm);
          if (abs_val > max_val) max_val = abs_val;
        }
        dist_matrix(i, j) = max_val;
      }
    }
  }

  return dist_matrix;
}

// ------------------------------
// [[Rcpp::export]]
NumericMatrix pseudo_dist_time(NumericVector z) {
  IntegerVector dim = z.attr("dim");
  int N = dim[0], T = dim[1], K1 = dim[2];
  double norm = N * K1;

  NumericMatrix dist_matrix(T, T);

  for (int t = 0; t < T; ++t) {
    for (int s = 0; s < T; ++s) {
      if (t == s) {
        dist_matrix(t, s) = 0.0;
      } else {
        double max_val = -1e10;
        for (int tau = 0; tau < T; ++tau) {
          if (tau == t || tau == s) continue;
          double sum = 0.0;
          for (int i = 0; i < N; ++i)
            for (int k = 0; k < K1; ++k)
              sum += (z[idx(i, t, k, N, T)] - z[idx(i, s, k, N, T)]) * z[idx(i, tau, k, N, T)];
          double abs_val = std::abs(sum / norm);
          if (abs_val > max_val) max_val = abs_val;
        }
        dist_matrix(t, s) = max_val;
      }
    }
  }

  return dist_matrix;
}

// ------------------------------
// [[Rcpp::export]]
NumericMatrix pseudo_dist_covariate(NumericVector z) {
  IntegerVector dim = z.attr("dim");
  int N = dim[0], T = dim[1], K1 = dim[2];
  double norm = N * T;

  NumericMatrix dist_matrix(K1, K1);

  for (int k = 0; k < K1; ++k) {
    for (int q = 0; q < K1; ++q) {
      if (k == q) {
        dist_matrix(k, q) = 0.0;
      } else {
        double max_val = -1e10;
        for (int u = 0; u < K1; ++u) {
          if (u == k || u == q) continue;
          double sum = 0.0;
          for (int i = 0; i < N; ++i)
            for (int t = 0; t < T; ++t)
              sum += (z[idx(i, t, k, N, T)] - z[idx(i, t, q, N, T)]) * z[idx(i, t, u, N, T)];
          double abs_val = std::abs(sum / norm);
          if (abs_val > max_val) max_val = abs_val;
        }
        dist_matrix(k, q) = max_val;
      }
    }
  }

  return dist_matrix;
}
