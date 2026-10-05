#include <Rcpp.h>
using namespace Rcpp;

// Helper: index for 3D array stored as vector (column-major, dims N x T x K)
// i = 0..N-1, t=0..T-1, k=0..K-1
inline int idx(int i, int t, int k, int N, int T) {
  return i + t * N + k * N * T;
}

// [[Rcpp::export]]
NumericVector group_demean_formula_cpp(
  NumericVector Z,            // N x T x (K+1)
  IntegerVector unit_group,   // length N
  IntegerVector time_group,   // length T
  IntegerVector covar_group   // length K+1
) {
  IntegerVector dim = Z.attr("dim");
  int N = dim[0], T = dim[1], Kp1 = dim[2];
  
  NumericVector Z_out(clone(Z)); // output same dims
  
  // Pre-count group sizes to speed up
  // unit_group: count units per group
  std::map<int,int> unit_count;
  for (int i = 0; i < N; ++i) unit_count[unit_group[i]]++;
  
  // time_group: count times per group
  std::map<int,int> time_count;
  for (int t = 0; t < T; ++t) time_count[time_group[t]]++;
  
  // covar_group: count covars per group
  std::map<int,int> covar_count;
  for (int k = 0; k < Kp1; ++k) covar_count[covar_group[k]]++;
  
  for (int i = 0; i < N; ++i) {
    int g_i = unit_group[i];
    for (int t = 0; t < T; ++t) {
      int m_t = time_group[t];
      for (int k = 0; k < Kp1; ++k) {
        int l_k = covar_group[k];
        
        // Compute averages:
        // \bar{z}_{g_i t k}
        double sum_gitk = 0.0;
        int cnt_gi = 0;
        for (int j = 0; j < N; ++j) {
          if (unit_group[j] == g_i) {
            sum_gitk += Z[idx(j,t,k,N,T)];
            cnt_gi++;
          }
        }
        double avg_gitk = cnt_gi > 0 ? sum_gitk / cnt_gi : 0.0;
        
        // \bar{z}_{i m_t k}
        double sum_imtk = 0.0;
        int cnt_mt = 0;
        for (int s = 0; s < T; ++s) {
          if (time_group[s] == m_t) {
            sum_imtk += Z[idx(i,s,k,N,T)];
            cnt_mt++;
          }
        }
        double avg_imtk = cnt_mt > 0 ? sum_imtk / cnt_mt : 0.0;
        
        // \bar{z}_{g_i m_t k}
        double sum_gimtk = 0.0;
        int cnt_gimt = 0;
        for (int j = 0; j < N; ++j) {
          if (unit_group[j] == g_i) {
            for (int s = 0; s < T; ++s) {
              if (time_group[s] == m_t) {
                sum_gimtk += Z[idx(j,s,k,N,T)];
                cnt_gimt++;
              }
            }
          }
        }
        double avg_gimtk = cnt_gimt > 0 ? sum_gimtk / cnt_gimt : 0.0;
        
        // \bar{z}_{i t l_k}
        double sum_itlk = 0.0;
        int cnt_lk = 0;
        for (int q = 0; q < Kp1; ++q) {
          if (covar_group[q] == l_k) {
            sum_itlk += Z[idx(i,t,q,N,T)];
            cnt_lk++;
          }
        }
        double avg_itlk = cnt_lk > 0 ? sum_itlk / cnt_lk : 0.0;
        
        // \bar{z}_{g_i t l_k}
        double sum_gitlk = 0.0;
        int cnt_gitlk = 0;
        for (int j = 0; j < N; ++j) {
          if (unit_group[j] == g_i) {
            for (int q = 0; q < Kp1; ++q) {
              if (covar_group[q] == l_k) {
                sum_gitlk += Z[idx(j,t,q,N,T)];
                cnt_gitlk++;
              }
            }
          }
        }
        double avg_gitlk = cnt_gitlk > 0 ? sum_gitlk / cnt_gitlk : 0.0;
        
        // \bar{z}_{i m_t l_k}
        double sum_imtlk = 0.0;
        int cnt_imtlk = 0;
        for (int s = 0; s < T; ++s) {
          if (time_group[s] == m_t) {
            for (int q = 0; q < Kp1; ++q) {
              if (covar_group[q] == l_k) {
                sum_imtlk += Z[idx(i,s,q,N,T)];
                cnt_imtlk++;
              }
            }
          }
        }
        double avg_imtlk = cnt_imtlk > 0 ? sum_imtlk / cnt_imtlk : 0.0;
        
        // \bar{z}_{g_i m_t l_k}
        double sum_gimtlk = 0.0;
        int cnt_gimtlk = 0;
        for (int j = 0; j < N; ++j) {
          if (unit_group[j] == g_i) {
            for (int s = 0; s < T; ++s) {
              if (time_group[s] == m_t) {
                for (int q = 0; q < Kp1; ++q) {
                  if (covar_group[q] == l_k) {
                    sum_gimtlk += Z[idx(j,s,q,N,T)];
                    cnt_gimtlk++;
                  }
                }
              }
            }
          }
        }
        double avg_gimtlk = cnt_gimtlk > 0 ? sum_gimtlk / cnt_gimtlk : 0.0;
        
        // Apply formula:
        double val = Z[idx(i,t,k,N,T)] 
             - avg_gitk - avg_imtk + avg_gimtk - (avg_itlk - avg_gitlk - avg_imtlk + avg_gimtlk);

                     
        
        Z_out[idx(i,t,k,N,T)] = val;
      }
    }
  }
  
  return Z_out;
}
