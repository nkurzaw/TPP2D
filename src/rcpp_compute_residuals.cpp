
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_compute_residuals(
  NumericVector temp_i,
  double zeta,
  double slope,
  double beta_max,
  NumericVector beta_0,
  NumericVector alpha,
  NumericVector log_conc,
  NumericVector log2_value) {
    int n = temp_i.size();
    NumericVector out(n);
    for(int i = 0; i < n; ++i){
      int temp_index = temp_i[i] - 1;
      out[i] = pow(
        (beta_0[temp_index] + 
           ((alpha[temp_index] * beta_max) / 
              (1 + exp( -slope * (log_conc[i] - zeta))))) - log2_value[i], 2);
      }
    return out;
    }
