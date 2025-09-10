#include <Rcpp.h>
#include <string>

// Finalize SNR from variance components and optional noise model
Rcpp::NumericVector estimateSNR(double totalVar,
                                double noiseRaw,
                                double eps = 1e-3,
                                const Rcpp::NumericVector& betas = Rcpp::NumericVector(),
                                const Rcpp::NumericVector& features = Rcpp::NumericVector(),
                                const std::string& noise_mode = "floor");

