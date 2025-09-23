#include <Rcpp.h>
#include <algorithm> // std::max
#include <cmath>     // std::log2
#include <limits>
#include <string>

//' @title Estimate a Time Series Signal to Noise Ratio (SNR) given
//' total and noise variance components and (optionally) a background Noise model.
//'
//' @description
//' Given Time series variance components (total variance, raw noise variance),
//' compute the final noise, signal, and SNR. For noise
//' estimation use one of:
//' \enumerate{
//'   \item **A Raw noise variance** e.g the result of \code{estimateNoise}
//'   \item **A Background noise model** based on linear predictors
//'   \item **A Flooring rule**: final noise = max(raw, background)
//' }
//' The Background noise model is a *general* linear noise baseline:
//' \code{baseline = sum_i beta[i] * feat[i]}.
//'
//' @param totalVar Total time series variance. Typically calculated with \code{estimateNoise}
//' @param noiseRaw noise variance estimate. Typically calculated with \code{estimateNoise}
//' @param eps Numeric. Lower bound applied to both noise and signal.
//' @param betas,features Coefficients (\eqn{betas}) and predictors (\eqn{features})
//'   for the background noise model. Ignored when \code{noise_mode = "raw"}.
//'   Both vectors must have the same length when used. The predictors form a
//'   *design vector*: include a constant \code{1} if the model has an intercept.
//'   For example, for \eqn{baseline = b0 + b1 * mean}, use
//'   \code{betas = c(b0, b1)}, \code{features = c(1, mean)}.
//' @param noise_mode Character scalar: one of \code{"raw"}, \code{"model"}, \code{"floor"}.
//'        - \code{"raw"}   → use \code{noise_raw}
//'        - \code{"model"} → use \code{baseline = sum(betas * features)}
//'        - \code{"floor"} → use \code{max(noise_raw, baseline)}
//'
//' @return Named numeric vector with elements:
//'         - snr (log2 scale)
//'         - signal
//'         - noise final estimate used for signal and snr calculation
//'         - noise baseline (sum(betas * features); NA if not applied)
//'         - noise raw
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector estimateSNR(double totalVar,
                                double noiseRaw,
                                double eps,
                                const Rcpp::NumericVector& betas,
                                const Rcpp::NumericVector& features,
                                const std::string& noise_mode = "floor") {
    Rcpp::NumericVector out(5, NA_REAL);
    out.names() = Rcpp::CharacterVector::create("snr", "signal", "noise", "baseline", "raw");

    if (!R_finite(totalVar)) {
        return out;
    }

    // For raw mode, both totalVar and noiseRaw must be finite
    if (noise_mode == "raw" && !R_finite(noiseRaw)) {
        return out;
    }

    // Compute baseline if needed
    double baseline = 0.0;
    const bool need_baseline = (noise_mode == "model" || noise_mode == "floor");

    if (need_baseline) {
        if (betas.size() != features.size() || betas.size() == 0) {
            Rcpp::stop("For noise_mode='%s', 'betas' and 'features' must be the same non-zero length.", noise_mode.c_str());
        }

        // Fail on misconfigured betas
        for (int i = 0; i < betas.size(); i++) {
            if (!R_finite(betas[i])) {
                Rcpp::stop("Non-finite value in 'betas' for noise_mode='%s'.", noise_mode.c_str());
            }
        }

        // In case of non-finite features return NA vector
        for (int i = 0; i < features.size(); i++) {
            if (!R_finite(features[i])) {
                return out;  // baseline required but cannot be computed
            }
        }

        for (int i = 0; i < betas.size(); i++) {
            baseline += betas[i] * features[i];
        }

        if (!R_finite(baseline)) {
            // cannot proceed if baseline required but not finite
            return out;
        }
    }

    // Select noise according to mode
    double noiseV;
    if (noise_mode == "raw") {
        noiseV = noiseRaw;
    } else if (noise_mode == "model") {
        noiseV = baseline;
    } else if (noise_mode == "floor") {
        if (!R_finite(noiseRaw)) {
            noiseV = baseline;
        } else {
            noiseV = std::max(noiseRaw, baseline);
        }
    } else {
        Rcpp::stop("Invalid noise_mode: '%s'. Use 'raw', 'model', or 'floor'.", noise_mode.c_str());
    }

    // Apply eps to noise and ensure finite
    if (!R_finite(noiseV)) {
        return out;
    }
    if (noiseV < eps) {
        noiseV = eps;
    }

    // Signal and SNR
    double signalV = totalVar - noiseV;
    if (!R_finite(signalV)) {
        return out;
    }
    if (signalV < eps) {
        signalV = eps;
    }

    const double snr = std::log2(signalV / noiseV);

    out[0] = snr;
    out[1] = signalV;
    out[2] = noiseV;
    out[3] = need_baseline ? baseline : NA_REAL;
    out[4] = noiseRaw;
    return out;
}
