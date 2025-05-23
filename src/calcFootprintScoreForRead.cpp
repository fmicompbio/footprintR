#include <Rcpp.h>

//' Calculate footprinting score for a defined footprint and a single read.
//'
//' This function calculates the footprinting scores corresponding to a
//' footprint in the form of a weight vector \code{wgt} for an individual
//' read. The score is based on a cross-correlation of the modification
//' probabilities in \code{pmod} (centered by subtracting \code{0.5}) with
//' \code{wgt}, weighted by the minimum of \code{minweight} and the elements
//' of \code{pmod}.
//'
//' @param pos Integer vector with positions (genomic coordinates) of modified
//'     bases.
//' @param pmod Numeric vector with modification probabilities of modified
//'     bases.
//' @param wgt Numeric vector with weights that define the footprint to score.
//'     Typically centered at zero.
//' @param minconf Numeric scalar giving the minimal confidence of a base to be
//'     included in the calculation.
//' @param minweight Numeric scalar giving the minimal weight for position.
//'
//' @returns
//' A \code{data.frame} with columns \code{pos}, \code{pmod} and \code{score},
//'     in (continuous) base space from \code{min(pos)} to \code{max(pos)}.
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame calcFootprintScoreForRead(Rcpp::IntegerVector pos,
                                          Rcpp::NumericVector pmod,
                                          Rcpp::NumericVector wgt,
                                          double minconf = 0.7,
                                          double minweight = 0.05) {
    // mask low confidence probabilities
    pmod = Rcpp::ifelse(Rcpp::abs(pmod - 0.5) > minconf - 0.5, pmod, NA_REAL);

    // center pmod
    pmod = pmod - 0.5;

    // expand to base space
    Rcpp::IntegerVector posall = Rcpp::seq(Rcpp::min(pos), Rcpp::max(pos));
    Rcpp::NumericVector pmodall = Rcpp::rep(NA_REAL, posall.size());
    pmodall[Rcpp::match(pos, posall) - 1] = pmod;
    Rcpp::LogicalVector hasNA = Rcpp::is_na(pmodall);

    // calculate scores
    Rcpp::NumericVector scores = Rcpp::rep(NA_REAL, posall.size());

    if (posall.size() >= wgt.size()) {
        double currentScore = 0.0;
        double currentTotalWeight = 0.0;
        double currentWeight = 0.0;
        for (size_t i = (unsigned)((wgt.size() - 1) / 2), s = 0;
             i < (unsigned)(posall.size() - (wgt.size() / 2));
             i++, s++) {

            currentScore = 0.0;
            currentTotalWeight = 0.0;
            currentWeight = 0.0;
            for (size_t j = s; j < s + wgt.size(); j++) {
                if (hasNA[j] == false) {
                    currentWeight = (pmodall[j] + 0.5 > minweight ?
                                         pmodall[j] + 0.5 :
                                         minweight);
                    currentScore += pmodall[j] * wgt[j - s] * currentWeight;
                    currentTotalWeight += currentWeight;
                }
            }
            scores[i] = (currentTotalWeight == 0.0 ? NA_REAL : currentScore / currentTotalWeight);
        }
    }

    // create return value
    Rcpp::DataFrame df = Rcpp::DataFrame::create(Rcpp::_["pos"] = posall,
                                                 Rcpp::_["pmod"] = pmodall,
                                                 Rcpp::_["score"] = scores);
    return(df);
}
