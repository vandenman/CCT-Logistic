#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector get_probs_ordered_logistic_cpp(const NumericVector thresholds, const double location, const double scale) {

  const int nc = thresholds.size() + 1;
  NumericVector prob(nc);

  double vals0, vals1;

  vals0 = 1.0 / (1.0 + exp((thresholds[0] - location) / scale));
  prob[0] = 1.0 - vals0;
  for (size_t i = 1; i < prob.size()-1; i++)
  {
    vals1 = 1.0 / (1.0 + exp((thresholds[i] - location) / scale));
    prob[i] = vals0 - vals1;
    vals0 = vals1;
  }

  prob[nc - 1] = vals1;

  return prob;
}


/*** R
thresholds <- 1:18
location <- 3
scale <- 2

get_probs_ordered_logistic_cpp(thresholds, location, scale)
CCTLogistic:::get_probs_ordered_logistic(thresholds, location, scale)

bch <- bench::mark(
  cpp = get_probs_ordered_logistic_cpp(thresholds, location, scale),
  R = CCTLogistic:::get_probs_ordered_logistic(thresholds, location, scale)
)
summary(bch, relative = TRUE)

*/
