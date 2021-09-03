data{
  int<lower = 1> n_observed;
  int<lower = 0> n_missing;
  int<lower = 1> nr;
  int<lower = 1> ni;
  int<lower = 1> np;

  // observed data
  real x[n_observed];
  int<lower = 1, upper = nr> idx_rater  [n_observed];
  int<lower = 1, upper = ni> idx_item   [n_observed];
  int<lower = 1, upper = np> idx_patient[n_observed];

  // missing data
  int<lower = 1, upper = nr> idx_rater_missing  [n_missing];
  int<lower = 1, upper = ni> idx_item_missing   [n_missing];
  int<lower = 1, upper = np> idx_patient_missing[n_missing];

  int prior_only;
  int debug;
}
transformed data{

  // constant hyperparameters, see Anders, Oravecz, and Batchelder (2014) p. 5
  real mu_log_lambda = 0.0;
  real mu_log_a      = 0.0;
  real mu_b          = 0.0;

  // QR trick
  matrix[nr, nr] A = diag_matrix(rep_vector(1,nr));
  matrix[nr, nr - 1] A_qr;
  for (i in 1:nr - 1) A[nr,i] = -1;
  A[nr, nr] = 0;
  A_qr = qr_Q(A)[ , 1:(nr-1)];

}
parameters{
  // See Anders, Oravecz, and Batchelder (2014) p. 5

  // hyperparameters
  real mu_lt;
  real mu_log_E;
  real<lower = 0> sd_lt;
  real<lower = 0> sd_log_lambda;
  real<lower = 0> sd_log_E;
  real<lower = 0> sd_log_a;
  real<lower = 0> sd_b;

  // parameters
  matrix [ni, np] lt;         // item location values
  matrix [ni, np] log_lambda; // item difficulty

  // see https://mc-stan.org/docs/2_18/stan-users-guide/parameterizing-centered-vectors.html
  vector [nr-1]     log_E_raw;  // rater competency
  vector [nr-1]     log_a_raw;  // rater scaling

  vector [nr]     b;          // rater shifting

}
transformed parameters{

  vector [nr]     log_E = A_qr * log_E_raw;      // rater competency
  vector [nr]     log_a = A_qr * log_a_raw;      // rater scaling

}
model{

  log_E_raw ~ normal(0, inv(sqrt(1 - inv(nr))));
  log_a_raw ~ normal(0, inv(sqrt(1 - inv(nr))));

  // hyperpriors
  mu_lt         ~ normal(0, 0.25);
  mu_log_E      ~ normal(0, .01);
  sd_lt         ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_log_lambda ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_log_E      ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_log_a      ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_b          ~ gamma(1, 1);//inv_gamma(0.01, 0.01);

  // priors
  for (p in 1:np) {
    lt[, p]         ~ normal(mu_lt,         sd_lt);
    log_lambda[, p] ~ normal(mu_log_lambda, sd_log_lambda);
  }
  // log_E      ~ normal(mu_log_E, sd_log_E);
  // log_a      ~ normal(mu_log_a, sd_log_a);
  b ~ normal(mu_b, sd_b);

  if (!prior_only) {
    // likelihood in long form
    for (o in 1:n_observed) {
      int r = idx_rater[o];
      int i = idx_item[o];
      int p = idx_patient[o];

      real location = exp(log_a[r]) * lt[i, p] + b[r];
      // real scale    = exp(2.0 * log_a[r] + log_E[r] + log_lambda[i, p]);
      real scale    = exp(log_E[r] + log_lambda[i, p]);

      // real location = lt[i, p] + b[r];
      // real scale    = exp(log_lambda[i, p]);

      // real location = lt[i, p];// + b[r];
      // real scale    = 1.0;

      x[o] ~ normal(location, scale);
    }
  }
  // marginalize over missing values -- yields a constant so left out?
  // for (m in 1:nMissing) {
  //   vector[nc] logProbs;
  //   int r = idxRaterMissing[m];
  //   int i = idxItemMissing[m];
  //   int p = idxPatientMissing[m];
  //   real scale = kappa[i, p] / zeta[r];
  //   real location = lt[i, p];
  //   vector[nt] delta = alpha[r] * lambda + beta[r];
  //   for (c in 1:nc) {
  //     logProbs[c] = ordered_logistic_lpmf(c | location / scale, delta ./ scale);
  //   }
  //   target += log_sum_exp(2*logProbs);
  //   //   log(probs) +
  //   //   categorical_lpmf(1:nc | probs)
  //   // );
  // }
}
generated quantities {

  real             mus[debug ? np : 0, debug ? ni : 0, debug ? nr : 0];
  real<lower = 0>  sds[debug ? np : 0, debug ? ni : 0, debug ? nr : 0];

  if (debug) {

    for (o in 1:n_observed) {
      int r = idx_rater[o];
      int i = idx_item[o];
      int p = idx_patient[o];

      mus[p, i, r] = exp(log_a[r]) * lt[i, p] + b[r];
      sds[p, i, r] = exp(log_E[r] + log_lambda[i, p]);
    }
  }

//   // impute missing values by sampling from the posterior predictive distribution
//   int<lower = 1, upper = nc> xPred[nMissing];
//   for (m in 1:nMissing) {
//     int r = idxRaterMissing[m];
//     int i = idxItemMissing[m];
//     int p = idxPatientMissing[m];
//     real scale = kappa[i, p] / zeta[r];
//     real location = lt[i, p];
//     vector[nt] delta = alpha[r] * lambda + beta[r];
//     xPred[m] = ordered_logistic_rng(location / scale, delta ./ scale);
//   }
}
