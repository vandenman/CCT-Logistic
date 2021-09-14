functions{
  matrix getQR(int K) {
    // see https://mc-stan.org/docs/2_18/stan-users-guide/parameterizing-centered-vectors.html
    matrix[K, K] A = diag_matrix(rep_vector(1, K));
    matrix[K, K-1] A_qr;
    for (i in 1:K-1) {
      A[K, i] = -1;
    }
    A[K, K] = 0;
    A_qr = qr_Q(A)[, 1:(K-1)];
    return A_qr;
  }
}
data{
  int<lower = 1> n_observed;
  int<lower = 0> n_missing;
  int<lower = 2> nr;
  int<lower = 2> ni;
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

  // constant hyperparameters, see Anders, Oravecz, and Batchelder (2014) p. 5
  real mu_log_lambda;
  real mu_log_a;
  real mu_b;

  int prior_only;
  int debug;
  int vary_lambda_across_patients;
}
transformed data{

  matrix[nr, nr - 1] QR_nr = getQR(nr);
  matrix[ni, ni - 1] QR_ni = getQR(ni);

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
  matrix [ni-1, np] lt_raw;         // item location values
  matrix [ni-1, vary_lambda_across_patients ? np : 1] log_lambda_raw; // item difficulty

  vector [nr-1]   log_a_raw; // rater scaling
  vector [nr]     log_E_raw; // rater competency
  vector [nr]     b_raw;     // rater shifting



}
transformed parameters{

  // noncentered to centered
  vector [nr] b = mu_b + sd_b * b_raw;
  vector [nr] log_E = mu_log_E + sd_log_E * log_E_raw;

  // vector [nr]     log_E = mu_log_E + sd_log_E * QR_nr * log_E_raw;      // rater competency
  vector [nr]     log_a = mu_log_a + sd_log_a * QR_nr * log_a_raw;      // rater scaling

  matrix [ni, vary_lambda_across_patients ? np : 1] log_lambda;
  matrix [ni, np] lt;

  for (p in 1:np) {
    lt[, p] = mu_lt + sd_lt * QR_ni * lt_raw[, p];
  }
  if (vary_lambda_across_patients) {
    for (p in 1:np) {
      log_lambda[, p] = mu_log_lambda + sd_log_lambda * QR_ni * log_lambda_raw[, p];
    }
  } else {
    log_lambda[, 1] = mu_log_lambda + sd_log_lambda * QR_ni * log_lambda_raw[, 1];
  }

}
model{

  // hyperpriors
  mu_lt         ~ normal(0, 0.25);
  mu_log_E      ~ normal(0, .01);
  sd_lt         ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_log_lambda ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_log_E      ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_log_a      ~ gamma(1, 1);//inv_gamma(0.01, 0.01);
  sd_b          ~ gamma(1, 1);//inv_gamma(0.01, 0.01);

  // priors
  to_vector(log_lambda_raw) ~ normal(0, inv(sqrt(1 - inv(ni))));
  to_vector(lt_raw)         ~ normal(0, inv(sqrt(1 - inv(ni))));

  // if (vary_lambda_across_patients) {
  //   for (p in 1:np) {
  //     log_lambda_raw[, p] ~ normal(0, inv(sqrt(1 - inv(ni))));
  //   }
  // } else {
  //   log_lambda_raw[, 1] ~ normal(0, inv(sqrt(1 - inv(ni))));
  // }
  // for (p in 1:np) {
  //
    // lt[, p]         ~ normal(mu_lt,         sd_lt);
    // log_lambda[, p] ~ normal(mu_log_lambda, sd_log_lambda);
//
    // log_lambda_raw[, p] ~ normal(0, inv(sqrt(1 - inv(ni))));
  //
  // }
  // log_E_raw ~ normal(0, inv(sqrt(1 - inv(nr))));
  log_a_raw ~ normal(0, inv(sqrt(1 - inv(nr))));
  log_E_raw ~ std_normal();
  b_raw     ~ std_normal();

  if (!prior_only) {
    // likelihood in long form
    for (o in 1:n_observed) {
      int r = idx_rater[o];
      int i = idx_item[o];
      int p = idx_patient[o];

      real location = exp(log_a[r]) * lt[i, p] + b[r];
      real scale    = exp(log_E[r] + log_lambda[i, vary_lambda_across_patients ? p : 1]);

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
      // sds[p, i, r] = exp(log_E[r]);
      // sds[p, i, r] = exp(log_E[r] + log_lambda[i, p]);
      sds[p, i, r] = exp(log_E[r] + log_lambda[i, vary_lambda_across_patients ? p : 1]);
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
