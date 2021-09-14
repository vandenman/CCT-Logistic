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
  vector rep_each(vector x, int K) {
    // from https://stackoverflow.com/a/60834449
    int N = rows(x);
    vector[N * K] y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
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

  // matrix[nr, nr - 1] QR_nr = getQR(nr);
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
  matrix [np, ni] lt;         // item location values
  // matrix [ni, vary_lambda_across_patients ? np : 1] log_lambda; // item difficulty
  vector [ni-1] log_lambda_raw; // item difficulty


  vector [nr] log_a; // rater scaling
  vector [nr] log_E; // rater competency
  vector [nr] b;     // rater shifting

}
transformed parameters{
  vector[ni] log_lambda = mu_log_lambda + sd_log_lambda * QR_ni * log_lambda_raw;
}
model{

  // hyperpriors
  mu_lt         ~ normal(0, 0.25);
  mu_log_E      ~ normal(0, .01);
  sd_lt         ~ gamma(1, 1);
  sd_log_lambda ~ gamma(1, 1);
  sd_log_E      ~ gamma(1, 1);
  sd_log_a      ~ gamma(1, 1);
  sd_b          ~ gamma(1, 1);

  // priors
  log_lambda_raw ~ normal(0, inv(sqrt(1 - inv(ni))));
  // to_vector(log_lambda) ~ normal(mu_log_lambda, sd_log_lambda);
  to_vector(lt)         ~ normal(mu_lt,    sd_lt);
  log_a                 ~ normal(mu_log_a, sd_log_a);
  log_E                 ~ normal(mu_log_E, sd_log_E);
  b                     ~ normal(mu_b,     sd_b);

  if (!prior_only) {
    // likelihood in long form
    vector[n_observed] locations = exp(log_a)[idx_rater] .* rep_each(to_vector(lt), nr) + b[idx_rater];
    vector[n_observed] scales    = exp(log_E[idx_rater] + log_lambda[idx_item]);
    x ~ normal(locations, scales);

    // for (o in 1:n_observed) {
    //   int r = idx_rater[o];
    //   int i = idx_item[o];
    //   int p = idx_patient[o];
    //
    //   real location = exp(log_a[r]) * lt[i, p] + b[r];
    //   // real scale    = exp(2.0 * log_a[r] + log_E[r] + log_lambda[i, p]);
    //   // real scale    = exp(log_E[r]);
    //   // real scale    = exp(log_E[r] + log_lambda[i, p]);
    //   real scale    = exp(log_E[r] + log_lambda[i]);
    //
    //   x[o] ~ normal(location, scale);
    // }
  }
}
generated quantities {

  // these 3d arrays don't work when not all combinations are observed
  // real             mus[debug ? np : 0, debug ? ni : 0, debug ? nr : 0];
  // real<lower = 0>  sds[debug ? np : 0, debug ? ni : 0, debug ? nr : 0];
  //
  // if (debug) {
    // for (o in 1:n_observed) {
    //   int r = idx_rater[o];
    //   int i = idx_item[o];
    //   int p = idx_patient[o];
    //
    //   mus[p, i, r] = exp(log_a[r]) * lt[i, p] + b[r];
    //   // sds[p, i, r] = exp(log_E[r]);
    //   // sds[p, i, r] = exp(log_E[r] + log_lambda[i, p]);
    //   sds[p, i, r] = exp(log_E[r] + log_lambda[i]);
    // }
  // }

  vector[debug ? n_observed : 0] mus;
  vector[debug ? n_observed : 0] sds;

  if (debug) {

    mus = exp(log_a[idx_rater]) .* rep_each(to_vector(lt), nr) + b[idx_rater];
    sds = exp(log_E[idx_rater] + log_lambda[idx_item]);

  }

}
