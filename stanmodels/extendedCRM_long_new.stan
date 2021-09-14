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

  // https://discourse.mc-stan.org/t/test-soft-vs-hard-sum-to-zero-constrain-choosing-the-right-prior-for-soft-constrain/3884/31
  vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  vector sum_to_zero_QR(vector x_raw, vector Q_r, int N) {
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i]  = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
}
data{
  // booleans
  int<lower = 0, upper = 1> prior_only;
  int<lower = 0, upper = 1> debug;
  int<lower = 0, upper = 1> vary_lambda_across_patients;
  int<lower = 0, upper = 1> vectorized;

  // sizes
  int<lower = 1> n_observed;
  int<lower = 0> n_missing;
  int<lower = 2> nr;
  int<lower = 2> ni;
  int<lower = 1> np;
  int<lower = 1> no_rater_groups;

  // observed data
  real x[n_observed];
  // indicator vectors
  int<lower = 1, upper = nr> idx_rater  [n_observed];
  int<lower = 1, upper = ni> idx_item   [n_observed];
  int<lower = 1, upper = np> idx_patient[n_observed];

  int<lower = 1, upper = no_rater_groups> idx_rater_group[nr];

  // faster than rep_each
  // int<lower = 1, upper = ni * np> idx_longer_lt[vectorized ? nr * ni * np : 0];
  int<lower = 1, upper = ni * np> idx_longer_lt[nr * ni * np];

  // missing data
  int<lower = 1, upper = nr> idx_rater_missing  [n_missing];
  int<lower = 1, upper = ni> idx_item_missing   [n_missing];
  int<lower = 1, upper = np> idx_patient_missing[n_missing];

  // constant hyperparameters, see Anders, Oravecz, and Batchelder (2014) p. 5
  real mu_log_lambda;
  vector[no_rater_groups] mu_log_a;
  // real mu_b;

  // constants governing the distributions over hyperparameters
  real a_sd_lt;
  real b_sd_lt;
  real a_sd_log_lambda;
  real b_sd_log_lambda;
  real a_sd_log_E;
  real b_sd_log_E;
  real a_sd_log_a;
  real b_sd_log_a;
  real a_sd_b;
  real b_sd_b;

}
transformed data{

  vector[2*ni] Q_ni = Q_sum_to_zero_QR(ni);
  real sd_log_lambda_raw = inv_sqrt(1 - inv(ni));

  vector[2*nr] Q_nr = Q_sum_to_zero_QR(nr);
  real sd_log_a_raw = inv_sqrt(1 - inv(nr));

  // matrix[nr, nr - 1] QR_nr = getQR(nr);
  // matrix[ni, ni - 1] QR_ni = getQR(ni);

}
parameters{
  // See Anders, Oravecz, and Batchelder (2014) p. 5

  // hyperparameters
  real mu_lt;
  real<lower = 0> sd_lt;

  real<lower = 0> sd_log_lambda;

  // vector           [no_rater_groups] mu_log_a;
  real<lower = 0> sd_log_a;
  // vector<lower = 0>[no_rater_groups] sd_log_a;
  vector           [no_rater_groups] mu_b;
  vector<lower = 0>[no_rater_groups] sd_b;
  vector           [no_rater_groups] mu_log_E;
  vector<lower = 0>[no_rater_groups] sd_log_E;

  // parameters
  matrix [np, ni] lt;         // item location values
  // matrix [ni, vary_lambda_across_patients ? np : 1] log_lambda; // item difficulty
  vector [ni - 1] log_lambda_raw; // item difficulty
  // vector [ni] log_lambda; // item difficulty

  vector [nr - 1] log_a_raw; // rater scaling
  // vector [nr] log_a; // rater scaling

  // vector [nr] log_E_raw; // rater competency
  // vector [nr]     b_raw; // rater shifting
  vector [nr] log_E; // rater competency
  vector [nr]     b;     // rater shifting

}
transformed parameters{
  // vector[ni] log_lambda = mu_log_lambda + sd_log_lambda * (QR_ni * log_lambda_raw);
  vector [ni] log_lambda = mu_log_lambda             + sd_log_lambda              * sum_to_zero_QR(log_lambda_raw, Q_ni, ni);
  // vector [nr] log_a      = mu_log_a[idx_rater_group] + sd_log_a[idx_rater_group] .* sum_to_zero_QR(log_a_raw, Q_nr, nr);
  vector [nr] log_a      = mu_log_a[idx_rater_group] + sd_log_a * sum_to_zero_QR(log_a_raw, Q_nr, nr);
  // vector [nr] log_E      = mu_log_E[idx_rater_group] + sd_log_E[idx_rater_group] .*                log_E_raw;
  // vector [nr] b          = mu_b    [idx_rater_group] + sd_b    [idx_rater_group] .*                    b_raw;
}
model{
profile("priors") {
  // hyperpriors
  mu_lt         ~ normal(0, 2); // 1 / sqrt(.25)
  // mu_log_a      ~ normal(0, 1);// 1 / sqrt(0.01)
  mu_b          ~ normal(0, 1);// 1 / sqrt(0.01)
  mu_log_E      ~ normal(0, 10);// 1 / sqrt(0.01)
  sd_lt         ~ gamma(a_sd_lt,         b_sd_lt);
  sd_log_lambda ~ gamma(a_sd_log_lambda, b_sd_log_lambda);
  sd_log_E      ~ gamma(a_sd_log_E,      b_sd_log_E);
  sd_log_a      ~ gamma(a_sd_log_a,      b_sd_log_a);
  sd_b          ~ gamma(a_sd_b,          b_sd_b);

  // priors
  log_lambda_raw ~ normal(0, sd_log_lambda_raw);
  log_a_raw      ~ normal(0, sd_log_a_raw);
  // log_lambda     ~ normal(mu_log_lambda, sd_log_lambda);
  to_vector(lt)  ~ normal(mu_lt, sd_lt);
  // log_a          ~ normal(mu_log_a, sd_log_a);
  // log_E          ~ normal(mu_log_E, sd_log_E);
  // b              ~ normal(mu_b,     sd_b);

  log_E          ~ normal(mu_log_E[idx_rater_group], sd_log_E[idx_rater_group]);
  b              ~ normal(mu_b[idx_rater_group],     sd_b[idx_rater_group]);
}
  // log_E_raw ~ std_normal();
  // b_raw     ~ std_normal();

profile("likelihood") {
  if (!prior_only) {

    if (vectorized) {

      vector[n_observed] locations = exp(log_a)[idx_rater] .* to_vector(lt)[idx_longer_lt] + b[idx_rater];
      vector[n_observed] scales    = exp(log_E[idx_rater] + log_lambda[idx_item]);
      target += normal_lpdf(x | locations, scales);

    } else {

      for (o in 1:n_observed) {
        int r = idx_rater[o];
        int i = idx_item[o];
        int p = idx_patient[o];

        real location = exp(log_a[r]) * lt[p, i] + b[r];
        real scale    = exp(log_E[r] + log_lambda[i]);

        x[o] ~ normal(location, scale);
      }
    }
  }
}
}
generated quantities {
profile("generated quantities") {
  // 3d arrays don't work when not all combinations are observed
  vector[debug ? n_observed : 0] locations;
  vector[debug ? n_observed : 0] scales;

  if (debug) {

    if (vectorized) {

      locations = exp(log_a)[idx_rater] .* to_vector(lt)[idx_longer_lt] + b[idx_rater];
      scales = exp(log_E[idx_rater] + log_lambda[idx_item]);

    } else {
      for (o in 1:n_observed) {
        int r = idx_rater[o];
        int i = idx_item[o];
        int p = idx_patient[o];

        locations[o] = exp(log_a[r]) * lt[p, i] + b[r];
        scales[o] = exp(log_E[r] + log_lambda[i]);
      }
    }
  }
}
}
