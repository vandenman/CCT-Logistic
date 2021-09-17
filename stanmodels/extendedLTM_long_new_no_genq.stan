functions {
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
  int<lower = 1> np;
  int<lower = 2> ni;
  int<lower = 2> nr;
  int<lower = 2> nc;
  int<lower = 1> no_rater_groups;

  // observed data
  int x[n_observed];

    // indicator vectors
  int<lower = 1, upper = nr> idx_rater  [n_observed];
  int<lower = 1, upper = ni> idx_item   [n_observed];
  int<lower = 1, upper = np> idx_patient[n_observed];

  int<lower = 1, upper = no_rater_groups> idx_rater_group[nr];

  // missing data
  int<lower = 1, upper = nr> idx_rater_missing  [n_missing];
  int<lower = 1, upper = ni> idx_item_missing   [n_missing];

  // constant hyperparameters, see Anders, Oravecz, and Batchelder (2014) p. 5
  real mu_log_lambda;
  vector[no_rater_groups] mu_log_a;

  // constants governing the distributions over hyperparameters
  real<lower = 0> a_sd_lt;
  real<lower = 0> b_sd_lt;
  real<lower = 0> a_sd_log_lambda;
  real<lower = 0> b_sd_log_lambda;
  real<lower = 0> a_sd_log_E;
  real<lower = 0> b_sd_log_E;
  real<lower = 0> a_sd_log_a;
  real<lower = 0> b_sd_log_a;
  real<lower = 0> a_sd_b;
  real<lower = 0> b_sd_b;

}
transformed data{

  vector[2*ni] Q_ni = Q_sum_to_zero_QR(ni);
  real sd_log_lambda_raw = inv_sqrt(1 - inv(ni));

  vector[2*nr] Q_nr = Q_sum_to_zero_QR(nr);
  real sd_log_a_raw = inv_sqrt(1 - inv(nr));

  vector[2*np] Q_np = Q_sum_to_zero_QR(np);
  real sd_lt_raw = inv_sqrt(1 - inv(np));

  vector[nc-1] lambda; // <- 1:nc;
  int nt; // no thresholds
  nt = nc - 1;

  for (t in 1:nt) {
    real tt = t; // Stan can only cast in this way and we need to cast to real to avoid integer division
    lambda[t] = log(tt / (nc - tt));
  }


}
parameters{

  // parameters
  // matrix [np, ni] lt;          // item location values
  matrix [np - 1, ni] lt_raw;         // item location values
  vector [ni - 1] log_lambda_raw; // item difficulty
  vector [nr - 1] log_a_raw;      // rater scaling
  vector [nr]     log_E;          // rater competency
  vector [nr]     b;              // rater shifting
  // vector [nr - 1] b_raw;       // rater shifting

  // hyperparameters
  real mu_lt;
  real<lower = 0> sd_lt;

  // TODO: log_a should also vary across raters!
  vector           [no_rater_groups] mu_b;
  vector<lower = 0>[no_rater_groups] sd_b;
  vector           [no_rater_groups] mu_log_E;
  vector<lower = 0>[no_rater_groups] sd_log_E;

  real<lower = 0>              sd_log_a;
  real<lower = 0>              sd_beta;
  real<lower = 0>              sd_kappa;
  real<lower = 0>              sd_log_lambda;

}
transformed parameters {

  // identification constraints
  vector [ni] log_lambda = mu_log_lambda             + sd_log_lambda          * sum_to_zero_QR(log_lambda_raw, Q_ni, ni);
  vector [nr] log_a      = mu_log_a[idx_rater_group] + sd_log_a               * sum_to_zero_QR(log_a_raw, Q_nr, nr);
  // vector [nr] b          = mu_b    [idx_rater_group] + sd_b[idx_rater_group] .* sum_to_zero_QR(b_raw,     Q_nr, nr);
  matrix [np, ni] lt;
  for (i in 1:ni) {
    lt[, i] = sd_lt * sum_to_zero_QR(lt_raw[, i], Q_np, np);
  }


}
model{
  // hyperpriors
  mu_lt         ~ normal(0, 2); // 1 / sqrt(.25)           // mean of item truths
  mu_b          ~ normal(0, 1);                            // mean of rater shifts
  mu_log_E      ~ normal(0, 10);// 1 / sqrt(0.01)          // mean of rater competence
  sd_lt         ~ gamma(a_sd_lt,         b_sd_lt);         // standard deviation of item truths
  sd_log_lambda ~ gamma(a_sd_log_lambda, b_sd_log_lambda); // standard deviation of item difficulty
  sd_b          ~ gamma(a_sd_b,          b_sd_b);          // standard deviation of rater shift
  sd_log_a      ~ gamma(a_sd_log_a,      b_sd_log_a);      // standard deviation of rater scale
  sd_log_E      ~ gamma(a_sd_log_E,      b_sd_log_E);      // standard deviation of rater competence

  // priors
  log_lambda_raw ~ normal(0, sd_log_lambda_raw);
  log_a_raw      ~ normal(0, sd_log_a_raw);
  // to_vector(lt)  ~ normal(mu_lt, sd_lt);
  to_vector(lt_raw)  ~ std_normal();
  log_E          ~ normal(mu_log_E[idx_rater_group], sd_log_E[idx_rater_group]);

  // b_raw          ~ normal(0, sd_log_a_raw); // sd_log_a_raw is not a mistake
  b              ~ normal(mu_b    [idx_rater_group],     sd_b[idx_rater_group]);

  // likelihood in long form
  for (o in 1:n_observed) {

    int p = idx_patient[o];
    int i = idx_item[o];
    int r = idx_rater[o];

    real scale = exp(log_lambda[i] - log_E[r]);
    real location = lt[p, i];
    vector[nt] delta = exp(log_a[r]) * lambda + b[r];

    x[o] ~ ordered_logistic(location / scale, delta ./ scale);

  }

  // // marginalize over missing values
  // for (m in 1:nMissing) {
  //   vector[nc] logProbs;
  //   int r = idxRaterMissing[m];
  //   int i = idxItemMissing[m];
  //   real scale = kappa[i] / zeta[r];
  //   real location = lt[i];
  //   vector[nt] delta = alpha[r] * lambda + beta[r];
  //   for (c in 1:nc) {
  //     logProbs[c] = ordered_logistic_lpmf(c | location / scale, delta ./ scale);
  //   }
  //   target += log_sum_exp(2*logProbs);
  // }
}
