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

real stan_log_inv_logit_diff(real x, real y) {
  // the manual defines log_inv_logit_diff but it's not accepted?
  return x - log1p_exp(x) + log1m_exp(y - x) - log1p_exp(y);
}
real stan_ordered_logistic_pmf(int x, real location, real scale, vector delta) {
  return ordered_logistic_lpmf(x | location / scale, delta ./ scale);
}
real ordered_logistic_pmf_simplified(int x, real location, real scale, real threshold_scale, real threshold_shift, vector default_thresholds, int nc) {
  if (x == 1) {
    return log1m_inv_logit((location - threshold_shift - threshold_scale * default_thresholds[x    ]) / scale);
  } else if (x == nc) {
    return log_inv_logit((location - threshold_shift - threshold_scale * default_thresholds[x - 1]) / scale);
  } else {
    // return stan_log_inv_logit_diff(
    //   (location - threshold_shift - threshold_scale * default_thresholds[x - 1]) / scale,
    //   (location - threshold_shift - threshold_scale * default_thresholds[x    ]) / scale
    // );
    return log_inv_logit_diff(
      (location - threshold_shift - threshold_scale * default_thresholds[x - 1]) / scale,
      (location - threshold_shift - threshold_scale * default_thresholds[x    ]) / scale
    );
  }
}

int signum(real x) {
  return x < 0 ? -1 : x > 0;
}

real sta_qELGW(real p, real location, real scale, real shape) {

  real direction = signum(shape);
  real shape_abs = fabs(shape);
  return location + direction * scale * log((1 - shape_abs * log(0.5 - direction * (p - 0.5))) ^ (1.0 / shape_abs) - 1);

}

real ordered_logistic_pmf_skew_simplified(int x, real location, real scale, real threshold_scale, real threshold_shift, real threshold_shape, vector default_threshold_probs, int nc) {
  if (x == 1) {
    return log1m_inv_logit((location - sta_qELGW(default_threshold_probs[1], threshold_shift, threshold_scale, threshold_shape)) / scale);
  } else if (x == nc) {
    return log_inv_logit((location - sta_qELGW(default_threshold_probs[nc - 1], threshold_shift, threshold_scale, threshold_shape)) / scale);
  } else {
    return log_inv_logit_diff(
      (location - sta_qELGW(default_threshold_probs[x - 1], threshold_shift, threshold_scale, threshold_shape)) / scale,
      (location - sta_qELGW(default_threshold_probs[x    ], threshold_shift, threshold_scale, threshold_shape)) / scale
    );
  }
}


}
data{
  // booleans
  int<lower = 0, upper = 1> prior_only;
  int<lower = 0, upper = 1> debug;
  int<lower = 0, upper = 1> vary_lambda_across_patients;
  int<lower = 0, upper = 1> vectorized;
  int<lower = 0, upper = 1> use_skew_logistic_thresholds;

  // sizes
  int<lower = 1> n_observed;
  int<lower = 0> n_missing;
  int<lower = 1> np;
  int<lower = 2> ni;
  int<lower = 2> nr;
  int<lower = 2> nc;
  int<lower = 1> no_rater_groups;
  int<lower = 1> rater_group_counts[no_rater_groups];

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

  vector[2*np] Q_np = Q_sum_to_zero_QR(np);
  real sd_lt_raw = inv_sqrt(1 - inv(np));

  vector[2*nr] Q_nr = Q_sum_to_zero_QR(nr);
  real sd_log_E_raw = inv_sqrt(1 - inv(nr));

  vector[2*no_rater_groups] Q_no_rater_groups = Q_sum_to_zero_QR(no_rater_groups);
  real sd_mu_log_E_raw = no_rater_groups > 1 ? 1.0 : inv_sqrt(1 - inv(no_rater_groups));

  vector[nc-1] default_threshold_probs;
  vector[nc-1] default_thresholds; // <- 1:nc;
  int nt = nc - 1; // no thresholds

  for (t in 1:nt) {
    real tt = t; // Stan can only cast in this way and we need to cast to real to avoid integer division
    default_thresholds[t] = log(tt / (nc - tt));
    default_threshold_probs[t] = tt / nc;
  }

}
parameters{

  // parameters
  matrix [np - 1, ni] lt_raw;         // item location values
  // TODO: should the layout be flipped (like lt_raw)?
  matrix [ni - 1, vary_lambda_across_patients ? np : 1] log_lambda_raw; // item difficulty
  vector [nr - 1] log_E_raw;          // rater competency

  ordered[nc - 1] free_thresholds[nr];

  // hyperparameters
  real mu_lt;
  real<lower = 0> sd_lt;

  real<lower = 0>              sd_log_lambda;

  vector<lower = 0>[no_rater_groups] sd_log_E;

  vector           [no_rater_groups - 1] mu_log_E_raw;
  vector           [no_rater_groups - 1] log_sd_log_E_raw;

}
transformed parameters {

  // identification constraints
  matrix [ni, vary_lambda_across_patients ? np : 1] log_lambda;
  matrix [np, ni] lt;
  vector [no_rater_groups] mu_log_E;
  vector [no_rater_groups] log_sd_log_E;
  vector [nr] log_E;
  if (no_rater_groups > 1) {
    mu_log_E     = sum_to_zero_QR(mu_log_E_raw,     Q_no_rater_groups, no_rater_groups);
    log_sd_log_E = sum_to_zero_QR(log_sd_log_E_raw, Q_no_rater_groups, no_rater_groups);

    // compute a correction so that the means are still zero after multiplying with the within group variances
    vector [nr] temp = exp(log_sd_log_E)[idx_rater_group] .* sum_to_zero_QR(log_E_raw, Q_nr, nr);

    // TODO: this requires that idx_rater_group is sorted, which is not yet guaranteed!
    // adapted from https://mc-stan.org/docs/2_27/stan-users-guide/ragged-data-structs-section.html
    int pos = 1;
    vector [no_rater_groups] temp_means;
    for (i in 1:no_rater_groups) {
      temp_means[i] = mean(segment(temp, pos, rater_group_counts[i]));
      pos += rater_group_counts[i];
    }
    log_E = (mu_log_E - temp_means)[idx_rater_group] + temp;

    // vector [no_rater_groups] temp_sums = rep_vector(0, no_rater_groups);
    // for (i in 1:nr) {
    //   temp_sums[idx_rater_group[i]] += temp[i];
    // }
    // log_E = (mu_log_E - temp_sums ./ rater_group_counts)[idx_rater_group] + temp;
  } else {
    log_E = sum_to_zero_QR(log_E_raw, Q_nr, nr);
  }

  for (i in 1:ni) {
    lt[, i] = sd_lt * sum_to_zero_QR(lt_raw[, i], Q_np, np);
  }
  if (vary_lambda_across_patients) {
    for (p in 1:np) {
      log_lambda[, p] = mu_log_lambda + sd_log_lambda * sum_to_zero_QR(log_lambda_raw[, p], Q_ni, ni);
    }
  } else {
    log_lambda[, 1] = mu_log_lambda + sd_log_lambda * sum_to_zero_QR(log_lambda_raw[, 1], Q_ni, ni);
  }
}
model{
  // hyperpriors
  mu_lt         ~ normal(0, 2); // 1 / sqrt(.25)           // mean of item truths
  // mu_log_E      ~ normal(0, 10);// 1 / sqrt(0.01)          // mean of rater competence
  sd_lt         ~ gamma(a_sd_lt,         b_sd_lt);         // standard deviation of item truths
  sd_log_lambda ~ gamma(a_sd_log_lambda, b_sd_log_lambda); // standard deviation of item difficulty
  // sd_log_E      ~ gamma(a_sd_log_E,      b_sd_log_E);      // standard deviation of rater competence

  // priors
  to_vector(log_lambda_raw) ~ normal(0, sd_log_lambda_raw);
  to_vector(lt_raw)         ~ normal(0, sd_lt_raw);
  log_E_raw ~ normal(0, sd_log_E_raw);

  if (no_rater_groups > 1) {
    mu_log_E_raw     ~ normal(0, sd_mu_log_E_raw);
    log_sd_log_E_raw ~ normal(0, sd_mu_log_E_raw);
  }

  // log_E ~ normal(mu_log_E[idx_rater_group], sd_log_E[idx_rater_group]);

  // should force thresholds for unobserved categories towards 0
  for (r in 1:nr)
    free_thresholds[r] ~ std_normal();

  // print("free_thresholds    = ", free_thresholds);
  // print("free_thresholds[1] = ", free_thresholds[1]);
  // print("free_thresholds[2] = ", free_thresholds[2]);
  // likelihood in long form
  for (o in 1:n_observed) {

    int p = idx_patient[o];
    int i = idx_item[o];
    int r = idx_rater[o];

    real location = lt[p, i];
    // probably doesn't work to add log_E because it's unidentified when free_thresholds becomes twice as wide
    // real scale = exp(log_lambda[i, vary_lambda_across_patients ? p : 1]);// - log_E[r]);
    real scale = exp(log_lambda[i, vary_lambda_across_patients ? p : 1] - log_E[r]);

    target += ordered_logistic_lpmf(x[o] | location / scale, free_thresholds[r] ./ scale);

    // target += ordered_logistic_pmf_simplified(
    //   x[o],
    //   location, scale, 1.0, 0.0,
    //   free_thresholds[r],
    //   nc
    // );

  }
}
generated quantities {

  matrix[debug ? nc : 0, debug ? n_observed : 0] log_probs;

  if (debug) {
    // likelihood in long form
    for (o in 1:n_observed) {

      int p = idx_patient[o];
      int i = idx_item[o];
      int r = idx_rater[o];

      real location = lt[p, i];
      // probably doesn't work to add log_E because it's unidentified when free_thresholds becomes twice as wide
      real scale = exp(log_lambda[i, vary_lambda_across_patients ? p : 1] - log_E[r]);

      for (c in 1:nc) {
        log_probs[c, o] = ordered_logistic_lpmf(c | location / scale, free_thresholds[r] ./ scale);
        // log_probs[c, o] = ordered_logistic_pmf_simplified(
        //   x[o],
        //   location, scale, 1.0, 0.0,
        //   free_thresholds[r],
        //   nc
        // );
      }
    }
  }
}
