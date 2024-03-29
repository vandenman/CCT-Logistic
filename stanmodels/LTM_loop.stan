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
  real shape_abs = abs(shape);
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

vector raw_to_constrained(vector raw_vec, vector mu, vector log_sd, vector qr, int len, int no_groups, array[] int group_idx, array[] int group_counts) {

  /*
    NOTE: the naive solution

    constrained = mu[group_idx] + exp(log_sd)[group_idx] .* sum_to_zero_QR(raw_vec, qr, len)

    violates sum to zero constrained because the result of

    exp(log_sd)[group_idx] .* sum_to_zero_QR(raw_vec, qr, len)

    need not sum to zero. The naive solution introduces unintended correlations between mu and log_sd.
    This function fixes this behavior.

    Assumed sizes:

    vector[len - 1]   raw_vec
    vector[no_groups] mu
    vector[no_groups] log_sd
    vector[2*len]     qr
    int               len
    int               no_groups
    int               group_idx[no_groups]
    int               group_counts[no_groups]
  */

  vector [len] temp = exp(log_sd)[group_idx] .* sum_to_zero_QR(raw_vec, qr, len);

  // TODO: this requires that the input idx_rater_group is sorted, which is not yet guaranteed!
  // adapted from https://mc-stan.org/docs/2_27/stan-users-guide/ragged-data-structs-section.html
  int pos = 1;
  vector [no_groups] temp_means;
  for (i in 1:no_groups) {
    temp_means[i] = mean(segment(temp, pos, group_counts[i]));
    pos += group_counts[i];
  }

  vector [len] constrained = (mu - temp_means)[group_idx] + temp;

  return constrained;
}

real lpdf_ltm_sliced_observation(array[] int slice, int start, int end,
  array[] int idx_patient,
  array[] int idx_item,
  array[] int idx_rater,
  array[] int idx_time_point,
  array[] int x,
  matrix lt,
  matrix offset_lt,
  vector log_lambda,
  vector log_E,
  vector log_a,
  vector b,
  array[] vector free_thresholds,
  vector threshold_shape,
  vector default_thresholds,
  vector default_threshold_probs,
  int nc,
  int no_time_points,
  int use_free_logistic_thresholds,
  int use_skew_logistic_thresholds
  ) {

  real result = 0.0;
  for (o in start:end) {

    int p = idx_patient[o];
    int i = idx_item[o];
    int r = idx_rater[o];

    real location = lt[p, i];
    // real scale = exp(log_lambda[i] - log_E[r]);
    real scale = exp(log_lambda[i]);
    real E = exp(log_E[r]);

    if (no_time_points == 2) {

      int t = idx_time_point[o];
      // (2 * t - 3) maps {1, 2} to {-1, 1}
      location += (2 * t - 3) * offset_lt[p, i];

    }

    if (use_free_logistic_thresholds) {

      result += ordered_logistic_lpmf(x[o] | location / scale, free_thresholds[r] ./ scale);
      // result += ordered_logistic_lpmf(x[o] | location / (E * scale), free_thresholds[r] ./ scale);

    } else {

      real threshold_scale = exp(log_a[r]);
      real threshold_shift = b[r];

      if (use_skew_logistic_thresholds) {

        result += ordered_logistic_pmf_skew_simplified(
          x[o],
          location, scale, threshold_scale, threshold_shift,
          threshold_shape[r],
          default_threshold_probs,
          nc
        );

      } else {

        // result += ordered_logistic_lpmf(x[o] | location / (E * scale), (threshold_scale .* default_thresholds + threshold_shift) / scale);
        result += ordered_logistic_pmf_simplified(
          x[o],
          location, scale, threshold_scale, threshold_shift,
          default_thresholds,
          nc
        );

      }
    }
  }
  return result;
}


}
data{
  // booleans
  int<lower = 0, upper = 1> prior_only;
  int<lower = 0, upper = 1> debug;
  // int<lower = 0, upper = 1> vectorized;
  // only one of these can be true
  int<lower = 0, upper = 1> use_skew_logistic_thresholds;
  int<lower = 0, upper = 1> use_free_logistic_thresholds;
  int<lower = 0, upper = 1> predict_missings;

  // fit logistic regression or not?
  int<lower = 0, upper = 1> fit_logistic;

  // logistic regression
  int<lower = 0, upper = fit_logistic> store_predictions;

  // sizes
  int<lower = 1> n_observed;
  int<lower = 0> n_missing;
  int<lower = 1> np;
  int<lower = 2> ni;
  int<lower = 2> nr;
  int<lower = 2> nc;
  int<lower = 1> no_rater_groups;
  array[no_rater_groups] int<lower = 1> rater_group_counts;
  int<lower = 1, upper = 2> no_time_points;

  // observed data
  array[n_observed] int<lower = 1, upper = nc> x;

  // indicator vectors
  array[n_observed] int<lower = 1, upper = nr> idx_rater;
  array[n_observed] int<lower = 1, upper = ni> idx_item;
  array[n_observed] int<lower = 1, upper = np> idx_patient;

  array[n_observed] int<lower = 1, upper = no_time_points>  idx_time_point;

  array[nr] int<lower = 1, upper = no_rater_groups> idx_rater_group;

  // logistic regression data
  int<lower = 0>                                                 np_log_reg;
  array[fit_logistic ? np_log_reg : 0] int<lower = 0, upper = 1> log_reg_outcomes;
  int<lower = 0>                                                 no_covariates;
  // this is not space efficient, but required for bernoulli_logit_glm
  matrix[fit_logistic ? np : 0, fit_logistic ? no_covariates : 0] design_matrix_covariates;

  // which observation in log_reg_outcomes corresponds to which patient?
  // basically lt[log_reg_np_idx, ] yields the lt that correspond to log_reg_outcomes
  array[np_log_reg] int<lower=1, upper=np> log_reg_np_idx;

  // missing data
  array[n_missing] int<lower = 1, upper = nr>             idx_rater_missing;
  array[n_missing] int<lower = 1, upper = ni>             idx_item_missing;
  array[n_missing] int<lower = 1, upper = np>             idx_patient_missing;
  array[n_missing] int<lower = 1, upper = no_time_points> idx_time_point_missing;

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

  // array [store_predictions ? np : 0] int ones_np;
  // if (store_predictions) {
  //   for (p in 1:np) {
  //     ones_np[p] = 1;
  //   }
  // }

  // only do this once.
  matrix[fit_logistic ? np_log_reg : 0, fit_logistic ? no_covariates : 0] design_matrix_covariates_observed = design_matrix_covariates[log_reg_np_idx, ];

}
parameters{

  // parameters
  matrix [np - 1, ni] lt_raw;     // item location values
  vector [ni - 1] log_lambda_raw; // item difficulty
  vector [nr - 1] log_E_raw;      // rater competency

  // hyperparameters
  real mu_lt;
  real<lower = 0> sd_lt;
  real<lower = 0> sd_log_lambda;

  vector           [no_rater_groups - 1] mu_log_E_raw;
  vector           [no_rater_groups - 1] log_sd_log_E_raw;

  // conditional parameters
  vector [use_free_logistic_thresholds ? 0 : nr - 1] log_a_raw;      // rater scaling
  vector [use_free_logistic_thresholds ? 0 : nr    ] b;              // rater shifting
  array[nr] ordered[use_free_logistic_thresholds ? nc - 1 : 0] free_thresholds;
  vector [use_skew_logistic_thresholds ? nr     : 0] threshold_shape; // rater shape

  // conditional hyperparameters
  vector           [use_free_logistic_thresholds ? 0 : no_rater_groups] mu_b;
  vector<lower = 0>[use_free_logistic_thresholds ? 0 : no_rater_groups] sd_b;
  vector           [use_free_logistic_thresholds ? 0 : no_rater_groups] log_sd_log_a_raw;

  // logistic regression
  vector[fit_logistic ? 1                                   : 0] log_reg_intercept;
  vector[fit_logistic ? no_covariates + ni * no_time_points : 0] log_reg_slopes;

  // time offsets
  matrix[no_time_points == 2 ? np : 0, no_time_points == 2 ? ni : 0] offset_lt;


}
transformed parameters {

  matrix [np, ni] lt;
  for (i in 1:ni) {
    lt[, i] = sd_lt * sum_to_zero_QR(lt_raw[, i], Q_np, np);
  }
  vector [ni] log_lambda = mu_log_lambda + sd_log_lambda * sum_to_zero_QR(log_lambda_raw, Q_ni, ni);

  vector [no_rater_groups] mu_log_E;
  vector [no_rater_groups] log_sd_log_E;
  vector [nr] log_E;
  if (no_rater_groups > 1) {

    mu_log_E     = sum_to_zero_QR(mu_log_E_raw,     Q_no_rater_groups, no_rater_groups);
    log_sd_log_E = sum_to_zero_QR(log_sd_log_E_raw, Q_no_rater_groups, no_rater_groups);
    log_E = raw_to_constrained(log_E_raw, mu_log_E, log_sd_log_E, Q_nr, nr, no_rater_groups, idx_rater_group, rater_group_counts);

  } else {

    mu_log_E[1]     = 0.0;
    log_sd_log_E[1] = 0.0;
    log_E = sum_to_zero_QR(log_E_raw, Q_nr, nr);

  }

  vector [use_free_logistic_thresholds ? 0 : no_rater_groups] log_sd_log_a;
  vector [use_free_logistic_thresholds ? 0 : nr             ] log_a;
  if (!use_free_logistic_thresholds) {
    if (no_rater_groups > 1) {

      log_sd_log_a = sum_to_zero_QR(log_sd_log_a_raw, Q_no_rater_groups, no_rater_groups);
      log_a = raw_to_constrained(log_a_raw, mu_log_a, log_sd_log_a, Q_nr, nr, no_rater_groups, idx_rater_group, rater_group_counts);

    } else {

      log_sd_log_a[1] = 0.0;
      log_a = sum_to_zero_QR(log_a_raw, Q_nr, nr);

    }
  }

}
model{
  // hyperpriors
  mu_lt         ~ normal(0, 2); // 1 / sqrt(.25)           // mean of item truths
  sd_lt         ~ gamma(a_sd_lt,         b_sd_lt);         // standard deviation of item truths
  mu_b          ~ normal(0, 1);                            // mean of rater shifts
  sd_b          ~ gamma(a_sd_b,          b_sd_b);          // standard deviation of rater shift

  sd_log_lambda ~ gamma(a_sd_log_lambda, b_sd_log_lambda); // standard deviation of item difficulty

  if (no_rater_groups > 1) {
    mu_log_E_raw     ~ normal(0, sd_mu_log_E_raw);
    log_sd_log_E_raw ~ normal(0, sd_mu_log_E_raw);

    if (!use_free_logistic_thresholds)
      log_sd_log_a_raw ~ normal(0, sd_mu_log_E_raw);

  }

  // priors
  to_vector(log_lambda_raw) ~ normal(0, sd_log_lambda_raw);
  to_vector(lt_raw)         ~ normal(0, sd_lt_raw);
  log_E_raw ~ normal(0, sd_log_a_raw); // not a typo

  if (use_free_logistic_thresholds) {

    for (r in 1:nr)
      free_thresholds[r] ~ std_normal();

  } else {

    log_a_raw ~ normal(0, sd_log_a_raw);
    b ~ normal(mu_b[idx_rater_group], sd_b[idx_rater_group]);

  }

  if (use_skew_logistic_thresholds)
    threshold_shape ~ std_normal();

  if (no_time_points == 2)
    to_vector(offset_lt) ~ std_normal();

  // with threads
  target += reduce_sum(
    lpdf_ltm_sliced_observation, idx_patient, 1,
    idx_patient,
    idx_item,
    idx_rater,
    idx_time_point,
    x,
    lt,
    offset_lt,
    log_lambda,
    log_E,
    log_a,
    b,
    free_thresholds,
    threshold_shape,
    default_thresholds,
    default_threshold_probs,
    nc,
    no_time_points,
    use_free_logistic_thresholds,
    use_skew_logistic_thresholds
  );

  // logistic regression
  if (fit_logistic) {
    // priors
    log_reg_intercept ~ std_normal();
    log_reg_slopes    ~ std_normal();

    matrix[np_log_reg, no_covariates + no_time_points * ni] log_reg_design_matrix;
    if (no_time_points == 1) {
      log_reg_design_matrix = append_col(design_matrix_covariates_observed, lt[log_reg_np_idx, ]);
    } else {
      log_reg_design_matrix = append_col(
        design_matrix_covariates_observed,
        append_col(
          lt[log_reg_np_idx, ] - offset_lt[log_reg_np_idx, ],
          lt[log_reg_np_idx, ] + offset_lt[log_reg_np_idx, ]
        )
      );
    }
    log_reg_outcomes ~ bernoulli_logit_glm(log_reg_design_matrix, rep_vector(log_reg_intercept[1], np_log_reg), log_reg_slopes);
  }

  // missing data in IFTE
  // TODO: this should also use reduce_sum!
  // for (o in 1:n_missing) {
  //
  //   int p = idx_patient_missing[o];
  //   int i = idx_item_missing[o];
  //   int r = idx_rater_missing[o];
  //
  //   real location = lt[p, i];
  //   real scale = exp(log_lambda[i] - log_E[r]);
  //
  //   if (no_time_points == 2) {
  //
  //     int t = idx_time_point[o];
  //     // (2 * t - 3) maps {1, 2} to {-1, 1}
  //     location += (2 * t - 3) * offset_lt[p, i];
  //
  //   }
  //
  //   int loop_limit = nc;//i == 1 ? 13 : nc;
  //   vector[loop_limit] missing_log_probs;
  //   if (use_free_logistic_thresholds) {
  //
  //     for (c in 1:loop_limit) {
  //       missing_log_probs[c] = ordered_logistic_lpmf(c | location / scale, free_thresholds[r] ./ scale);
  //     }
  //
  //   } else {
  //
  //     real threshold_scale = exp(log_a[r]);
  //     real threshold_shift = b[r];
  //
  //     if (use_skew_logistic_thresholds) {
  //
  //       for (c in 1:loop_limit) {
  //         missing_log_probs[c] = ordered_logistic_pmf_skew_simplified(
  //           c,
  //           location, scale, threshold_scale, threshold_shift,
  //           threshold_shape[r],
  //           default_threshold_probs,
  //           nc
  //         );
  //       }
  //
  //     } else {
  //
  //       for (c in 1:loop_limit) {
  //         missing_log_probs[c] = ordered_logistic_pmf_simplified(
  //           c,
  //           location, scale, threshold_scale, threshold_shift,
  //           default_thresholds,
  //           nc
  //         );
  //       }
  //     }
  //   }
  //   print("log_sum_exp(missing_log_probs)");
  //   print(log_sum_exp(missing_log_probs));
  //   target += log_sum_exp(missing_log_probs); // TODO: should this be 2 * logProbs?
  // }


}
generated quantities {

  matrix<upper = 0>[debug ? nc : 0, debug ? n_observed : 0] log_probs;

  // this is incredibly memory intense!
  // TODO: it would be MUCH cleaner to somehow reuse lpdf_ltm_sliced_observation here!
  if (debug) {

    // likelihood in long form
    for (o in 1:n_observed) {

      int p = idx_patient[o];
      int i = idx_item[o];
      int r = idx_rater[o];

      real location = lt[p, i];
      real scale = exp(log_lambda[i] - log_E[r]);

      if (no_time_points == 2) {

        int t = idx_time_point[o];
        // (2 * t - 3) maps {1, 2} to {-1, 1}
        location += (2 * t - 3) * offset_lt[p, i];

      }

      if (use_free_logistic_thresholds) {

        for (c in 1:nc)
          log_probs[c, o] = ordered_logistic_lpmf(c | location / scale, free_thresholds[r] ./ scale);

      } else {

        real threshold_scale = exp(log_a[r]);
        real threshold_shift = b[r];

        if (use_skew_logistic_thresholds) {
          for (c in 1:nc) {
            log_probs[c, o] = ordered_logistic_pmf_skew_simplified(
              c, location, scale, threshold_scale, threshold_shift,
              threshold_shape[r],
              default_threshold_probs,
              nc
            );
          }
        } else {
          for (c in 1:nc) {
            log_probs[c, o] = ordered_logistic_pmf_simplified(
              c, location, scale, threshold_scale, threshold_shift,
              default_thresholds,
              nc
            );
          }
        }
      }
    }
  }

  vector<lower = 0, upper = 1>[fit_logistic * store_predictions ? np : 0] log_reg_predictions;
  if (fit_logistic * store_predictions) {

    // NOTE: the code below makes predictions for both the observed log_reg_outcomes and the unobserved log_reg_outcomes.
    for (p in 1:np) {
      vector[1] alpha;
      if (no_time_points == 1) {
        alpha = log_reg_intercept +
          dot_product(design_matrix_covariates[p, ], log_reg_slopes[1:no_covariates]) +
          dot_product(lt[p, ], log_reg_slopes[(no_covariates + 1):(no_covariates + ni)]);
      } else {
        alpha = log_reg_intercept +
          dot_product(design_matrix_covariates[p, ], log_reg_slopes[1:no_covariates]) +
          dot_product(lt[p, ] - offset_lt[p, ], log_reg_slopes[(no_covariates      + 1):(no_covariates +   ni)]) +
          dot_product(lt[p, ] + offset_lt[p, ], log_reg_slopes[(no_covariates + ni + 1):(no_covariates + 2*ni)]);
      }
      log_reg_predictions[p] = exp(bernoulli_logit_lpmf(1 | alpha));
      // bernoulli_logit_glm_lpmf(ones_np | append_col(design_matrix_covariates, lt), log_reg_intercept, log_reg_slopes);
    }
  }

  // array[fit_logistic * store_slopes ? np : 0, fit_logistic * store_slopes ? ni : 0, fit_logistic * store_slopes ? no_time_points : 0] real log_reg_derived_slopes;
  // if (fit_logistic * store_slopes) {
  //   for (p in 1:np) {
  //     if (no_time_points == 1) {
  //       log_reg_derived_slopes[p, , 1] = lt[p, ] .* log_reg_slopes[(no_covariates + 1):(no_covariates + ni)]
  //     } else {
  //       log_reg_derived_slopes[p, , 1] = (lt[p, ] - offset_lt[p, ]) .* log_reg_slopes[(no_covariates      + 1):(no_covariates +   ni)]
  //       log_reg_derived_slopes[p, , 2] = (lt[p, ] - offset_lt[p, ]) .* log_reg_slopes[(no_covariates + ni + 1):(no_covariates + 2*ni)]
  //     }
  //   }
  // }

  vector<lower = 0, upper = nc>[predict_missings && n_missing > 0 ? n_missing : 0] predicted_missings;
  if (predict_missings && n_missing > 0) {
    for (o in 1:n_missing) {

      int p = idx_patient_missing[o];
      int i = idx_item_missing[o];
      int r = idx_rater_missing[o];

      real location = lt[p, i];
      real scale = exp(log_lambda[i] - log_E[r]);

      if (no_time_points == 2) {

        int t = idx_time_point_missing[o];
        // (2 * t - 3) maps {1, 2} to {-1, 1}
        location += (2 * t - 3) * offset_lt[p, i];

      }

      int loop_limit = nc;//i == 1 ? 13 : nc;
      vector[loop_limit] missing_log_probs;
      if (use_free_logistic_thresholds) {

        for (c in 1:loop_limit) {
          missing_log_probs[c] = ordered_logistic_lpmf(c | location / scale, free_thresholds[r] ./ scale);
        }

      } else {

        real threshold_scale = exp(log_a[r]);
        real threshold_shift = b[r];

        if (use_skew_logistic_thresholds) {

          for (c in 1:loop_limit) {
            missing_log_probs[c] = ordered_logistic_pmf_skew_simplified(
              c,
              location, scale, threshold_scale, threshold_shift,
              threshold_shape[r],
              default_threshold_probs,
              nc
            );
          }

        } else {

          for (c in 1:loop_limit) {
            missing_log_probs[c] = ordered_logistic_pmf_simplified(
              c,
              location, scale, threshold_scale, threshold_shift,
              default_thresholds,
              nc
            );
          }
        }
      }
      predicted_missings[o] = categorical_logit_rng(missing_log_probs);
    }
  }

}
