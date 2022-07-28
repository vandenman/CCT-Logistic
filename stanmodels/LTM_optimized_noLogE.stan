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

  // TODO: begin of stuff to delete
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

real lpdf_ltm_sliced_observation(array[] int slice, int start, int end,
  array[] int idx_patient,
  array[] int idx_item,
  array[] int idx_rater,
  array[] int idx_time_point,
  array[] int x,
  matrix lt,
  matrix offset_lt,
  matrix log_lambda,
  vector log_E,
  vector log_a,
  vector b,
  array[] vector free_thresholds,
  vector threshold_shape,
  vector default_thresholds,
  vector default_threshold_probs,
  int nc,
  int no_time_points,
  int vary_lambda_across_patients,
  int use_free_logistic_thresholds,
  int use_skew_logistic_thresholds
  ) {

  real result = 0.0;
  for (o in start:end) {

    int p = idx_patient[o];
    int i = idx_item[o];
    int r = idx_rater[o];

    real location = lt[p, i];
    real scale = exp(log_lambda[i, vary_lambda_across_patients ? p : 1] - log_E[r]);

    if (no_time_points == 2) {

      int t = idx_time_point[o];
      // (2 * t - 3) maps {1, 2} to {-1, 1}
      location += (2 * t - 3) * offset_lt[p, i];

    }

    if (use_free_logistic_thresholds) {

      result += ordered_logistic_lpmf(x[o] | location / scale, free_thresholds[r] ./ scale);

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

// TODO: end of stuff to delete

vector sign_vec(vector x) {
  vector [rows(x)] y;
  for (i in 1:rows(x)) {
    y[i] = (x[i] > 0 ? 1.0 : -1.0);
  }
  return y;
}

vector fast_log_lik_common(
  matrix lt, vector log_lambda, vector log_E, vector delta0, vector delta1,
  matrix offset_lt, vector sign_offset_lt,
  int nobs, int no_time_points,
  array[] int idx_item_patient_vectorized,
  array[] int idx_item,
  array[] int idx_rater) {

  vector[nobs] locations = to_vector(lt)[idx_item_patient_vectorized];
  vector[nobs] scales    = exp(log_lambda[idx_item] - log_E[idx_rater]);

  if (no_time_points == 2) {
    locations += sign_offset_lt .* to_vector(offset_lt)[idx_item_patient_vectorized];
  }

  vector [nobs] d0 = (locations - delta0) ./ scales;
  vector [nobs] d1 = (locations - delta1) ./ scales;

  // could be simplified to
  // log_inv_logit_diff(d0, d1)
  // but see https://github.com/stan-dev/math/issues/2796

  vector [nobs] v0 = inv_logit(d0);
  vector [nobs] v1 = inv_logit(d1);

  print("locations");
  print(locations);
  print("scales");
  print(scales);
  print("delta0");
  print(delta0);
  print("delta1");
  print(delta1);
  print("v0");
  print(v0);
  print("v1");
  print(v1);
  print("log(v0 - v1)");
  print(log(v0 - v1));

  return log(v0 - v1);

}

vector log_inv_logit_diff_patched2(vector l, vector s, vector d0, vector d1) {
  // correctly handles Inf and -Inf values in d0 and d1, but the gradient can't handle that anyway
return log1m_exp((d0 - d1) ./ s) - log1p_exp((d0 - l) ./ s) - log1p_exp((l - d1) ./ s);
  }

vector fast_log_lik_common2(
  matrix lt, vector log_lambda, //vector log_E,
  matrix offset_lt, vector sign_offset_lt,
  int nobs, int no_time_points,
  array[] int idx_item_patient_vectorized,
  array[] int idx_item,
  array[] int idx_rater,
  vector delta_1, vector delta_r1, vector delta_r2, vector delta_nc,
  array[] int idx_is_one, array[] int idx_is_nc, array[] int idx_other) {

  vector[nobs] locations = to_vector(lt)[idx_item_patient_vectorized];
  vector[nobs] scales    = exp(log_lambda[idx_item]);// - log_E[idx_rater]);

  if (no_time_points == 2) {
    locations += sign_offset_lt .* to_vector(offset_lt)[idx_item_patient_vectorized];
  }

  vector[nobs] result;
  result[idx_is_one] = log1m_inv_logit(           (locations[idx_is_one] - delta_1) ./ scales[idx_is_one]);
  result[idx_other]  = log_inv_logit_diff_patched2(locations[idx_other], scales[idx_other], delta_r1, delta_r2);
  result[idx_is_nc]  = log_inv_logit(             (locations[idx_is_nc] - delta_nc) ./ scales[idx_is_nc]);

  return result;

  // return
  //   sum(log1m_inv_logit(           (locations[idx_is_one] - delta_1) ./ scales[idx_is_one])) +
    //   sum(log_inv_logit_diff_patched2(locations[idx_other], scales[idx_other], delta_r1, delta_r2)) +
    //   sum(log_inv_logit(             (locations[idx_is_nc] - delta_nc) ./ scales[idx_is_nc]));

}

vector fast_loglik_logistic(
  vector log_a, vector b, matrix lt, /*vector log_E,*/ vector log_lambda, matrix offset_lt, vector sign_offset_lt,
  int nobs, int nc, int no_time_points,
  array[] int scores0, array[] int scores1,
  array[] int idx_item_patient_vectorized,
  array[] int idx_item, array[] int idx_rater,
  array[] int idx_is_one, array[] int idx_is_nc, array[] int idx_other,
  vector gammas
) {

  // vector [nc - 1] gammas;
  // for (t in 1:(nc - 1)) {
    //   real tt = t; // Stan can only cast in this way and we need to cast to real to avoid integer division
    //   gammas[t] = log(tt / (nc - tt));
    // }

  vector[num_elements(idx_is_one)] delta_1  = exp(log_a[idx_rater[idx_is_one]]) * gammas[1]      + b[idx_rater[idx_is_one]];
  vector[num_elements(idx_is_nc )] delta_nc = exp(log_a[idx_rater[idx_is_nc ]]) * gammas[nc - 1] + b[idx_rater[idx_is_nc ]];

  vector[num_elements(idx_other)]  delta_r1 = exp(log_a[idx_rater[idx_other]]) .* gammas[scores0] + b[idx_rater[idx_other]];
  vector[num_elements(idx_other)]  delta_r2 = exp(log_a[idx_rater[idx_other]]) .* gammas[scores1] + b[idx_rater[idx_other]];

  return fast_log_lik_common2(
    lt, log_lambda, /*log_E,*/ offset_lt, sign_offset_lt, nobs, no_time_points, idx_item_patient_vectorized, idx_item, idx_rater,
    delta_1, delta_r1, delta_r2, delta_nc,
    idx_is_one, idx_is_nc, idx_other
  );

}

vector fast_loglik_free(
  matrix lt, /*/*vector log_E,*/ vector log_lambda, matrix free_thresholds, matrix offset_lt, vector sign_offset_lt,
  int nobs, int nc, int no_time_points,
  array[] int scores_by_rater_r1, array[] int scores_by_rater_r2,
  array[] int scores_by_rater_1, array[] int scores_by_rater_nc,
  array[] int idx_item_patient_vectorized,
  array[] int idx_item, array[] int idx_rater,
  array[] int idx_is_one, array[] int idx_is_nc, array[] int idx_other
) {

  vector[num_elements(idx_is_one)] delta_1  = to_vector(free_thresholds)[scores_by_rater_1];
  vector[num_elements(idx_is_nc )] delta_nc = to_vector(free_thresholds)[scores_by_rater_nc];

  vector[num_elements(idx_other)]  delta_r1 = to_vector(free_thresholds)[scores_by_rater_r1];
  vector[num_elements(idx_other)]  delta_r2 = to_vector(free_thresholds)[scores_by_rater_r2];

  return fast_log_lik_common2(
    lt, log_lambda, /*log_E,*/ offset_lt, sign_offset_lt, nobs, no_time_points, idx_item_patient_vectorized, idx_item, idx_rater,
    delta_1, delta_r1, delta_r2, delta_nc,
    idx_is_one, idx_is_nc, idx_other
  );

}

vector sta_qELGW_vec(vector p, vector shifts, vector directions, vector scales, vector shapes) {
  return shifts + directions .* scales .* log((1 - shapes .* log(0.5 - directions .* (p - 0.5))) ^ (1.0 ./ shapes) - 1);
}

vector fast_loglik_skew(
  vector log_a, vector b, vector shape, matrix lt, /*vector log_E,*/ vector log_lambda, matrix offset_lt, vector sign_offset_lt,
  int nobs, int nc, int no_time_points,
  array[] int scores0, array[] int scores1,
  array[] int idx_item_patient_vectorized,
  array[] int idx_item, array[] int idx_rater,
  array[] int idx_is_one, array[] int idx_is_nc, array[] int idx_other,
  vector threshold_probs
) {

  // vector [nc - 1] threshold_probs;
  // for (i in 1:(nc-1)) {
    //   real j = i; // cast to real to avoid integer division
    //   threshold_probs[i] = j / nc;
    // }

  vector [nobs] directions = sign_vec(shape)[idx_rater];
  vector [nobs] shapes     =      abs(shape)[idx_rater];

  // vector sta_qELGW_vec(vector p, vector shifts, vector directions, vector scales, vector shapes)
  vector[num_elements(idx_is_one)] delta_1  = sta_qELGW_vec(rep_vector(threshold_probs[1],      num_elements(idx_is_one)), b[idx_rater[idx_is_one]], directions[idx_is_one], exp(log_a[idx_rater[idx_is_one]]), shapes[idx_is_one]);
  vector[num_elements(idx_is_nc )] delta_nc = sta_qELGW_vec(rep_vector(threshold_probs[nc - 1], num_elements(idx_is_nc)),  b[idx_rater[idx_is_nc ]], directions[idx_is_nc ], exp(log_a[idx_rater[idx_is_nc ]]), shapes[idx_is_nc ]);

  vector[num_elements(idx_other)]  delta_r1 = sta_qELGW_vec(threshold_probs[scores0],                                      b[idx_rater[idx_other]],  directions[idx_other],  exp(log_a[idx_rater[idx_other]]),  shapes[idx_other]);
  vector[num_elements(idx_other)]  delta_r2 = sta_qELGW_vec(threshold_probs[scores1],                                      b[idx_rater[idx_other]],  directions[idx_other],  exp(log_a[idx_rater[idx_other]]),  shapes[idx_other]);

  return fast_log_lik_common2(
    lt, log_lambda, /*log_E,*/ offset_lt, sign_offset_lt, nobs, no_time_points, idx_item_patient_vectorized, idx_item, idx_rater,
    delta_1, delta_r1, delta_r2, delta_nc,
    idx_is_one, idx_is_nc, idx_other
  );

}

}
data{
  // booleans
  int<lower = 0, upper = 1> prior_only;
  int<lower = 0, upper = 1> debug;
  int<lower = 0, upper = 1> vary_lambda_across_patients;
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
  array[n_observed] int<lower = 1, upper = nr>             idx_rater  ;
  array[n_observed] int<lower = 1, upper = ni>             idx_item   ;
  array[n_observed] int<lower = 1, upper = np>             idx_patient;
  array[n_observed] int<lower = 1, upper = no_time_points> idx_time_point;

  array[nr] int<lower = 1, upper = no_rater_groups> idx_rater_group;

  // new stuff
  int<lower = 0, upper = n_observed> n_is_one;
  int<lower = 0, upper = n_observed> n_is_nc;
  int<lower = 0, upper = n_observed> n_is_other;
  array[n_is_one]   int<lower = 1, upper = n_observed> idx_is_one;
  array[n_is_nc]    int<lower = 1, upper = n_observed> idx_is_nc;
  array[n_is_other] int<lower = 1, upper = n_observed> idx_other;
  array[!use_free_logistic_thresholds ? n_is_other : 0] int scores0;
  array[!use_free_logistic_thresholds ? n_is_other : 0] int scores1;
  array[ use_free_logistic_thresholds ? n_is_other : 0] int scores_by_rater_r1;
  array[ use_free_logistic_thresholds ? n_is_other : 0] int scores_by_rater_r2;
  array[ use_free_logistic_thresholds ? n_is_one   : 0] int scores_by_rater_1;
  array[ use_free_logistic_thresholds ? n_is_nc    : 0] int scores_by_rater_nc;

  // logistic regression data
  int<lower = 0>                                                  np_log_reg;
  array [fit_logistic ? np_log_reg : 0] int<lower = 0, upper = 1> log_reg_outcomes;
  int<lower = 0>                                                  no_covariates;
  // this is not space efficient, but required for bernoulli_logit_glm
  matrix[fit_logistic ? np : 0, fit_logistic ? no_covariates : 0] design_matrix_covariates;

  // which observation in log_reg_outcomes corresponds to which patient?
    // basically lt[log_reg_np_idx, ] yields the lt that correspond to log_reg_outcomes
  array [np_log_reg] int<lower=1, upper=np> log_reg_np_idx;

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
  // real sd_mu_log_E_raw = no_rater_groups > 1 ? 1.0 : inv_sqrt(1 - inv(no_rater_groups));

  vector[nc-1] default_threshold_probs;
  vector[nc-1] default_thresholds; // <- 1:nc;
  int nt = nc - 1; // no thresholds

  for (t in 1:nt) {
    real tt = t; // Stan can only cast in this way and we need to cast to real to avoid integer division
    default_thresholds[t] = log(tt / (nc - tt));
    default_threshold_probs[t] = tt / nc;
  }

  array[store_predictions ? np : 0] int ones_np;
  if (store_predictions) {
    for (p in 1:np) {
      ones_np[p] = 1;
    }
  }

  array[n_observed] int <lower = 1, upper = np*ni> idx_item_patient_vectorized;
  array[n_missing]  int <lower = 1, upper = np*ni> idx_item_patient_vectorized_missing;
  for (o in 1:n_observed) idx_item_patient_vectorized[o]         = idx_patient[o]         + np * (idx_item[o]         - 1);
  for (o in 1:n_missing)  idx_item_patient_vectorized_missing[o] = idx_patient_missing[o] + np * (idx_item_missing[o] - 1);


  vector [no_time_points == 1 ? 0 : n_observed] sign_offset_lt;
  if (no_time_points == 2) {
    for (o in 1:n_observed) {
      // (2 * t - 3) maps {1, 2} to {-1, 1}
      sign_offset_lt[o] = 2 * idx_time_point[o] - 3;
    }
  }

  vector [n_missing > 0 && no_time_points == 1 ? 0 : n_missing] sign_offset_lt_missing;
  if (n_missing > 0 && no_time_points == 2) {
    for (o in 1:n_missing) {
      // (2 * t - 3) maps {1, 2} to {-1, 1}
      sign_offset_lt_missing[o] = 2 * idx_time_point_missing[o] - 3;
    }
  }

  // array[use_free_logistic_thresholds ? 0 : n_observed] int scores0;
  // array[use_free_logistic_thresholds ? 0 : n_observed] int scores1;
  // if (!use_free_logistic_thresholds) {
    //   for (o in 1:n_observed) {
      //     scores0[o] = x[o] - 1;
      //     scores1[o] = x[o];
      //   }
    // }
  //
    // array [use_free_logistic_thresholds ? n_observed : 0] int scores0_by_rater;
  // array [use_free_logistic_thresholds ? n_observed : 0] int scores1_by_rater;
  // if (use_free_logistic_thresholds) {
    //   for (o in 1:n_observed) {
      //     scores0_by_rater[o] = (x[o] - 1) * nr + idx_rater[o];
      //     scores1_by_rater[o] =  x[o]      * nr + idx_rater[o];
      //   }
    // }

  // only do this once.
  matrix[fit_logistic ? np_log_reg : 0, fit_logistic ? no_covariates : 0] design_matrix_covariates_observed = design_matrix_covariates[log_reg_np_idx, ];

}
parameters{

  // parameters
  matrix [np - 1, ni] lt_raw;         // item location values
  vector [ni - 1] log_lambda_raw; // item difficulty
  // vector [nr - 1] log_E_raw;          // rater competency

  // hyperparameters
  real mu_lt;
  real<lower = 0> sd_lt;
  real<lower = 0> sd_log_lambda;

  // vector           [no_rater_groups - 1] mu_log_E_raw;
  // vector           [no_rater_groups - 1] log_sd_log_E_raw;

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

  vector [ni] log_lambda;
  matrix [np, ni] lt;
  for (i in 1:ni) {
    lt[, i] = sd_lt * sum_to_zero_QR(lt_raw[, i], Q_np, np);
  }
  log_lambda = mu_log_lambda + sd_log_lambda * sum_to_zero_QR(log_lambda_raw, Q_ni, ni);

  // vector [no_rater_groups] mu_log_E;
  // vector [no_rater_groups] log_sd_log_E;
  // vector [nr] log_E;
  // if (no_rater_groups > 1) {
  //
  //   mu_log_E     = sum_to_zero_QR(mu_log_E_raw,     Q_no_rater_groups, no_rater_groups);
  //   log_sd_log_E = sum_to_zero_QR(log_sd_log_E_raw, Q_no_rater_groups, no_rater_groups);
  //   log_E = raw_to_constrained(log_E_raw, mu_log_E, log_sd_log_E, Q_nr, nr, no_rater_groups, idx_rater_group, rater_group_counts);
  //
  // } else {
  //
  //   mu_log_E[1]     = 0.0;
  //   log_sd_log_E[1] = 0.0;
  //   log_E = sum_to_zero_QR(log_E_raw, Q_nr, nr);
  //
  // }

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

  matrix[use_free_logistic_thresholds ? nr : 0, use_free_logistic_thresholds ? nc - 1 : 0] free_thresholds_mat;
  if (use_free_logistic_thresholds) {
    for (r in 1:nr)
      free_thresholds_mat[r, ] = free_thresholds[r]';
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
    // mu_log_E_raw     ~ normal(0, sd_mu_log_E_raw);
    // log_sd_log_E_raw ~ normal(0, sd_mu_log_E_raw);

    if (!use_free_logistic_thresholds)
      log_sd_log_a_raw ~ normal(0, sd_log_a_raw);//sd_mu_log_E_raw);

  }

  // priors
  to_vector(log_lambda_raw) ~ normal(0, sd_log_lambda_raw);
  to_vector(lt_raw)         ~ normal(0, sd_lt_raw);
  // log_E_raw ~ normal(0, sd_log_a_raw); // not a typo

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

  // CCT
  if (use_free_logistic_thresholds) {
    target += sum(fast_loglik_free(
      lt, /*log_E,*/ log_lambda, free_thresholds_mat, offset_lt, sign_offset_lt,
      n_observed, nc, no_time_points,
      scores_by_rater_r1, scores_by_rater_r2, scores_by_rater_1, scores_by_rater_nc,
      idx_item_patient_vectorized,
      idx_item, idx_rater,
      idx_is_one, idx_is_nc, idx_other
    ));
  } else if (use_skew_logistic_thresholds) {
    target += sum(fast_loglik_skew(
      log_a, b, threshold_shape, lt, /*log_E,*/ log_lambda, offset_lt, sign_offset_lt,
      n_observed, nc, no_time_points,
      scores0, scores1,
      idx_item_patient_vectorized,
      idx_item, idx_rater,
      idx_is_one, idx_is_nc, idx_other,
      default_threshold_probs
    ));
  } else {
    target += sum(fast_loglik_logistic(
      log_a, b, lt, /*log_E,*/ log_lambda, offset_lt, sign_offset_lt,
      n_observed, nc, no_time_points,
      scores0, scores1,
      idx_item_patient_vectorized,
      idx_item, idx_rater,
      idx_is_one, idx_is_nc, idx_other,
      default_thresholds
    ));
  }

  // CCT missing data
  if (n_missing > 0) {

    for (c in 1:nc) {

      int n_is_one_missing   = c == 1            ? n_missing : 0;
      int n_is_nc_missing    =           c == nc ? n_missing : 0;
      int n_is_other_missing = c != 1 && c != nc ? n_missing : 0;
      array[n_is_one_missing]   int idx_is_one_missing;
      array[n_is_nc_missing]    int idx_is_nc_missing;
      array[n_is_other_missing] int idx_other_missing;
      array[!use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_missing : 0] int scores0_missing;
      array[!use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_missing : 0] int scores1_missing;
      array[ use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_missing : 0] int scores_by_rater_r1_missing;
      array[ use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_missing : 0] int scores_by_rater_r2_missing;
      array[ use_free_logistic_thresholds && c == 1            ? n_is_one_missing   : 0] int scores_by_rater_1_missing;
      array[ use_free_logistic_thresholds &&           c == nc ? n_is_nc_missing    : 0] int scores_by_rater_nc_missing;

      if (c == 1)       for (o in 1:n_missing) idx_is_one_missing[o] = o;
      else if (c == nc) for (o in 1:n_missing) idx_is_nc_missing[o]  = o;
      else              for (o in 1:n_missing) idx_other_missing[o]  = o;


      if (use_free_logistic_thresholds) {

        if (c == 1)       for (o in 1:n_missing) scores_by_rater_1_missing[o]  = (c - 1) * nr + idx_rater_missing[o];
        else if (c == nc) for (o in 1:n_missing) scores_by_rater_nc_missing[o] = (c - 2) * nr + idx_rater_missing[o];
        else {
          for (o in 1:n_missing) {
            scores_by_rater_r1_missing[o] = (c - 2) * nr + idx_rater_missing[o];
            scores_by_rater_r2_missing[o] = (c - 1) * nr + idx_rater_missing[o];
          }
        }

        target += sum(fast_loglik_free(
          lt, /*log_E,*/ log_lambda, free_thresholds_mat, offset_lt, sign_offset_lt_missing,
          n_missing, nc, no_time_points,
          scores_by_rater_r1_missing, scores_by_rater_r2_missing, scores_by_rater_1_missing, scores_by_rater_nc_missing,
          idx_item_patient_vectorized_missing,
          idx_item_missing, idx_rater_missing,
          idx_is_one_missing, idx_is_nc_missing, idx_other_missing
        ));
      } else {

        if (c != 1 && c != nc) {
          for (o in 1:n_missing) {
            scores0_missing[o] = c - 1;
            scores1_missing[o] = c;
          }
        }

        if (use_skew_logistic_thresholds) {
          target += sum(fast_loglik_skew(
            log_a, b, threshold_shape, lt, /*log_E,*/ log_lambda, offset_lt, sign_offset_lt_missing,
            n_missing, nc, no_time_points,
            scores0_missing, scores1_missing,
            idx_item_patient_vectorized_missing,
            idx_item_missing, idx_rater_missing,
            idx_is_one_missing, idx_is_nc_missing, idx_other_missing,
            default_threshold_probs
          ));
        } else {
          target += sum(fast_loglik_logistic(
            log_a, b, lt, /*log_E,*/ log_lambda, offset_lt, sign_offset_lt_missing,
            n_missing, nc, no_time_points,
            scores0_missing, scores1_missing,
            idx_item_patient_vectorized_missing,
            idx_item_missing, idx_rater_missing,
            idx_is_one_missing, idx_is_nc_missing, idx_other_missing,
            default_thresholds
          ));
        }
      }
    }
  }

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

}
generated quantities {

  matrix<upper = 0>[debug ? nc : 0, debug ? n_observed : 0] log_probs;

  // this is incredibly memory intense!
  if (debug) {

    for (c in 1:nc) {

      int n_is_one_debug   = c == 1            ? n_observed : 0;
      int n_is_nc_debug    =           c == nc ? n_observed : 0;
      int n_is_other_debug = c != 1 && c != nc ? n_observed : 0;
      array[n_is_one_debug]   int idx_is_one_debug;
      array[n_is_nc_debug]    int idx_is_nc_debug;
      array[n_is_other_debug] int idx_other_debug;
      array[!use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_debug : 0] int scores0_debug;
      array[!use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_debug : 0] int scores1_debug;
      array[ use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_debug : 0] int scores_by_rater_r1_debug;
      array[ use_free_logistic_thresholds && c != 1 && c != nc ? n_is_other_debug : 0] int scores_by_rater_r2_debug;
      array[ use_free_logistic_thresholds && c == 1            ? n_is_one_debug   : 0] int scores_by_rater_1_debug;
      array[ use_free_logistic_thresholds &&           c == nc ? n_is_nc_debug    : 0] int scores_by_rater_nc_debug;

      if (c == 1)       for (o in 1:n_observed) idx_is_one_debug[o] = o;
      else if (c == nc) for (o in 1:n_observed) idx_is_nc_debug[o]  = o;
      else              for (o in 1:n_observed) idx_other_debug[o]  = o;

      if (use_free_logistic_thresholds) {

        if (c == 1)       for (o in 1:n_observed) scores_by_rater_1_debug[o]  = (c - 1) * nr + idx_rater[o];
        else if (c == nc) for (o in 1:n_observed) scores_by_rater_nc_debug[o] = (c - 2) * nr + idx_rater[o];
        else {
          for (o in 1:n_observed) {
            scores_by_rater_r1_debug[o] = (c - 2) * nr + idx_rater[o];
            scores_by_rater_r2_debug[o] = (c - 1) * nr + idx_rater[o];
          }
        }

        log_probs[c, ] = fast_loglik_free(
          lt, /*log_E,*/ log_lambda, free_thresholds_mat, offset_lt, sign_offset_lt,
          n_observed, nc, no_time_points,
          scores_by_rater_r1_debug, scores_by_rater_r2_debug, scores_by_rater_1_debug, scores_by_rater_nc_debug,
          idx_item_patient_vectorized,
          idx_item, idx_rater,
          idx_is_one_debug, idx_is_nc_debug, idx_other_debug
        )';

  } else {

    if (c != 1 && c != nc) {
      for (o in 1:n_observed) {
        scores0_debug[o] = c - 1;
        scores1_debug[o] = c;
      }
    }

    if (use_skew_logistic_thresholds) {
      log_probs[c, ] = fast_loglik_skew(
        log_a, b, threshold_shape, lt, /*log_E,*/ log_lambda, offset_lt, sign_offset_lt,
        n_observed, nc, no_time_points,
        scores0_debug, scores1_debug,
        idx_item_patient_vectorized,
        idx_item, idx_rater,
        idx_is_one_debug, idx_is_nc_debug, idx_other_debug,
        default_threshold_probs
      )';
        } else {
          log_probs[c, ] = fast_loglik_logistic(
            log_a, b, lt, /*log_E,*/ log_lambda, offset_lt, sign_offset_lt,
            n_observed, nc, no_time_points,
            scores0_debug, scores1_debug,
            idx_item_patient_vectorized,
            idx_item, idx_rater,
            idx_is_one_debug, idx_is_nc_debug, idx_other_debug,
            default_thresholds
          )';
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
  // This is not efficient, but probably also not what drives the total runtime
  for (o in 1:n_missing) {

    int p = idx_patient_missing[o];
    int i = idx_item_missing[o];
    int r = idx_rater_missing[o];

    real location = lt[p, i];
    real scale = exp(log_lambda[i]);// - log_E[r]);

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
