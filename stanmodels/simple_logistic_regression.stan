data{
  int<lower = 0>                 np;
  int<lower = 0>                 no_covariates;
  int<lower = 0, upper = 1>      log_reg_outcomes[np];
  matrix[np, no_covariates]      design_matrix;
}
parameters{
  real log_reg_intercept;
  vector[no_covariates] log_reg_slopes;
}
model{
  // priors
  log_reg_intercept ~ std_normal();
  log_reg_slopes    ~ std_normal();

  log_reg_outcomes ~ bernoulli_logit_glm(design_matrix, log_reg_intercept, log_reg_slopes);
}
