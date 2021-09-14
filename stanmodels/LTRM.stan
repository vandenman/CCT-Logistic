functions {
#include /functions.stan
}
data{
  int<lower = 1> nr;
  int<lower = 1> ni;
  int<lower = 2> nc;
  int x[nr, ni];
}
transformed data{

  real muAlpha = 1.0; // this is restriced to 1 anyway
  real muKappa = 1.0; // this is restriced to 1 anyway
  real muZeta  = 1.0; // this is restriced to 1 anyway
  vector[nc-1] lambda; // <- 1:nc;
  int nt; // no thresholds
  nt = nc - 1;

  for (t in 1:nt) {
    real tt = t; // Stan can only cast in this way and we need to cast to real to avoid integer division
    lambda[t] = log(tt / (nc - tt));
  }
}
parameters{
  vector[ni]            lt;       // item truths
  vector<lower = 0>[nr] tmpAlpha; // temporary rater scales
  vector[nr]            beta;     // rater shifts
  vector<lower = 0>[ni] tmpKappa; // temporary item difficulty
  vector<lower = 0>[nr] tmpZeta;  // temporary rater competence

  // hyperparameters
  real            muLt;
  real            muBeta;
  real<lower = 0> sdLt;
  real<lower = 0> sdBeta;
  real<lower = 0> sdAlpha;
  real<lower = 0> sdKappa;
  real<lower = 0> sdZeta;

}
transformed parameters {

  // identification constraints
  vector<lower = 0>[nr] alpha = tmpAlpha ./ mean(tmpAlpha);
  vector<lower = 0>[ni] kappa = tmpKappa ./ mean(tmpKappa);
  vector<lower = 0>[nr] zeta  = tmpZeta  ./ mean(tmpZeta);

}
model{
  // hyperpriors
  muLt    ~ normal(0, 10);            // mean of item truths
  muBeta  ~ normal(0, 10);            // mean of rater shifts
  sdLt    ~ gamma(0.01, 1 / 0.01);    // standard deviation of item truths
  sdBeta  ~ gamma(0.01, 1.0 / 0.01);  // standard deviation of rater shift
  sdAlpha ~ gamma(0.01, 1.0 / 0.01);  // standard deviation of rater scale
  sdKappa ~ gamma(0.01, 1.0 / 0.01);  // standard deviation of item difficulty
  sdZeta  ~ gamma(0.01, 1.0 / 0.01);

  // priors
  lt    ~ normal(muLt, sdLt);     // item truths
  beta  ~ normal(muBeta, sdBeta); // rater shifts
  alpha ~ gamma(square(muAlpha / sdAlpha), square(sdAlpha) / muAlpha); // rater scale
  kappa ~ gamma(square(muKappa / sdKappa), square(sdKappa) / muKappa); // item difficulty
  zeta  ~ gamma(square(muZeta  / sdZeta),  square(sdZeta)  / muZeta);  // rater ability

  // likelihood
  for (r in 1:nr) {
    vector[nt] delta; // create biased thresholds
    delta = alpha[r] * lambda + beta[r];
    for (i in 1:ni) {
      real scale = kappa[i] / zeta[r];
      x[r, i] ~ ordered_logistic(lt[i] / scale, delta ./ scale);
    }
  }
}
