functions {
#include /functions.stan
}
data{
  int<lower = 1> nr;
  int<lower = 1> ni;
  int<lower = 2> nc;
  int<lower = 1> np;
  int<lower = 1> nl;
  int<lower = 1> ncatR;
  int<lower = 1> ncatP;
  int x[nr, ni, np];

  int<lower = 1, upper = nl> latentMembership[ni];
  int<lower = 1, upper = ni> firstItemOfLv[nl];

  int<lower = 1, upper = ncatR> designR[nr];
  int<lower = 1, upper = ncatP> designP[np];

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
  matrix           [ni, np]    lt;        // item truths
  vector<lower = 0>[nr]        tmpAlpha;  // temporary rater scales
  vector           [nr]        beta;      // rater shifts
  matrix<lower = 0>[ni, np]    tmpKappa;  // temporary item difficulty
  vector<lower = 0>[nr]        tmpZeta;   // temporary rater competence
  matrix           [nl, np]    eta;       // factor scores
  vector<lower = 0>[ni]        reg;       // regression from latent variables on item difficulties

  // captures differences between groups of raters
  vector           [ncatR]     muBeta;
  // captures differences between patients
  matrix           [nl, ncatP] tmpEta;

  // hyperparameters
  real<lower = 0>              sdAlpha;
  real<lower = 0>              sdBeta;
  real<lower = 0>              sdZeta;
  vector<lower = 0>[np]        sdKappa;
  vector<lower = 0>[np]        sdLt;

}
transformed parameters {
  // identification constraints
  vector<lower = 0>[nr]     alpha = tmpAlpha ./ mean(tmpAlpha);
  vector<lower = 0>[nr]     zeta  = tmpZeta  ./ mean(tmpZeta);
  matrix<lower = 0>[ni, np] kappa;
  matrix           [nl, ncatP] muEta;

  for (p in 1:np) {
    kappa[, p] = tmpKappa[, p] ./ mean(tmpKappa[, p]);
  }

  for (l in 1:nl) {
    muEta[l ,]  = (tmpEta[l, ] - mean(tmpEta[l, ]));
  }
}
model{

  // hyperpriors
  to_vector(muEta) ~ normal(0, 10);  // mean of the latent variables per patient category
  muBeta  ~ normal(0, 10);           // mean of each rater category for threshold shift
  sdLt    ~ gamma(0.01, 1 / 0.01);   // standard deviation of item truths
  sdBeta  ~ gamma(0.01, 1.0 / 0.01); // standard deviation of rater shift
  sdAlpha ~ gamma(0.01, 1.0 / 0.01); // standard deviation of rater scale
  sdKappa ~ gamma(0.01, 1.0 / 0.01); // standard deviation of item difficulty
  sdZeta  ~ gamma(0.01, 1.0 / 0.01); // standard deviation of rater competence
  reg     ~ normal(0, 10);           // slopes from latent variable

  // priors on rater components
  alpha ~ gamma(square(muAlpha / sdAlpha), square(sdAlpha) / muAlpha); // rater scales
  beta  ~ normal(muBeta[designR], sdBeta);                             // rater shifts
  zeta  ~ gamma(square(muZeta / sdZeta), square(sdZeta) / muZeta);      // rater competence

  // factor structure over items
  to_vector(eta) ~ normal(to_vector(muEta[, designP]), 1); // latent variables
  for (p in 1:np) { // ideally this is vectorized over items!!!
    lt[, p]    ~ normal(reg .* (eta[, p])[latentMembership], sdLt[p]);
    kappa[, p] ~ gamma(square(muKappa / sdKappa[p]), square(sdKappa[p]) / muKappa);
  }

  // likelihood
  for (r in 1:nr) {
    vector[nt] delta; // create biased thresholds
    delta = alpha[r] * lambda + beta[r];
    for (i in 1:ni) {
      for (p in 1:np) {
        real scale = kappa[i, p] / zeta[r];
        x[r, i, p] ~ ordered_logistic(lt[i, p] / scale, delta ./ scale);
      }
    }
  }
}
