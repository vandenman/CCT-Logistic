data{
  int<lower = 1> n_arrays;
  int<lower = 1> total_length;
  array[n_arrays] int<lower = 1, upper = total_length> starts;
  array[n_arrays] int<lower = 1, upper = total_length> stops;
}
parameters{
  vector[total_length] y;
}
transformed parameters {
  // based on https://mc-stan.org/docs/reference-manual/ordered-vector.html
  real log_abs_det = 0.0;
  vector[total_length] x;
  for (i in 1:n_arrays) {
    x[starts[i]] = y[starts[i]];
    for (k in starts[i]+1:stops[i]) {
      x[k] = x[k-1] + exp(y[k]);
      log_abs_det += y[k];
    }
  }
}
model {

  x ~ std_normal();
  target += log_abs_det;

}
