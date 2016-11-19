data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  real y[N];
  int<lower=1, upper=J> group[N];
}
parameters {
  vector[K] phi;              // mu, log_sigma_y, log_sigma_alpha
  vector[J] alpha;
}
transformed parameters {
  real mu;
  real<lower=0> sigma_y;
  real<lower=0> sigma_alpha;
  mu = phi[1];
  sigma_y = exp(phi[2]);
  sigma_alpha = exp(phi[3]);
}
model {
  real y_pred[N];

  alpha ~ normal(0, sigma_alpha);
  for (n in 1:N) {
    y_pred[n] = mu + alpha[group[n]];
  }
  y ~ normal(y_pred, sigma_y);
  // adjust for change of vars:
  // log | d/d(log_sigma_y) exp(log_sigma_y) | = log_sigma_y
  // log | d/d(log_sigma_alpha) exp(log_sigma_alpha) | = log_sigma_alpha
  target += log(sigma_y);
  target += log(sigma_alpha);
}
