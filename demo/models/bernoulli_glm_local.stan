data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;  // fixed as K = 2
  int y[N];
  int<lower=1, upper=J> group[N];
  vector[K] phi;              // mu, log_sigma_alpha
}
parameters {
  vector[J] alpha;
}
transformed parameters {
  real mu;
  real<lower=0> sigma_alpha;
  mu = phi[1];
  sigma_alpha = exp(phi[2]);
}
model {
  real y_pred[N];

  alpha ~ normal(0, sigma_alpha);
  for (n in 1:N) {
    y_pred[n] = mu + alpha[group[n]];
  }
  y ~ bernoulli_logit(y_pred);
  // adjust for change of vars:
  // log | d/d(log_sigma_alpha) exp(log_sigma_alpha) | = log_sigma_alpha
  target += log(sigma_alpha);
}
