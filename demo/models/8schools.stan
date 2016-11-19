data {
  int<lower=0> J;          // number of schools
  real y[J];               // estimated treatment effect (school j)
  real<lower=0> sigma[J];  // std err of effect estimate (school j)

  int<lower=1> K;          // number of hyperparameters
}
parameters {
  vector[K] phi;           // mu, log_tau
  real eta[J];
}
transformed parameters {
  real mu;
  real<lower=0> tau;
  real theta[J];

  mu = phi[1];
  tau = exp(phi[2]);
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
  // adjust for change of vars:
  // log | d/d(log_tau) exp(log_tau) | = log_tau
  target += log(tau);
}
