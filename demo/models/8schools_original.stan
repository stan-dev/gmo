data {
  int<lower=0> J;          // number of schools
  real y[J];               // estimated treatment effect (school j)
  real<lower=0> sigma[J];  // std err of effect estimate (school j)
}
parameters {
  real mu;
  real eta[J];
  real<lower=0> tau;
}
transformed parameters {
  real theta[J];

  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
