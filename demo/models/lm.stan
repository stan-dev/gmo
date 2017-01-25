data {
  int<lower=1> N;
  int<lower=1> K;
  real y[N];
  matrix[N, K] X;
}
parameters {
  vector[K] phi;
  real<lower=0> sigma_y;
}
model {
  y ~ normal(X*phi, sigma_y);
}
