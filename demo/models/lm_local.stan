data {
  int<lower=1> N;
  int<lower=1> K;
  real y[N];
  matrix[N, K] X;
  vector[K] phi;
}
parameters {
  real<lower=0> sigma_y;
}
model {
  y ~ normal(X*phi, sigma_y);
}
