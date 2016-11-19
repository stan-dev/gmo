data {
  # TODO
  /*int<lower=0> J;          // number of schools*/
  /*int<lower=0> N_minibatch;                    */
  real<lower=0> J;          // number of schools
  real<lower=0> N_minibatch;
  real y[N_minibatch];               // estimated treatment effect (school j)
  real<lower=0> sigma[N_minibatch];  // std err of effect estimate (school j)

  int<lower=1> K;          // number of hyperparameters
  int<lower=0,upper=1> GMO_FLAG;
  vector[K * GMO_FLAG] fixed_phi;
}
transformed data {
  # TODO do i need this?
  real scale_factor;
  scale_factor = J / N_minibatch;
}
parameters {
  vector[K * (1 - GMO_FLAG)] phi;           // mu, log_tau
  real eta[J];
}
transformed parameters {
  real mu;
  real<lower=0> tau;
  real theta[J];

  if (GMO_FLAG) {
    mu = fixed_phi[1];
    tau = exp(fixed_phi[2]);
  }
  else {
    mu = phi[1];
    tau = exp(phi[2]);
  }
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
  // adjust for change of vars:
  // log | d/d(log_tau) exp(log_tau) | = log_tau
  target += log(tau);
  // scale joint density by N / N_minibatch
  // TODO i think this should really be N/N_minibatch * log p(x,z)
  target += log(scale_factor);
}
