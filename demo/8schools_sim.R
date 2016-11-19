#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)

data_sim <- function(J, K, mu=1, tau=2) {
  set.seed(42)
  sigma <- rep(2, J)
  eta <- rnorm(J, 0, 1);
  theta <- mu + tau * eta
  y <- rep(NA, J)
  for (j in 1:J) {
    y[j] <- rnorm(1, theta[j], sigma[j])
  }
  return(list(J=J, K=K, mu=mu, tau=tau, theta=theta, sigma=sigma, y=y))
}

J <- 8
K <- 2
data_full <- data_sim(J, K)
data <- list(J=J, K=K, y=data_full$y, sigma=data_full$sigma)
#fit.true <- c(data_full$mu, log(data_full$tau))

fit.gmo <- gmo("models/8schools.stan", "models/8schools_local.stan", data=data)
fit.stan <- stan("models/8schools.stan", data=data)

print(fit.gmo)
print(fit.stan, pars="phi")
print(compute_z_score(fit.gmo, fit.stan))

library(microbenchmark)
benchmark <- microbenchmark(
  gmo=gmo("models/8schools.stan", "models/8schools_local.stan", data=data),
  nuts=stan("models/8schools.stan", data=data),
  vb=vb(stan_model("models/8schools.stan"), data=data),
  opt=optimizing(stan_model("models/8schools.stan"), data=data),
  times=5L
)
print(benchmark)
