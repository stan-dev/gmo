#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)
library(nlme)
library(microbenchmark)

data_sim <- function(N, D) {
  K <- D + 1
  X <- matrix(rnorm(N*D), ncol=D)
  X <- cbind(1, X)
  phi <- rep(5, K)
  y <- as.vector(X %*% phi + rnorm(N))
  return(list(N=N, K=K, phi=phi, X=X, y=y))
}

## exploring

set.seed(42)
N <- 2000
D <- 50
data_full <- data_sim(N, D)
data <- list(N=N, K=data_full$K, X=data_full$X, y=data_full$y)
data_lm <- data.frame(y=data_full$y, X=data_full$X[, 1:D+1]) # don't include intercept

fit.gmo <- gmo("models/lm.stan", "models/lm_local.stan",
                data=data, inner_iter=1L, iter=1L)
fit.opt <- optimizing(stan_model("models/lm.stan"), data=data,
                      as_vector=FALSE, draws=1L*5L, constrained=FALSE)
fit.lm <- lm(y ~ ., data=data_lm)
fit.lme <- gls(y ~ ., data=data_lm)

print(fit.gmo$par)
print(fit.opt)
print(fit.lm)
print(fit.lme)

## Benchmarks

set.seed(42)
N <- 10000
D <- 1000
data_full <- data_sim(N, D)
data <- list(N=N, K=data_full$K, X=data_full$X, y=data_full$y)
data_lm <- data.frame(y=data_full$y, X=data_full$X[, 1:D+1]) # don't include intercept

nreps <- 5
fit.gmo.time <- rep(NA, nreps)
fit.gmo2.time <- rep(NA, nreps)
fit.opt.time <- rep(NA, nreps)
fit.lm.time <- rep(NA, nreps)
fit.lme.time <- rep(NA, nreps)
for (t in 1:nreps) {
  # Benchmark GMO
  start <- Sys.time()
  fit.gmo <- gmo("models/lm.stan", "models/lm_local.stan",
                  data=data, inner_iter=1L, iter=1L)
  end <- Sys.time()
  fit.gmo.time[t] <- as.numeric(end - start)

  # Benchmark GMO with other inner_iter configs
  start <- Sys.time()
  fit.gmo <- gmo("models/lm.stan", "models/lm_local.stan",
                  data=data, inner_iter=10L, iter=1L)
  end <- Sys.time()
  fit.gmo2.time[t] <- as.numeric(end - start)

  # Benchmark stan-optimize
  start <- Sys.time()
  fit.opt <- optimizing(stan_model("models/lm.stan"), data=data,
                        as_vector=FALSE, draws=1L*5L, constrained=FALSE)
  end <- Sys.time()
  fit.opt.time[t] <- as.numeric(end - start)

  # Benchmark lm
  start <- Sys.time()
  fit.lm <- lm(y ~ ., data=data_lm)
  end <- Sys.time()
  fit.lm.time[t] <- as.numeric(end - start)

  # Benchmark lme
  start <- Sys.time()
  fit.lme <- gls(y ~ ., data=data_lm)
  end <- Sys.time()
  fit.lme.time[t] <- as.numeric(end - start)
}
print(fit.gmo.time)
print(fit.gmo2.time)
print(fit.opt.time)
print(fit.lm.time)
print(fit.lme.time)

# N = 10000; D=1000
#> print(fit.gmo.time) # inner_iter=1L, draws=5L
#[1] 2.138441 1.660907 1.738494 2.178171 1.803406
#> print(fit.gmo2.time) # inner_iter=10L, draws=5L
#[1] 4.579119 4.580811 5.023243 6.996433 4.821303
#> print(fit.opt.time)
#[1] 2.010840 2.163910 2.349114 2.529385 2.347591
#> print(fit.lm.time)
#[1] 5.852461 6.109148 9.949945 6.245936 6.042907
#> print(fit.lme.time)
#[1]  6.150059  6.679910 10.822603  6.808194  6.676801
# Note although GMO is actually faster here than stan-optimize, this is
# because it runs stan-optimize on just the local model, which is 1
# parameter rather than the total of 1002.
