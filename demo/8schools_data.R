#!/usr/bin/env Rscript
# Shows how to use a single Stan program with branch logic
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)

data <- list(J = 8,
             K = 2,
             y = c(28,  8, -3,  7, -1,  1, 18, 12),
             sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit.gmo <- gmo("models/8schools_data.stan", data = data)
print(fit.gmo)

