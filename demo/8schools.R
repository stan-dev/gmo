#!/usr/bin/env Rscript
# Benchmarks based on
# 2.6 GHz, Intel Core i5
# Commit: 47b12887210892e69227320bd1d8c05893fee5b4
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)

data <- list(J = 8,
             K = 2,
             y = c(28,  8, -3,  7, -1,  1, 18, 12),
             sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit.gmo <- gmo("models/8schools.stan", "models/8schools_local.stan", data=data)
fit.stan <- stan("models/8schools.stan", data=data)

print(fit.gmo)
# TODO
#fit.lme <- lmer(data$y ~ (1|data$group), REML=FALSE)
print(fit.stan, pars="phi")
print(compute_z_score(fit.gmo, fit.stan))

fit.gmo$.cond_infer(c(fit.gmo$data, list(fixed_phi=fit.gmo$par)))

library(microbenchmark)
benchmark <- microbenchmark(
  gmo=gmo("models/8schools.stan", "models/8schools_local.stan", data=data),
  nuts=stan("models/8schools.stan", data=data),
  vb=vb(stan_model("models/8schools.stan"), data=data),
  opt=optimizing(stan_model("models/8schools.stan"), data=data),
  times=5L
)
print(benchmark)
## Unit: milliseconds
##  expr       min        lq      mean    median        uq       max neval
##   gmo  484.5622  490.6744  630.4749  555.4390  687.3984  934.3006     5
##  nuts 4364.8530 4374.0095 4407.5686 4394.1122 4438.9821 4465.8863     5
##    vb  258.2433  283.7053  311.1285  285.5008  351.7344  376.4589     5
##   opt  180.9541  184.2089  187.5173  186.7766  191.4394  194.2075     5
