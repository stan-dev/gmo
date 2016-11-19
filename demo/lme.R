#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)
library(lme4)

data_sim <- function(n, p, q, mu=1, sigma_y=2, sigma_alpha=2) {
  # TODO p is not used for design matrix
  group <- sort(rep(1:q, length=n))
  alpha <- rnorm(q, 0, sigma_alpha)
  y_pred <- mu + alpha[group]
  y <- rnorm(n, y_pred, sigma_y)
  return(list(n=n, q=q, p=p, mu=mu, sigma_y=sigma_y, sigma_alpha=sigma_alpha, group=group, alpha=alpha, y=y))
}

set.seed(42)
n <- 200
p <- 3
q <- 10
data_full <- data_sim(n, p, q)
data <- list(n=n, p=p, q=q, y=data_full$y, group=data_full$group)
#fit.true <- c(data_full$mu, log(data_full$sigma_y), log(data_full$sigma_alpha))

fit.gmo <- gmo("models/hlm.stan", "models/hlm_local.stan", data=data,
               inner_iter=100L, tol=1e-5)
fit.stan <- stan("models/hlm.stan", data=data)
fit.lme <- lmer(data$y ~ (1|data$group), REML=FALSE)
fit.lme <- extract_lme_params(fit.lme)
fit.lme[2:3] <- log(fit.lme[2:3])

print(fit.gmo)
print(extract_stan_params(fit.stan))
print(fit.lme)
print(compute_z_score(fit.gmo, fit.stan))
print(compute_z_score(fit.lme, fit.stan))
