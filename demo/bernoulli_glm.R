#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)
library(lme4)
library(microbenchmark)

data_sim <- function(N, J, mu=1, sigma_alpha=1) {
  K <- 2
  group <- sort(rep(1:J, length=N))
  alpha <- rnorm(J, 0, sigma_alpha)
  y_pred <- mu + alpha[group]
  y <- rbinom(N, 1, 1 / (1 + exp(-y_pred)))
  return(list(N=N, J=J, K=K, mu=mu, sigma_alpha=sigma_alpha, group=group, alpha=alpha, y=y))
}

##
# Benchmarking hierarchical logistic regression as # data points and #
# groups increase.
##

set.seed(42)

Ns <- c(2e2, 2e3, 1e4, 2e4)
Js <- c(1e1, 1e2, 5e2, 1e3)
benchmark <- list()
for (ntrial in 1:length(Ns)) {
  N <- Ns[ntrial]
  J <- Js[ntrial]
  data_full <- data_sim(N, J)
  data <- list(N=N, J=J, K=data_full$K, y=data_full$y, group=data_full$group)
  data_glmer <- data.frame(y=data$y, group=data$group)

  benchmark[[ntrial]] <- microbenchmark(
    nuts=stan("models/bernoulli_glm.stan", data=data),
    gmo=gmo("models/bernoulli_glm.stan", "models/bernoulli_glm_local.stan", data=data),
    opt=optimizing(stan_model("models/bernoulli_glm.stan"), data=data),
    glmer=glmer(y ~ (1|group), family="binomial", data=data_glmer),
    times=5L
  )
  print(benchmark[[ntrial]])
}
print(benchmark)
#[[1]]
#Unit: milliseconds
#  expr        min         lq       mean     median         uq        max neval
#  nuts 2622.01993 2830.50469 2865.80767 2942.29926 2960.66631 2973.54816     5
#   gmo  431.58659  476.98951  679.80280  570.76744  651.67474 1267.99570     5
#   opt  191.97145  194.98640  198.60582  199.07556  200.16966  206.82603     5
# glmer   46.26443   47.29814   58.44586   48.78835   53.38333   96.49504     5

#[[2]]
#Unit: milliseconds
#  expr        min         lq       mean     median         uq        max neval
#  nuts 17401.9923 17558.0543 18785.2996 18449.2936 19371.7965 21145.3614     5
#   gmo   652.7802   769.3413   893.2921   940.1495  1048.1690  1056.0204     5
#   opt   198.3568   472.8487   433.3641   492.5632   493.3170   509.7349     5
# glmer   123.7474   126.5184   130.4702   127.0052   133.1286   141.9511     5

#[[3]]
#Unit: milliseconds
#  expr        min         lq       mean     median         uq        max neval
#  nuts 75389.0363 77537.0069 81689.8771 78528.4142 84802.9434 92191.9848     5
#   gmo  1775.1070  5086.5572  5718.4593  6932.5627  7187.7950  7610.2744     5
#   opt   198.6876   199.9336   209.6445   209.3155   217.0524   223.2334     5
# glmer   741.7271   769.7428   786.1662   773.6309   780.5135   865.2167     5

#[[4]]
#Unit: milliseconds
#  expr         min         lq       mean     median         uq        max neval
#  nuts 193512.8733 204496.443 238477.019 205641.526 256332.828 332401.422     5
#   gmo  16577.8806  25892.483  30009.314  30411.099  37063.772  40101.338     5
#   opt   3242.3529   3410.763   3638.031   3450.882   3662.200   4423.957     5
# glmer    970.2502   1086.679   1199.461   1160.609   1285.677   1494.088     5

##
# Analyzing GMO when dim(alpha) is high.
##

set.seed(42)
N <- 2000
J <- 50
data_full <- data_sim(N, J)
data <- list(N=N, J=J, K=data_full$K, y=data_full$y, group=data_full$group)

fit.gmo <- gmo("models/bernoulli_glm.stan", "models/bernoulli_glm_local.stan",
                data=data, inner_iter=100L, tol=1e-5)
fit.stan <- stan("models/bernoulli_glm.stan", data=data)
#fit.opt <- optimizing(stan_model("models/bernoulli_glm.stan"), data=data)
data_glmer <- data.frame(y=data$y, group=data$group)
fit.glmer <- glmer(y ~ (1|group), family="binomial", data=data_glmer)

print(c(fit.gmo$par[1], exp(fit.gmo$par[2])))
print(fit.stan, pars=c("mu", "sigma_alpha"))
#print(fit.opt)
print(fit.glmer)

# Note: Both alpha and phi are approximately Gaussian posteriors
mu_draws <- extract(fit.stan)$phi[,1]
hist(mu_draws)
alpha_one_draws <- extract(fit.stan)$alpha[,1]
hist(alpha_one_draws)
