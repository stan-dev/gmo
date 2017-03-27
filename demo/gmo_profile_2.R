#!/usr/bin/env Rscript
# This is a benchmark profiling the bottlenecks of GMO, investigating
# why its implementation might be slower than its speed in principle.
# In particular, stan-optimize should be the main bottleneck, and each
# outer iteration should be roughly the time of one stan-optimize (the
# rest of the code should be relatively slow.) See conclusion notes in
# gmo.R.
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

## How many inner iterations (number of PLS calls) does glmer do?
print(fit.glmer@optinfo$feval)
## [1] 57

## What is the relationship between GMO inner iterations and
## stan-optimize?

set.seed(42)
N <- 10000
J <- 1000
data_full <- data_sim(N, J)
data <- list(N=N, J=J, K=data_full$K, y=data_full$y, group=data_full$group)

nreps <- 5
outer_iter <- rep(NA, nreps)
fit.gmo.time <- rep(NA, nreps)
fit.opt.time <- rep(NA, nreps)
for (t in 1:nreps) {
  # Benchmark GMO
  start <- Sys.time()
  fit.gmo <- gmo("models/bernoulli_glm.stan", "models/bernoulli_glm_local.stan",
                  data=data, inner_iter=100L)
  end <- Sys.time()
  fit.gmo.time[t] <- as.numeric(end - start)
  outer_iter[t] <- fit.gmo$eval_iter

  # Benchmark stan-optimize
  start <- Sys.time()
  fit.opt <- optimizing(stan_model("models/bernoulli_glm.stan"), data=data,
                        draws=100L * 5L)
  end <- Sys.time()
  fit.opt.time[t] <- as.numeric(end - start)
}
# How many outer iterations of GMO?
print(outer_iter)
# How many stan-optimizes for one GMO?
print(fit.gmo.time / fit.opt.time)
# How many stan-optimizes for each outer iteration?
print((fit.gmo.time / fit.opt.time) / outer_iter)

## setting: N=2000, J=50
 [1]  5  9  2  3  3  5  3  3  3  6  3  4  3  5  3 10  3  6  4  3  3  4  8  5  3

 [1]  7.429917 11.608673  4.252572  5.372748  6.497529  7.641569  5.306172
 [8]  4.985353  5.014481  8.418750  5.642028  3.030309  5.047641  7.622994
[15]  4.867438 13.438799  4.824375  9.821888  6.164957  6.057974  4.583134
[22]  5.979284 11.989245  8.211815  4.987972

 [1] 1.4859834 1.2898525 2.1262862 1.7909159 2.1658430 1.5283138 1.7687241
 [8] 1.6617844 1.6714937 1.4031250 1.8806760 0.7575773 1.6825471 1.5245989
[15] 1.6224792 1.3438799 1.6081250 1.6369813 1.5412393 2.0193246 1.5277113
[22] 1.4948210 1.4986557 1.6423630 1.6626572

print(mean((fit.gmo.time / fit.opt.time) / outer_iter))
[1] 1.613438

## setting: N=2000; J=1000
 [1] 23  9  4  5  7 12  4 17  5 18  6 12  7  3  4  6 10 10  9  4  6  6 20  5  5

 [1] 25.755481  8.957079  5.710978  5.721429  7.377425 14.061240  4.477241
 [8] 17.692113  6.193133 20.178277  6.645798 13.740446  7.601800  4.037732
[15]  5.010207  6.392143 10.205888  9.188929  8.099182  5.020506  6.755861
[22]  4.372215 21.815727  5.676842  6.135908

 [1] 1.1198035 0.9952310 1.4277445 1.1442857 1.0539179 1.1717700 1.1193101
 [8] 1.0407125 1.2386265 1.1210154 1.1076330 1.1450372 1.0859714 1.3459107
[15] 1.2525517 1.0653572 1.0205888 0.9188929 0.8999091 1.2551264 1.1259768
[22] 0.7287024 1.0907863 1.1353684 1.2271816

print(mean((fit.gmo.time / fit.opt.time) / outer_iter))
[1] 1.113496
# the more parameters, the closer it is to being 1-1?

 # N=2000; J=25
 [1] 3 3 4 3 3 3 8 3 4 3 3 3 3 4 3 3 3 3 3 3 4 3 3 4 5
> # How many stan-optimizes for one GMO?
> print(fit.gmo.time / fit.opt.time)
 [1]  5.457160  5.485563  5.196448  5.127766  6.487417  4.938206 14.475979
 [8]  5.354501  6.869711  5.492943  5.253335  6.409005  5.181836  6.949680
[15]  4.977357  5.180133  5.505145  5.255705  5.059860  2.189427  6.294726
[22]  5.563503  5.120851  6.566188  7.870227
> # How many stan-optimizes for each outer iteration?
> print((fit.gmo.time / fit.opt.time) / outer_iter)
 [1] 1.8190534 1.8285212 1.2991119 1.7092552 2.1624724 1.6460688 1.8094973
 [8] 1.7848337 1.7174278 1.8309811 1.7511115 2.1363350 1.7272786 1.7374199
[15] 1.6591191 1.7267109 1.8350482 1.7519018 1.6866200 0.7298091 1.5736816
[22] 1.8545009 1.7069503 1.6415470 1.5740454


## setting: N=10000; J=50

 [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3

 [1] 9.034366 8.249696 7.402007 9.092147 8.731370 9.099034 8.257271 8.712976
 [9] 9.580059 8.582806 8.432407 8.536784 8.897127 8.372369 8.813382 9.484358
[17] 8.594349 8.428439 9.188817 8.985175 9.097787 8.600278 9.380408 8.866178
[25] 8.563014

 [1] 3.011455 2.749899 2.467336 3.030716 2.910457 3.033011 2.752424 2.904325
 [9] 3.193353 2.860935 2.810802 2.845595 2.965709 2.790790 2.937794 3.161453
[17] 2.864783 2.809480 3.062939 2.995058 3.032596 2.866759 3.126803 2.955393
[25] 2.854338

print(mean((fit.gmo.time / fit.opt.time) / outer_iter))
[1] 2.919768
