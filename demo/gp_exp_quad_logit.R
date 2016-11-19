#!/usr/bin/env Rscript
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)

data <- read.table("data/wells.dat", sep=" ", header=TRUE)
data <- data[complete.cases(data), ]
D <- 1
N <- 100
K <- 2
x <- array(as.numeric(data$arsenic[1:N]), c(N, 1))
y <- as.numeric(data$switch[1:N])
dat <- list(D=D, N=N, K=K, x=x, y=y)

fit.gmo <- gmo("models/gp_exp_quad_logit.stan", "models/gp_exp_quad_logit_local.stan", dat=dat, eta=0.1, iter=100L, tol=1e-15)
fit.stan <- stan("models/gp_exp_quad_logit.stan")
gpmodel <- stan_model(file = "models/gp_exp_quad_logit.stan")
fit.advi <- vb(gpmodel, data=dat)

print(fit.gmo)
print(fit.stan, pars="phi")
print(fit.advi, pars="phi")
print(compute_z_score(fit.gmo, fit.stan))
print(compute_z_score(fit.gmo, fit.advi))

# Plot latent function p(alpha | y, phi).
nsamples <- nrow(fit.gmo$sims)
dat <- data.frame(
  x=as.vector(x),
  eta=as.vector(fit.gmo$sims),
  group=rep(1:N, nsamples))

#library(ggplot)
##ggplot(dat, aes(x=x, y=eta, group=group)) + geom_line()
#ggplot(dat, aes(x=x, y=eta, group=1)) + geom_line()

# Plot mean of latent function p(alpha | y, phi).
dat <- data.frame(x=as.vector(x), eta=colMeans(fit.gmo$sims))
ggplot(dat, aes(x=x, y=eta)) + geom_line()


fit.gmo <- gmo_approx("models/gp_exp_quad_logit.stan", "models/gp_exp_quad_logit_local.stan", dat=dat, eta=0.1, iter=100L, tol=1e-15)



# Stan conditional on GMO parameters.
data <- c(dat, list(phi=fit.gmo$par))
fit.stan <- stan("models/gp_exp_quad_logit_local.stan", data=data)

sims <- extract(fit.stan, permuted = FALSE)
#sims <- sims[, 1, 1:100] # eta's
sims <- sims[, 1, 103:202] # f's
ggplot_dat <- data.frame(x=as.vector(x), f=colMeans(sims))
ggplot(ggplot_dat, aes(x=x, y=f)) + geom_line()

# Full Bayes.
fit.stan <- stan("models/gp_exp_quad_logit.stan", data=dat, chains=1)
