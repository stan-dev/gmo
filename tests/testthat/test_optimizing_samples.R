library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)

data <- list(J = 8,
             K = 2,
             y = c(28,  8, -3,  7, -1,  1, 18, 12),
             sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

local_file <- "models/8schools_local.stan"
M <- 2
m <- 1
draws <- 10
phi <- c(2,5)
seed <- 42
alpha <- "random"

# Check conditional approximation.
g_alpha <- optimizing(stan_model(local_file),
                      data=c(data, list(phi=phi)),
                      seed=seed, init=alpha,
                      as_vector=FALSE,
                      draws=M*draws, constrained=FALSE)

# Check extraction of samples.
alpha_sims <- g_alpha$theta_tilde[(m-1)*draws + 1:draws, ]
