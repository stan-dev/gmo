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
g_alpha <- vb(stan_model(local_file),
              data=c(data, list(phi=phi)),
              seed=seed, init=alpha,
              output_samples=M*draws)

samples <- extract(g_alpha)
alpha <- list()
for (i in 1:(length(samples)-1)) {
  dims <- length(dim(samples[[i]]))
  if (dims == 1) {
    alpha[[i]] <- mean(samples[[i]])
  } else if (dims == 2) {
    alpha[[i]] <- apply(samples[[i]], 2, mean)
  } else {
    stop("")
  }
}
names(alpha) <- names(samples)[1:(length(samples)-1)]

# Check conditional approximation with warm start.
g_alpha <- vb(stan_model(local_file),
              data=c(data, list(phi=phi)),
              seed=seed, init=alpha,
              output_samples=M*draws)
