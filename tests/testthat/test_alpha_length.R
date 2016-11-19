library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(gmo)

get_stan_params <- function(object) {
  # Get vector of parameters which are not in the transformed block.
  stopifnot(is(object, "stanfit"))
  params <- grep("context__.vals_r", fixed = TRUE, value = TRUE,
                 x = strsplit(get_cppcode(get_stanmodel(object)), "\n")[[1]])
  params <- sapply(strsplit(params, "\""), FUN = function(x) x[[2]])
  params <- intersect(params, object@model_pars)
  return(params)
}

count_params <- function(object, params) {
  # Count number of parameters in object that belong to the params
  # vector.
  stopifnot(is(object, "stanfit"))
  return(as.integer(sum(sapply(object@par_dims[params], prod))))
}

data <- list(J = 8,
             K = 2,
             y = c(28,  8, -3,  7, -1,  1, 18, 12),
             sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

file <- "models/8schools.stan"
local_file <- "models/8schools_local.stan"
phi <- c(2,5)

full_model <- suppressMessages(stan(file, data = data, iter = 0, chains = 1))
local_model <- stan_model(local_file)
g_alpha <- vb(local_model, data=c(data,phi), iter=1)

pars <- get_stan_params(g_alpha)
count_params(g_alpha, pars)
