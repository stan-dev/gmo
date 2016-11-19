#' A Reference Class for optimization.
#'
#' @field eta
#'   Double, constant scale factor for learning rate.
#' @field grad_prev
#'   A numeric vector, previous gradient value.
#' @field I_est
#'   A numeric vector, an estimate of the inverse Fisher information.
Opt <- setRefClass("opt",
  fields=list(
    eta="numeric",
    grad_prev="numeric",
    I_est="numeric"
  ), methods=list(
    initialize = function(eta) {
      eta <<- eta
      grad_prev <<- 0
      I_est <<- 0
    },
    update_params = function(par, grad, iter) {
      "
      Update parameters using Nesterov's accelerated gradient with
      RMSprop.
      "
      I_est <<- 0.9*I_est + 0.1*grad^2
      lr <- eta/(sqrt(I_est) + 1)*1/sqrt(iter)
      grad_prev <<- 0.9*grad_prev + lr*grad

      par <- par + grad_prev
      par[par > 25] <- 25 # TODO this is an arbitrary truncation
      par[par < -25] <- -25
      return(par)
    }
  )
)

#' A Reference Class for diagnosing convergence.
#'
#' @field par_prev
#'   A numeric vector, previous checked parameter value.
#' @field fn_prev
#'   Double, previous checked function value.
#' @field grad_prev
#'   A numeric vector, previous checked gradient value.
#' @field tol_rel
#'   Double, relative tolerance for signaling convergence.
Diagnostic <- setRefClass("diagnostic",
  fields=list(
    par_prev="numeric",
    fn_prev="numeric",
    tol_param="numeric",
    tol_obj="numeric",
    tol_rel_obj="numeric",
    tol_grad="numeric",
    tol_rel_grad="numeric"
  ), methods=list(
    initialize = function(tol) {
      "
      Initializes previous checked values at 0.

      @param tol
        Tolerance threshold used for all tolerance checks.
      "
      par_prev <<- 0
      fn_prev <<- 0
      tol_param <<- tol
      tol_obj <<- tol
      tol_rel_obj <<- tol
      tol_grad <<- tol
      tol_rel_grad <<- tol
      # Defaults used in Stan's L-BFGS.
      #tol_param <<- 1e-8
      #tol_obj <<- 1e-12
      #tol_rel_obj <<- 1e+4 * .Machine$double.eps
      #tol_grad <<- 1e-8
      #tol_rel_grad <<- 1e+3 * .Machine$double.eps
      # A toned down version of Stan's defaults.
      #tol_param <<- 1e-3
      #tol_obj <<- 1e-7
      #tol_rel_obj <<- 1e+9 * .Machine$double.eps
      #tol_grad <<- 1e-3
      #tol_rel_grad <<- 1e+8 * .Machine$double.eps
    },
    check_converge = function(par, fn, grad) {
      "
      Check if the difference in parameters, function value, or
      gradient values reached below tolerances.
      "
      par_check      <- norm(par - par_prev, type="2")
      fn_check       <- abs(fn - fn_prev)
      rel_fn_check   <- fn_check/max(abs(fn), abs(fn_prev), 1)
      grad_check     <- norm(grad, type="2")
      # This is Stan's relative grad check for LBFGS but without the
      # Hessian.
      rel_grad_check <- (grad %*% grad)/max(max(abs(grad)), 1)

      par_prev  <<- par
      fn_prev   <<- fn

      par_flag      <- par_check      < tol_param    & is.finite(par_check)
      fn_flag       <- fn_check       < tol_obj      & is.finite(fn_check)
      rel_fn_flag   <- rel_fn_check   < tol_rel_obj  & is.finite(rel_fn_check)
      grad_flag     <- grad_check     < tol_grad     & is.finite(grad_check)
      rel_grad_flag <- rel_grad_check < tol_rel_grad & is.finite(rel_grad_check)
      return(c(par_flag, fn_flag, rel_fn_flag, grad_flag, rel_grad_flag))
    }
  )
)

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

#' alpha_mean
#'
#' Take alpha from mean of samples. (This is used for warm starts.)
#'
#' @export
alpha_mean <- function(g_alpha) {
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
  return(alpha)
}

#' @export
compute_z_score <- function(x, ...) UseMethod("compute_z_score")

#' compute_z_score
#'
#' Calculate z-score of parameter estimates from a vector using the
#' the true posterior mean and standard deviation (based on MCMC).
#'
#' @export
compute_z_score.default <- function(par, fit.stan) {
  fit.stan <- extract(fit.stan)
  par_stan <- apply(fit.stan$phi, 2, mean)
  par_stan_sd <- apply(fit.stan$phi, 2, sd)
  z_scores <- (par - par_stan) / par_stan_sd
  return(z_scores)
}

#' compute_z_score
#'
#' Calculate z-score of parameter estimates from GMO using the
#' the true posterior mean and standard deviation (based on MCMC).
#'
#' @param gmo_obj
#'   A gmo object with estimated hyperparameters.
#' @param stanfit
#'   A stanfit object with MCMC samples.
#'
#' @return
#' For i=1,...,length(phi),
#'   z_i = ( (phi^{GMO}_i - mean^{MCMC}(phi_i) ) / sd^{MCMC}(phi_i)
#'
#' @export
compute_z_score.gmo <- function(fit.gmo, fit.stan) {
  fit.stan <- extract(fit.stan)
  par_stan <- apply(fit.stan$phi, 2, mean)
  par_stan_sd <- apply(fit.stan$phi, 2, sd)
  z_scores <- (fit.gmo$par - par_stan) / par_stan_sd
  return(z_scores)
}

#' extract_stan_params
#'
#' Returns vector of parameters.
#' @export
extract_stan_params <- function(fit.stan) {
  fit.stan <- extract(fit.stan)
  par_stan <- apply(fit.stan$phi, 2, mean)
  return(par_stan)
}

#' extract_lme_params
#'
#' Returns vector of fixed effects and variance components.
#' @export
extract_lme_params <- function(fit.lme) {
  return(c(attributes(fit.lme)$beta,
           as.data.frame(VarCorr(fit.lme))$sdcor))
}
