#' Uncertainty estimate of marginal posterior parameters
#'
#' Calculate a covariance matrix of the marginal posterior
#' parameters, based on an unbiased estimate of the inverse of the
#' observed Fisher information.
#'
#' @param fit.gmo
#'   An object of reference class \code{"gmo"}.
#'
#' @export
est_covariance <- function(fit.gmo, ...) {
  #density_sims <- fit.gmo$.log_p(alpha_sims, m, g_flag=TRUE)
  #if (draws < 25) {
  #  log_r <- density_sims$log_p - density_sims$log_g
  #} else {
  #  log_r <- psislw(density_sims$log_p - density_sims$log_g)$lw_smooth
  #}
  #max_log_r <- max(log_r)
  #r <- exp(log_r - max_log_r)

  ## Note that weighted.mean normalizes the importance ratios
  #log_p_phi <- max_log_r + log(mean(r))
  #grad_log_p_phi <- apply(density_sims$grad_log_p, 2, weighted.mean, r)

  #grad_log_joint_sims
  #hess_log_joint_sims
  #log_p_phi <- mean()
  #grad_log_p_phi
  #hess_log_p_phi
  #return(diag(length(phi)))
  return(matrix())
}
