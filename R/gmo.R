#' Gradient-based marginal optimization
#'
#' Maximize the marginal posterior with respect to specified
#' parameters, with nuisance parameters marginalized out.
#'
#' @param file
#'   A character string or a connection that R supports specifying the
#'   Stan model specification in Stan's modeling language.
#' @param local_file
#'   A character string or a connection that R supports specifying the
#'   Stan model specification in Stan's modeling language.
#' @param full_model
#'   If provided, an object of class 'stanfit' that makes it unnecessary
#'   to pass 'file' or 'local_file'
#' @param data
#'   A named ‘list’ or ‘environment’ providing the data for the model
#'   or a character vector for all the names of objects used as data.
#'   See \code{\link[rstan]{stan}} for more details.
#' @param method
#'   A character string naming the conditional inference:
#'   "laplace"
#' @param init
#'   A numeric vector of length the number of hyperparameters.
#' @param draws
#'   A positive integer, number of draws to calculate stochastic gradient.
#' @param iter
#'   A positive integer, the maximum number of iterations.
#' @param inner_iter
#'   A positive integer, the number of iterations after each
#'   conditional inference.
#' @param cond_iter
#'   A positive integer, the maximum number of iterations for the
#'   conditional inference. Default is to run until convergence.
#' @param eta
#'   Double, constant scale factor for learning rate.
#' @param tol
#'   Double, tolerance for signaling convergence.
#' @param seed
#'   The seed, a positive integer, for random number generation of
#'   Stan. The default is generated from 1 to the maximum integer
#'   supported by R so fixing the seed of R's random number generator
#'   can essentially fix the seed of Stan.  When multiple chains are
#'   used, only one seed is needed, with other chains' seeds being
#'   generated from the first chain's seed to prevent dependency among
#'   the random number streams for the chains. When a seed is
#'   specified by a number, ‘as.integer’ will be applied to it.  If
#'   ‘as.integer’ produces ‘NA’, the seed is generated randomly.  We
#'   can also specify a seed using a character string of digits, such
#'   as ‘"12345"’, which is converted to integer.
#'
#' @return
#' An object of reference class \code{"gmo"}. It is a list containing
#' the following components:
#' \item{par}{a vector of optimized parameters}
#' \item{cov}{estimated covariance matrix at \code{par}}
#' \item{sims}{\code{draws * inner_iter} many samples from the last
#'     approximation to the conditional posterior, p(alpha | y, phi)}
#'
#' @import methods
#' @importFrom rstan stan optimizing vb sampling constrain_pars log_prob grad_log_prob stan_model get_stanmodel
#' @importFrom loo psislw
#' @export
gmo <- function(file, local_file, full_model, data,
  method=c("laplace"), init="random",
  draws=5L, iter=100L, inner_iter=10L, cond_iter=NA, eta=1,
  tol=1e-3, seed=1234L) {
  if (missing(full_model)) {
    full_model <- suppressMessages(
                  stan(file, data = c(data, list(GMO_FLAG = FALSE, fixed_phi = double())),
                       chains = 0, iter = 1))
  }
  else stopifnot(is(full_model, "stanfit"))
  local_model <- if (!missing(local_file)) stan_model(local_file) else
                 get_stanmodel(full_model)
  obj <- GMO$new(
    calc_log_p="exact",
    full_model=full_model,
    local_model=local_model,
    data=data,
    method=match.arg(method),
    init=init,
    draws=draws,
    iter=iter,
    inner_iter=inner_iter,
    cond_iter=structure(cond_iter, class="integer"),
    eta=eta,
    tol=tol,
    seed=seed
  )
  obj$run()
  return(obj)
}

#' Approximate gradient-based marginal optimization
#'
#' Maximize a lower bound to the marginal posterior with respect to
#' specified parameters, with nuisance parameters marginalized out.
#'
#' @param file
#'   A character string or a connection that R supports specifying the
#'   Stan model specification in Stan's modeling language.
#' @param local_file
#'   A character string or a connection that R supports specifying the
#'   Stan model specification in Stan's modeling language.
#' @param full_model
#'   If provided, an object of class 'stanfit' that makes it unnecessary
#'   to pass 'file' or 'local_file'
#' @param data
#'   A named ‘list’ or ‘environment’ providing the data for the model
#'   or a character vector for all the names of objects used as data.
#'   See \code{\link[rstan]{stan}} for more details.
#' @param method
#'   A character string naming the conditional inference:
#'   "laplace", "vb", "sampling"
#' @param init
#'   A numeric vector of length the number of hyperparameters.
#' @param draws
#'   A positive integer, number of draws to calculate stochastic gradient.
#' @param iter
#'   A positive integer, the maximum number of iterations.
#' @param inner_iter
#'   A positive integer, the number of iterations after each
#'   conditional inference.
#' @param cond_iter
#'   A positive integer, the maximum number of iterations for the
#'   conditional inference. Default is to run until convergence.
#' @param eta
#'   Double, constant scale factor for learning rate.
#' @param tol
#'   Double, tolerance for signaling convergence.
#' @param seed
#'   The seed, a positive integer, for random number generation of
#'   Stan. The default is generated from 1 to the maximum integer
#'   supported by R so fixing the seed of R's random number generator
#'   can essentially fix the seed of Stan.  When multiple chains are
#'   used, only one seed is needed, with other chains' seeds being
#'   generated from the first chain's seed to prevent dependency among
#'   the random number streams for the chains. When a seed is
#'   specified by a number, ‘as.integer’ will be applied to it.  If
#'   ‘as.integer’ produces ‘NA’, the seed is generated randomly.  We
#'   can also specify a seed using a character string of digits, such
#'   as ‘"12345"’, which is converted to integer.
#'
#' @return
#' An object of reference class \code{"gmo"}. It is a list containing
#' the following components:
#' \item{par}{a vector of optimized parameters}
#' \item{cov}{estimated covariance matrix at \code{par}}
#' \item{sims}{\code{draws * inner_iter} many samples from the last
#'     approximation to the conditional posterior, p(alpha | y, phi)}
#'
#' @import methods
#' @importFrom rstan stan optimizing vb sampling constrain_pars log_prob grad_log_prob stan_model
#' @export
gmo_approx <- function(file, local_file, full_model, data,
  method=c("laplace", "vb", "sampling"), init="random",
  draws=5L, iter=100L, inner_iter=10L, cond_iter=NA, eta=1,
  tol=1e-3, seed=1234L) {
  if (missing(full_model)) {
    full_model <- suppressMessages(
                  stan(file, data = c(data, list(GMO_FLAG = FALSE, fixed_phi = double())),
                       chains = 0, iter = 1))
  }
  else stopifnot(is(full_model, "stanfit"))
  local_model <- if (!missing(local_file)) stan_model(local_file) else
                 get_stanmodel(full_model)
  obj <- GMO$new(
    calc_log_p="approx",
    full_model=full_model,
    local_model=local_model,
    data=data,
    method=match.arg(method),
    init=init,
    draws=draws,
    iter=iter,
    inner_iter=inner_iter,
    cond_iter=structure(cond_iter, class="integer"),
    eta=eta,
    tol=tol,
    seed=seed
  )
  obj$run()
  return(obj)
}

GMO <- setRefClass("gmo",
  fields=list(
    par="numeric",
    cov="matrix",
    sims="matrix",
    g_alpha="stanfit",
    alpha="list",
    num_par="integer",
    num_alpha="integer",
    full_model="stanfit",
    two_models="logical",
    local_model="stanmodel",
    data="ANY",
    method="character",
    draws="integer",
    iter="integer",
    eval_iter="integer",  # number of iterations to converge
    inner_iter="integer",
    cond_iter="integer",
    eta="numeric",
    tol="numeric",
    seed="integer",
    .cond_infer="function",
    .sample="function",
    .calc_log_p="function",
    .log_p="function"
  ), methods=list(
    initialize = function(calc_log_p, full_model, local_model, data, method,
      init, draws, iter, inner_iter, cond_iter, eta, tol, seed) {
      if (identical(init, "random")) {
        # Note that in stan data, par must be the parameter phi.
        # We could generalize this to count all parameters which
        # are in full_model but not the local model, although this
        # requires running g_alpha.
        num_par <<- as.integer(full_model@par_dims$phi)
        par <<- rnorm(num_par, 0, 0.01)
      } else {
        par <<- init
        num_par <<- length(init)
      }
      alpha <<- structure("random", class="list")
      full_model <<- full_model
      local_model <<- local_model
      two_models <<- !identical(local_model, get_stanmodel(full_model))
      data <<- c(data, list(GMO_FLAG = TRUE))
      method <<- method
      draws <<- draws
      iter <<- iter
      inner_iter <<- inner_iter
      cond_iter <<- cond_iter
      eta <<- eta
      tol <<- tol
      seed <<- seed

      # Initialize functions that depend on static arguments.
      # This avoids having to run an if-else chain inside the function
      # itself, which would be called at each iteration that the
      # function is called.
      #
      # When implementing GMO in Stan's C++, this would be done
      # automatically using metaprogramming tricks.
      if (method == "laplace") {
        if (is.na(cond_iter)) {
          .cond_infer <<- function(data) {
            sink(file="/dev/null", type=c("output", "message"))
            fit <- optimizing(local_model, data=data,
                              seed=seed, init=alpha,
                              as_vector=FALSE,
                              draws=inner_iter*draws, constrained=FALSE)
            closeAllConnections()
            g_alpha <<- structure(fit, class="stanfit")
            alpha <<- g_alpha$par
          }
        } else {
          .cond_infer <<- function(data) {
            sink(file="/dev/null", type=c("output", "message"))
            fit <- optimizing(local_model, data=data,
                              seed=seed, init=alpha,
                              as_vector=FALSE,
                              draws=inner_iter*draws, constrained=FALSE,
                              iter=cond_iter)
            closeAllConnections()
            g_alpha <<- structure(fit, class="stanfit")
            alpha <<- g_alpha$par
          }
        }
        .sample <<- function(m) {
          alpha_sims <- g_alpha$theta_tilde[(m-1)*draws + 1:draws, ]
          if (draws == 1) {
            alpha_sims <- matrix(alpha_sims, nrow=1)
          } else if (class(alpha_sims) != "matrix") {
            # note the latter case can occur if inner_iter=1 and we
            # end up collapsing a draws x 1 matrix into a draws vector
            alpha_sims <- matrix(alpha_sims, ncol=1)
          }
          return(alpha_sims)
        }
        .log_p <<- function(alpha_sims, m, g_flag) {
          return(.log_p_laplace(alpha_sims, m, g_flag))
        }
      } else if (method =="vb") {
        if (is.na(cond_iter)) {
          .cond_infer <<- function(data) {
            sink(file="/dev/null", type=c("output", "message"))
            g_alpha <<- vb(local_model, data=data,
                           seed=seed, init=alpha,
                           output_samples=inner_iter*draws)
            closeAllConnections()
            alpha <<- structure("random", class="list")
            #alpha <- alpha_mean(g_alpha)
          }
        } else {
          .cond_infer <<- function(data) {
            sink(file="/dev/null", type=c("output", "message"))
            g_alpha <<- vb(local_model, data=data,
                           seed=seed, init=alpha,
                           output_samples=inner_iter*draws,
                           iter=cond_iter)
                           #adapt_engaged=FALSE, # TODO use only the first time
            closeAllConnections()
            alpha <<- structure("random", class="list")
            #alpha <- alpha_mean(g_alpha)
          }
        }
        .sample <<- function(m) {
          # TODO avoid if-else chain within the function
          if (length(num_alpha) == 0) {
            pars <- get_stan_params(g_alpha)
            num_alpha <<- count_params(g_alpha, pars)
          }
          alpha_sims <- matrix(
                          unlist(attributes(g_alpha)$sim$samples[[1]][1:num_alpha]),
                          ncol=num_alpha)[(m-1)*draws + 1:draws, ]
          return(alpha_sims)
        }
        .log_p <<- function(alpha_sims, m, g_flag) {
          return(.log_p_vb(alpha_sims, m, g_flag))
        }
      } else if (method == "sampling") {
        .cond_infer <<- function(data) {
        # For sampling, cond_iter is always equal to 2*inner_iter*draws.
          sink(file="/dev/null", type=c("output", "message"))
          g_alpha <<- sampling(local_model, data=data,
                               iter=2*inner_iter*draws, chains=1,
                               seed=seed, init=alpha)
          closeAllConnections()
          alpha <<- structure("random", class="list")
          #alpha <<- alpha_mean(g_alpha)
        }
        .sample <<- function(m) {
          if (length(num_alpha) == 0) {
            pars <- get_stan_params(g_alpha)
            num_alpha <<- count_params(g_alpha, pars)
          }
          alpha_sims <- extract(g_alpha,
                                permuted=FALSE)[(m-1)*draws + 1:draws,,1:num_alpha]
          if (draws == 1) {
            alpha_sims <- matrix(alpha_sims, nrow=1)
          }
          return(alpha_sims)
        }
        .log_p <<- function(alpha_sims, m, g_flag) {
          return(.log_p_sampling(alpha_sims, m, g_flag))
        }
      } else {
        stop("Conditional inference method not valid.")
      }
      if (calc_log_p == "exact") {
        .calc_log_p <<- function(alpha_sims, m) {
          density_sims <- .log_p(alpha_sims, m, g_flag=TRUE)
          if (draws < 25) {
            log_r <- density_sims$log_p - density_sims$log_g
          } else {
            log_r <- psislw(density_sims$log_p - density_sims$log_g)$lw_smooth
          }
          max_log_r <- max(log_r)
          r <- exp(log_r - max_log_r)
          # Note that weighted.mean normalizes the importance ratios
          return(list(fn=max_log_r + log(mean(r)),
                      grad=apply(density_sims$grad_log_p, 2, weighted.mean, r)))
        }
      } else {
        .calc_log_p <<- function(alpha_sims, m) {
          # This outputs only the energy; we drop the entropy term for
          # computational reasons and because we assess the density
          # value only for convergence diagnostics.
          density_sims <- .log_p(alpha_sims, m, g_flag=FALSE)
          return(list(fn=mean(density_sims$log_p),
                      grad=apply(density_sims$grad_log_p, 2, mean)))
        }
      }
    },
    run = function() {
      diagnostic <- Diagnostic$new(tol)
      opt <- Opt$new(eta)
      for (tee in 1:iter) {
        print(sprintf("Iteration: %s", tee))
        print(par)

        if (two_models) .cond_infer(c(data, list(phi=par)))
        else .cond_infer(c(data, list(fixed_phi=par)))
        for (m in 1:inner_iter) {
          alpha_sims <- .sample(m)
          log_p <- .calc_log_p(alpha_sims, m)
          par <<- opt$update_params(par, log_p$grad, (tee-1)*inner_iter + m)
        }

        flags <- diagnostic$check_converge(par, log_p$fn, log_p$grad)
        if (any(flags)) {
          print("Optimization terminated normally:")
          print(.get_code_string(flags))
          eval_iter <<- tee
          cov <<- est_covariance(par)
          sims <<- .collect_alpha_sims()
          return()
        }
      }
      print("Maximum number of iterations hit, may not be at an optima")
      eval_iter <<- iter
      cov <<- est_covariance(par)
      sims <<- .collect_alpha_sims()
      return()
    },
    .collect_alpha_sims = function() {
      alpha_sims <- .sample(1)
      if (inner_iter > 1) {
        for (m in 2:inner_iter) {
          alpha_sims <- rbind(alpha_sims, .sample(m))
        }
      }
      return(alpha_sims)
    },
    .log_p_laplace = function(alpha_sims, m, g_flag=TRUE) {
      grad_log_p_sims <- array(NA, c(draws, num_par))
      sink(file="/dev/null")
      for (s in 1:draws) {
        grad_log_p_sims[s, ] <- grad_log_prob(full_model,
                                        c(par, alpha_sims[s, ]),
                                        adjust_transform=FALSE)[1:num_par]
      }
      closeAllConnections()
      log_p_sims <- g_alpha$log_p[(m-1)*draws + 1:draws]
      if (g_flag) {
        log_g_sims <- g_alpha$log_g[(m-1)*draws + 1:draws]
      } else {
        log_g_sims <- NA
      }
      return(list(grad_log_p=grad_log_p_sims,
                  log_p=log_p_sims,
                  log_g=log_g_sims))
    },
    .log_p_vb = function(alpha_sims, m, g_flag=TRUE) {
      # vb()/sampling() cannot output unconstrained samples via
      # constrained=FALSE. Further, to unconstrain_pars() requires a
      # list of each alpha's parameter name and data structure.
      #
      # To address this, we constrain phi, calculate log_prob on
      # constrained space without adjusting ("log_prob") and with
      # adjusting ("log_prob_Jt").
      #
      # Then
      # 2*"log_prob" - "log_prob_Jt" = "log_prob" - Jt
      # is log_prob on the constrained space adjusting to the
      # minus log determinant of the Jacobian. This is the same as
      # log_prob on the unconstrained space.
      par_const <- constrain_pars(full_model,
                                  upars=c(par,
                                          rep(0, ncol(alpha_sims))))[[1]]
      grad_log_p_sims <- array(NA, c(draws, num_par))
      log_p_sims <- rep(NA, draws)
      sink(file="/dev/null")
      for (s in 1:draws) {
        grad_log_p_sims[s, ] <- 2*grad_log_prob(full_model,
                                          c(par_const, alpha_sims[s, ]),
                                          adjust_transform=FALSE)[1:num_par] -
                            grad_log_prob(full_model,
                                          c(par_const, alpha_sims[s, ]),
                                          adjust_transform=TRUE)[1:num_par]
        log_p_sims[s] <- 2*log_prob(full_model,
                                    c(par_const, alpha_sims[s, ]),
                                    adjust_transform=FALSE) -
                           log_prob(full_model,
                                    c(par_const, alpha_sims[s, ]),
                                    adjust_transform=TRUE)
      }
      closeAllConnections()
      if (g_flag) {
        # TODO we require (preferably unconstrained)
        # log_g_sims = log_q outputted via variational.log_q
        # See aed4c42c4beebed0d43d1e1b868ef795a760fe9e in stan-dev/stan.
        # Ideally, we could also get (unconstrained)
        # log_p_sims = lp__ outputted via model_.template log_prob<false, false>
        # See ba4acfad356a85d2ebde79d82f6bcdced7034ca7 in stan-dev/stan.
        log_g_sims <- attributes(g_alpha)$sim$samples[[1]]$lp__
      } else {
        log_g_sims <- NA
      }
      return(list(grad_log_p=grad_log_p_sims,
                  log_p=log_p_sims,
                  log_g=log_g_sims))
    },
    .log_p_sampling = function(alpha_sims, m, g_flag=TRUE) {
      return(.log_p_vb(alpha_sims, m, g_flag))
      ## TODO
      ## using the awful harmonic mean estimator
      #grad_log_p_sims <- array(NA, c(draws, num_par))
      #log_p <- rep(NA, draws)
      #for (s in 1:draws) {
      #  grad_log_p_sims[s, ] <- grad_log_prob(full_model,
      #                                  c(par, alpha_sims[s, ]),
      #                                  adjust_transform=FALSE)[1:num_par]
      #  log_p_sims[s] <- log_prob(full_model,
      #                            c(par, alpha_sims[s, ]),
      #                            adjust_transform=FALSE)
      #}
    },
    .get_code_string = function(flags) {
      if (flags[1]) {
        return("  Convergence detected: absolute parameter change was below tolerance")
      } else if (flags[2]) {
        return("  Convergence detected: absolute change in objective function was below tolerance")
      } else if (flags[3]) {
        return("  Convergence detected: relative change in objective function was below tolerance")
      } else if (flags[4]) {
        return("  Convergence detected: gradient norm is below tolerance")
      } else {
        return("  Convergence detected: relative gradient magnitude is below tolerance")
      }
    }
  )
)
