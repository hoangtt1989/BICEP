#' Profile a confidence interval for a one-dimensional mean using either REL or EL
#' 
#' @param data_mat a matrix of the data with observations in rows and variables in columns.
#' @param q an integer specifying the number of detected outliers. Default is \code{floor(.1 * n)}.
#' @param conf_level the confidence level.
#' @param test_thresh a string specifying the calibration type. Only the chi-squared approximation is currently implemented.
#' @param upper_break the initial guess for the upper limit for root finding. Default uses a value based on the magnitude of the MLE.
#' @param lower_break the initial guess for the lower limit for root finding.
#' @param upper_divide divide the MLE by this much for \code{upper_break}.
#' @param upper_increase increase the initial guess for the upper limit by this amount.
#' @param left_algorithm the algorithm for the left confidence limit.
#' @param right_algorithm the algorithm for the right confidence limit.
#' @param ... additional arguments passed to \code{REL_fixed}.
#' @param verbose a boolean to allow console output.
#' @return A list with confidence interval attributes.
#' @export
REL_mean_CI <- function(data_mat, q = NULL, conf_level = .05, test_thresh = 'chisq', 
                        upper_break = NULL, lower_break = 1e-6, 
                        upper_divide = 5, upper_increase = 1e-1, left_algorithm = c('REL', 'EL'), right_algorithm = c('REL', 'EL'), ..., verbose = F) {
  left_algorithm <- left_algorithm[1]
  right_algorithm <- right_algorithm[1]
  ##get critical value
  if (test_thresh == 'chisq') {
    test_thresh <- stats::qchisq(1 - conf_level, df = 1)
  }
  ##
  ##get mean value
  MLE <- mean(data_mat)
  ##
  ##make a function for root finding
  left_root_fun <- function(inp) {
    test_val <- MLE + inp
    Z_mat <- data_mat - test_val
    ret <- switch(left_algorithm,
                  REL = REL_fixed(Z_mat, q = q, ...),
                  EL = emplik_concord(Z_mat))
    ret$logelr_stat - test_thresh
  }
  right_root_fun <- function(inp) {
    test_val <- MLE + inp
    Z_mat <- data_mat - test_val
    ret <- switch(right_algorithm,
                  REL = REL_fixed(Z_mat, q = q, ...),
                  EL = emplik_concord(Z_mat))
    ret$logelr_stat - test_thresh
  }
  ##
  #####make sure there is a sign change in the interval
  lower_init_fit <- right_root_fun(lower_break)
  if (lower_init_fit > 0) {
    stop('lower_break must produce a test statistic below the critical value - decrease it')
  }
  if (!is.null(upper_break)) {
    upper_init_fit <- right_root_fun(upper_break)
    if (upper_init_fit < 0) {
      stop('upper_break must produce a test statistic above the critical value - increase it')
    }
  } else { ## user did not supply upper_break - we will try to find one using the magnitude of beta
    upper_break <- abs(MLE) / upper_divide
    upper_init_fit <- right_root_fun(upper_break)
    while (upper_init_fit < 1e-3) {
      if (verbose) {
        message('increasing upper_break to produce a test statistic above the critical value')
      }
      upper_break <- upper_break + upper_increase
      upper_init_fit <- right_root_fun(upper_break)
    }
  } 
  #####
  ###use root finding to get the interval
  ##upper
  start_time <- Sys.time()
  upper_int <- stats::uniroot(right_root_fun, lower = lower_break, upper = upper_break, check.conv = F)
  end_time <- Sys.time()
  upper_time <- end_time - start_time
  ##
  ##lower
  ###make sure there is a sign change
  upper_break <- -upper_int$root
  while(left_root_fun(upper_break) < 0) {
    if (verbose) {
      message('decreasing upper_break to produce a test statistic above the critical value (for lower interval)')
    }
    upper_break <- upper_break - upper_increase
  }
  ###
  start_time <- Sys.time()
  lower_int <- stats::uniroot(left_root_fun, lower = upper_break, upper = -lower_break, check.conv = F)
  end_time <- Sys.time()
  lower_time <- end_time - start_time
  ##
  ###
  return(list(upper_int = MLE + upper_int$root, lower_int = MLE + lower_int$root, mean_val = MLE, time = as.double(upper_time + lower_time, unit = 'secs')))
}
##################


#' Compute REL for a fixed estimate.
#' 
#' @param Z_mat a matrix of the score equations with observations in rows and variables in columns.
#' @param q an integer specifying the number of detected outliers. Default is \code{floor(.1 * n)}.
#' @param outlier_loop_arg control arguments passed to the \code{delta} optimization; see \code{\link{optim_outlier_control}}.
#' @param wts_init an initial guess for the weights vector. Default is \code{1/n}.
#' @param vary_penalty a string specifying the type of penalty variation.
#' @param reuse_delta use the last estimate of \code{delta} in the current outer loop iteration.
#' @param RB_mu the \code{mu} parameter for the residual balancing scheme.
#' @param RB_tau the \code{tau} parameter for the residual balancing scheme.
#' @param RB2_ksi the \code{ksi} parameter for the \code{RB2} scheme.
#' @param outer_eps absolute tolerance required for outer loop convergence.
#' @param outer_rel_eps relative tolerance required for outer loop convergence.
#' @param dual_step initial penalty parameter.
#' @param outer_maxit maximum number of outer loop iterations.
#' @param wts_beta_rep the number of repetitions of the block coordinate descent update within each outer loop iteration.
#' @param outer_tol_type the type of tolerance check for the outer loop.
#' @param mirror_arg control arguments passed to the mirror descent optimization; see \code{\link{optim_mirror_control}}.
#' @param verbose a boolean to allow console output.
#' @return A list with estimates and convergence information.
#' @export
REL_fixed <- function(Z_mat, q = NULL,
                           outlier_loop_arg = list(), wts_init = NULL, 
                           vary_penalty = c('RB', 'RB2', 'none'), reuse_delta = F,
                           RB_mu = 10, RB_tau = 2, RB2_ksi = 2, outer_eps = 1e-8, outer_rel_eps = 1e-4, dual_step = 2, 
                           outer_maxit = 1000, wts_beta_rep = 1,
                           outer_tol_type = c('primal_dual', 'fval', 'gamma', 'primal'), 
                           mirror_arg = optim_mirror_control(line_search = T), verbose = F) {
  ##checking inputs
  outer_tol_type <- outer_tol_type[1]
  vary_penalty <- vary_penalty[1]
  primal_dual_flag <- ifelse(vary_penalty %in% c('RB', 'RB2') || outer_tol_type %in% c('primal', 'primal_dual'), T, F)
  gamma_diff_flag <- ifelse(outer_tol_type == 'gamma', T, F)
  if (!(outer_tol_type %in% c('gamma', 'fval', 'primal', 'primal_dual'))) {
    stop('supply a valid outer tolerance type')
  }
  if (!(vary_penalty %in% c('RB', 'RB2', 'none'))) {
    stop('supply a valid penalty varying type')
  }
  ##
  ##initial parameters
  tiny <- .Machine$double.eps
  tiny2 <- tiny^(.25)
  tiny_med <- tiny^(.5)
  prob_size <- dim(Z_mat)
  n <- prob_size[1]
  if (is.null(q)) {
    q <- floor(.1 * n)
  }
  q <- max(0, min(n, floor(q)))
  p <- prob_size[2]
  if (is.null(wts_init)) {
    wts_init <- rep(1, n)
  } else {
    if (length(wts_init) != n) {
      stop('wts_init must be a vector with n elements')
    }
  }
  wts_new <- wts_init
  delta_init <- rep(0, n)
  delta_curr <- delta_init
  delta_new <- delta_curr
  wts_curr <- wts_new
  gamma_curr <- rep(0, p)
  outer_converged <- F
  ##
  ##primal_dual stopping
  if (outer_tol_type == 'primal_dual') {
    primal_eps <- outer_eps * sqrt(p)
  }
  # dual_eps <- outer_eps * sqrt(n) ## dimension of dual residual changes depending on update order
  ##
  ##initialize holders
  outer_fval <- rep(NA, outer_maxit)
  mirror_nits <- rep(NA, outer_maxit)
  beta_outlier_nits <- rep(NA, outer_maxit)
  beta_outlier_tol <- rep(NA, outer_maxit)
  rho_vec <- rep(NA, outer_maxit)
  ##
  ###initialize function/gradient arguments
  fun_grad_arg <- list(gamma = gamma_curr, Z_mat = Z_mat, dual_step = dual_step)
  ##
  outer_fval_curr <- sum(log(wts_new))
  ####
  #################main loop
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    # print(paste('Outer loop', j))
    if (verbose) {
      if (! j %% 50) {
        message(paste('Iteration', j))
      }
    }
    ##random ordering of wts/beta update
    dual_resid_elt <- sqrt(n) ##dual residual has n elements
    if (outer_tol_type == 'primal_dual') {
      dual_eps <- outer_eps * dual_resid_elt
    }
    ##
    ##wts and beta updates can be repeated
    for (wts_beta_iter in 1:wts_beta_rep) {
      #######mirror descent
      A_max <- n * dual_step * max(abs(tcrossprod(Z_mat))) ##t(Xtr_R) %*% Xtr_R
      mirror_step <- 1 / ((A_max + max(abs(Z_mat %*% gamma_curr))) * max(abs(1 - delta_new)))
      mirror_res <- do.call('mirror_descent_REL', c(list(wts = wts_new, delta = delta_new, mirror_step = mirror_step), fun_grad_arg, mirror_arg), quote = T)
      wts_new <- as.vector(mirror_res$wts)
      #######
      #######outlier update
      if (reuse_delta) {
        delta_inp <- delta_new
      } else {
        delta_inp <- delta_init
      }
      outlier_res <- do.call('prox_grad_REL', c(list(wts_new, Z_mat, gamma_curr, dual_step, delta_inp, q), outlier_loop_arg), quote = T)
      delta_new <- outlier_res$delta_opt
      #######
    }
    ##
    mirror_nits[j] <- mirror_res$iter
    beta_outlier_nits[j] <- outlier_res$iter
    beta_outlier_tol[j] <- outlier_res$tol
    if (verbose) {
      if (!mirror_res$converged) {
        warning(paste('Mirror descent did not converge at iteration', j))
      }
      if (!outlier_res$converged) {
        warning(paste('REL update did not converge at iteration', j))
      }
    }
    delta_diff <- (delta_curr - delta_new)
    #######gamma update
    gamma_new <- gamma_curr + dual_step * crossprod(Z_mat, wts_new * (1 - delta_new))
    
    if (gamma_diff_flag) {
      gamma_sc <- max(abs(gamma_curr))
      gamma_sc <- ifelse(gamma_sc < .Machine$double.eps, .Machine$double.eps^(.25), gamma_sc)
      gamma_diff <- max(abs(gamma_new - gamma_curr)) / gamma_sc
    }
    #######
    #######dual residual
    if (primal_dual_flag) {
      dual_resid <- sqrt( sum( (dual_step * (tcrossprod(Z_mat) %*% (wts_new * (delta_diff))) * (delta_curr - 1))^2 ) )
      dual_resid <- ifelse(dual_resid < tiny, tiny_med, dual_resid)
      dual_resid_scale <- sqrt( sum( ((Z_mat %*% gamma_new) * (1 - delta_curr))^2 ) ) 
      dual_resid_scale <- ifelse(dual_resid_scale < tiny, tiny_med, dual_resid)
    }
    #######
    #######primal residual
    if (primal_dual_flag) {
      primal_resid <- sqrt( sum( crossprod(Z_mat, wts_new * (1 - delta_new))^2 ) )
      primal_resid_scale <- 1
    }
    #######
    ##update delta, gamma, wts
    delta_curr <- as.vector(delta_new)
    wts_curr <- as.vector(wts_new)
    gamma_curr <- gamma_new
    ##
    ####update fun_arg and grad_arg
    fun_grad_arg$gamma <- gamma_curr
    ####
    ###function value update
    outer_fval_new <- sum(log(wts_curr))
    outer_fval[j] <- outer_fval_new
    if (outer_tol_type == 'fval') {
      outer_fval_sc <- max(abs(outer_fval_curr))
      outer_fval_sc <- ifelse(outer_fval_sc < tiny, tiny2, outer_fval_sc)
      outer_fval_diff <- abs(outer_fval_new - outer_fval_curr) / abs(outer_fval_sc)
    }
    outer_fval_curr <- outer_fval_new
    ###
    #######vary penalty parameter
    if (vary_penalty == 'RB') {
      dual_step <- RB_vary(dual_step, primal_resid, dual_resid, RB_mu, RB_tau)
    }
    if (vary_penalty == 'RB2') {
      dual_step <- RB2_vary(dual_step, primal_resid, primal_resid_scale, dual_resid, dual_resid_scale, RB2_ksi, RB_tau, RB_mu)
    }
    fun_grad_arg$dual_step <- dual_step
    rho_vec[j] <- dual_step
    #######
    #######stopping
    if (outer_tol_type == 'gamma') {
      outer_tol <- gamma_diff
    } else if (outer_tol_type == 'fval') {
      outer_tol <- outer_fval_diff
    } else if (outer_tol_type %in% c('primal', 'primal_dual')) {
      outer_tol <- primal_resid
    }
    if (outer_tol_type == 'primal_dual') {
      primal_eps_scale <- primal_eps + outer_rel_eps * primal_resid_scale / n
      dual_eps_scale <- dual_eps + outer_rel_eps * dual_resid_scale
      if (primal_resid / n < primal_eps_scale && dual_resid / n < dual_eps_scale) {
        outer_converged <- T
        break
      }
    } else {
      if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
        outer_converged <- T
        break
      }
    }
    #######
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  #################
  if (!outer_converged) {
    warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
  }
  
  ####cleaning up output
  rho_vec <- rho_vec[1:j]
  outer_fval <- outer_fval[1:j]
  mirror_nits <- mirror_nits[1:j]
  beta_outlier_nits <- beta_outlier_nits[1:j]
  beta_outlier_tol <- beta_outlier_tol[1:j]
  logelr <- outer_fval_curr
  if (!primal_dual_flag) {
    primal_resid <- 'not calculated'
    dual_resid <- 'not calculated'
  }
  ####
  
  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
              gamma = gamma_curr, delta_opt = delta_new, outlier_idx = which(delta_new != 0), 
              outer_converged = outer_converged, time = as.double(tot_time, 'secs'), 
              outer_fval = -outer_fval, primal_resid = primal_resid, dual_resid = dual_resid, rho = rho_vec,
              outer_tol = outer_tol, outlier_nits = beta_outlier_nits, outlier_tol = beta_outlier_tol,
              mirror_nits = mirror_nits, outer_nits = j))
}
##################

