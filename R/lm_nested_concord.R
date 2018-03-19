##################function for CEL linear regression nested
#' Run the nested algorithm for linear regression composite tests.
#' 
#' @param beta_test a vector of candidate test values.
#' @param X a design matrix with observations in rows and variables in columns.
#' @param y a vector of the responses.
#' @param F_reparam the \code{F} matrix for the affine reparameterization.
#' @param s_hat_reparam the \code{s} vector for the affine reparameterization.
#' @param gamma_init a vector of initial values for the dual, or a \code{"random"} to specify sampling from \code{\link[stats]{rnorm}}. Default is zero vector.
#' @param outer_eps absolute tolerance required for outer loop convergence.
#' @param outer_maxit maximum number of outer loop iterations.
#' @param inner_opt_type optimization type for the inner loop.
#' @param inner_owen_arg control arguments passed to the dual formulation (\code{owen}) optimization for the inner loop; see \code{\link{optim_owen_inner_control}}.
#' @param inner_lbfgs_arg a list of arguments passed to \code{\link[lbfgs]{lbfgs}}.
#' @param beta_newton_arg control arguments passed to the damped Newton optimization for the inner loop; see \code{\link{optim_newton_control}}.
#' @param outer_tol_type the type of tolerance checking for the outer loop.
#' @param verbose a boolean to allow console output.
#' @export
CEL_lm_nested <- function(beta_test, X, y, F_reparam, s_hat_reparam, gamma_init = NULL,
                          outer_eps = 1e-8, outer_maxit = 1000, inner_opt_type = c('owen', 'LBFGS'), 
                          inner_owen_arg = optim_owen_inner_control(), inner_lbfgs_arg = list(invisible = 1), beta_newton_arg = optim_newton_control(),
                          outer_tol_type = c('fval', 'gval'), verbose = F) {
  ##checking inputs
  outer_tol_type <- outer_tol_type[1]
  inner_opt_type <- inner_opt_type[1]
  if (!(outer_tol_type %in% c('fval', 'gval'))) {
    stop('supply a valid outer tolerance type')
  }
  ##
  ##compile owen's functions
  owen_EL <- compiler::cmpfun(emplik_concord) ##inner loop
  ##
  ##compile gradient/hessian functions
  lm_hess_z <- compiler::cmpfun(lm_hess_z)
  lm_part_z <- compiler::cmpfun(lm_part_z)
  lm_fval_z <- compiler::cmpfun(lm_fval_z)
  lm_fval_gamma <- compiler::cmpfun(lm_fval_gamma)
  lm_part_gamma <- compiler::cmpfun(lm_part_gamma)
  ##
  ##initial parameters
  prob_size <- dim(X)
  n <- prob_size[1]
  p <- prob_size[2]
  if (is.null(gamma_init)) {
    gamma_init2 <- rep(0, p)
  } else if (gamma_init == 'random'){
    gamma_init2 <- stats::rnorm(p)
  } else {
    if (length(gamma_init) != p) {
      stop('gamma_init must have p elements')
    }
    gamma_init2 <- gamma_init
  }
  beta_curr <- beta_test
  z_reparam_curr <- crossprod(F_reparam, beta_test - s_hat_reparam)
  R_curr <- as.vector(y - X %*% beta_test)
  outer_converged <- F
  ##
  ##initialize holders
  outer_fval <- rep(NA, outer_maxit)
  inner_nits <- rep(NA, outer_maxit)
  ##
  ##first inner loop iteration
  score_curr <- R_curr * X
  if (inner_opt_type == 'owen') {
    inner_res <- do.call('emplik_concord', c(list(score_curr, lam = gamma_init2), inner_owen_arg), quote = T)
    gamma_curr <- inner_res$lambda
    logelr_curr <- inner_res$logelr
    # logelr_curr <- sum(mllog(inner_res$wts, 1/n, der = 0))
  } else if (inner_opt_type == 'LBFGS') {
    inner_res <- do.call('lbfgs', c(list(lm_fval_gamma, lm_part_gamma, gamma_init2, score_mat = score_curr), inner_lbfgs_arg), quote = T)
    gamma_curr <- inner_res$par
    logelr_curr <- sum(mllog(1 + score_curr %*% gamma_curr, 1/n, der = 0))
  }
  ##
  ##set up function, gradient, hessian arguments
  fun_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
  grad_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
  hess_arg <- list(X = X, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, y = y, gamma = gamma_curr)
  ##
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    # print(j)
    if (verbose) {
      if (! j %% 5) {
        message(paste('Iteration', j))
      }
    }
    ##beta optimization
    damped_res <- do.call('damped_newton_step', c(list(z_reparam_curr, lm_fval_z, lm_part_z, lm_hess_z, fun_arg, grad_arg, hess_arg), beta_newton_arg), quote = T)
    if (!damped_res$ls_converged) {
      warning(paste('Line search unsuccessful at j =', j))
    }
    x_curr <- damped_res$x_curr
    fval_curr <- damped_res$fval_curr
    ##
    ##update values
    z_reparam_curr <- x_curr
    beta_curr <- F_reparam %*% z_reparam_curr + s_hat_reparam
    R_new <- as.vector(y - X %*% beta_curr)
    score_new <- R_new * X
    outer_fval_curr <- fval_curr
    # outer_fval_new <- damped_res$fval_new
    ##
    ####inner loop
    if (!is.null(gamma_init)) {
      if (gamma_init == 'random') {
        gamma_init2 <- stats::rnorm(p)
      }
    }
    if (inner_opt_type == 'owen') {
      inner_res <- do.call('emplik_concord', c(list(score_new, lam = gamma_init2), inner_owen_arg), quote = T)
      gamma_curr <- inner_res$lambda
      inner_nits[j] <- inner_res$nits
      outer_fval_new <- -inner_res$logelr
      # outer_fval_new <- -sum(mllog(inner_res$wts, 1/n, der = 0))
    } else if (inner_opt_type == 'LBFGS') {
      inner_res <- do.call('lbfgs', c(list(lm_fval_gamma, lm_part_gamma, gamma_init2, score_mat = score_new), inner_lbfgs_arg), quote = T)
      gamma_curr <- inner_res$par
      inner_nits[j] <- inner_res$convergence
      outer_fval_new <- -sum(mllog(1 + score_new %*% gamma_curr, 1/n, der = 0))
    }
    outer_fval[j] <- outer_fval_new
    fun_arg$gamma <- gamma_curr
    grad_arg$gamma <- gamma_curr
    hess_arg$gamma <- gamma_curr
    ####
    #####stopping
    outer_tol <- switch(outer_tol_type,
                        fval = abs(outer_fval_new - outer_fval_curr)/abs(outer_fval_curr),
                        gval = sum(damped_res$gval_new^2))
    if (all(outer_tol < outer_eps) & !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
      outer_converged <- T
      break
    }
    #####
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  if (!outer_converged) {
    warning(paste('Outer loop did not converge, current tolerance is', outer_tol))
  }
  
  ####cleaning up output
  outer_fval <- outer_fval[1:j]
  inner_nits <- inner_nits[1:j]
  logelr <- -outer_fval_new
  wts <- 1 / (1 + score_new %*% gamma_curr) / n
  sum_wts <- sum(wts)
  primal_resid <- sqrt(sum((crossprod(score_new, wts))^2))
  if (abs(sum_wts - 1) > 1e-6) { ##output warning if weights don't sum to 1
    warning('weights do not sum to 1 - convex hull constraint may not be satisfied')
  }
  ####
  
  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts, sum_wts = sum_wts,
              gamma = gamma_curr, beta_opt = as.vector(beta_curr), primal_resid = primal_resid, 
              outer_converged = outer_converged, time = as.double(tot_time, units = 'secs'),
              outer_fval = outer_fval, outer_tol = outer_tol, inner_nits = inner_nits, outer_nits = j))
}
##################
