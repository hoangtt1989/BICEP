glm_beta_opt <- function(beta_new, wts_new, gamma_curr, y, X, F_reparam, s_hat_reparam, dual_step, family, beta_opt_type, prox_fun, prox_fun_arg, lbfgs_arg, prox_optim_arg) {
  if (beta_opt_type == 'closed_form') { ##linear model and closed form solution
    z_hat_new <- closed_beta_lm_AL_cpp(wts_new, gamma_curr, y, X, F_reparam, s_hat_reparam, dual_step) ## use a cpp function for the closed form solution
    converged <- T
    beta_new <- as.vector(F_reparam %*% z_hat_new + s_hat_reparam)
  } else if (beta_opt_type == 'LBFGS') { ## if family is not gaussian and optimization type is not closed_form, use LBFGS to optimize beta (affine reparam)
    z_hat_curr <- as.vector(crossprod(F_reparam, beta_new - s_hat_reparam))
    lbfgs_res <- do.call('lbfgs', c(list(fun_reparam_glm_AL, grad_reparam_glm_AL, z_hat_curr, 
                                         F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, 
                                         X = X, y = y, wts = wts_new, gamma = gamma_curr, dual_step = dual_step, family = family), lbfgs_arg), quote = T)
    z_hat_new <- as.vector(lbfgs_res$par)
    converged <- lbfgs_res$convergence >= 0
    beta_new <- as.vector(F_reparam %*% z_hat_new + s_hat_reparam)
  } else if (beta_opt_type == 'proximal_gradient') { ## proximal gradient descent for other cases
    prox_res <- do.call('proximal_gradient', c(list(fun_glm_AL, grad_glm_AL, beta_new, X = X, y = y, wts = wts_new, gamma = gamma_curr, dual_step = dual_step, family = family, prox_fun = prox_fun, prox_fun_arg = prox_fun_arg), prox_optim_arg), quote = T)
    beta_new <- prox_res$par
    converged <- prox_res$converged
  } else {
    stop('invalid beta optimization type')
  }
  link_fun <- switch(family,
                     gaussian = base::identity,
                     binomial = binomial_link,
                     poisson = exp)
  RX_new <- RX_canonical(X, y, beta_new, link_fun)
  return(list(beta_new = beta_new, RX_new = RX_new, converged = converged))
}
