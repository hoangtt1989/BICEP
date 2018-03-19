RX_mean <- function(X, beta) {
  sweep(X, 2, beta)
}

fun_reparam_mean <- function(z_hat, F_reparam, s_hat_reparam, X, wts, gamma, dual_step) {
  p <- ncol(X)
  n <- nrow(X)
  one_mat <- matrix(1, nrow = n, ncol = p)
  reparam_term <- as.vector(F_reparam %*% z_hat + s_hat_reparam)
  inner_term <- crossprod(RX_mean(X, reparam_term), wts)
  linear_term <- crossprod(gamma, inner_term)
  quad_term <- .5 * dual_step * sum( (inner_term)^2 )
  return(as.numeric(linear_term + quad_term))
}

grad_reparam_mean <- function(z_hat, F_reparam, s_hat_reparam, X, wts, gamma, dual_step) {
  p <- ncol(X)
  n <- nrow(X)
  one_mat <- matrix(1, nrow = n, ncol = p)
  reparam_term <- as.vector(F_reparam %*% z_hat + s_hat_reparam)
  linear_term <- -crossprod(F_reparam, crossprod(one_mat, wts) * gamma)
  quad_term <- -dual_step * crossprod(F_reparam, crossprod(one_mat, wts) * ((t(X) - reparam_term * t(one_mat)) %*% wts))
  return(as.numeric(linear_term + quad_term))
}

fun_mean <- function(beta, X, wts, gamma, dual_step) {
  inner_term <- crossprod(RX_mean(X, beta), wts)
  linear_term <- crossprod(gamma, inner_term)
  quad_term <- .5 * dual_step * sum( (inner_term)^2 )
  return(as.numeric(linear_term + quad_term))
}

grad_mean <- function(beta, X, wts, gamma, dual_step) {
  p <- ncol(X)
  n <- nrow(X)
  one_mat <- matrix(1, nrow = n, ncol = p)
  linear_term <- -crossprod(one_mat, wts) * gamma
  beta <- as.vector(beta)
  quad_term <- -dual_step * ((t(X) - beta * t(one_mat)) %*% wts) * crossprod(one_mat, wts)
  return(as.numeric(linear_term + quad_term))
}

mean_opt_wrapper <- function(beta_new, wts_new, gamma_curr, X, F_reparam, s_hat_reparam, dual_step, beta_opt_type, prox_fun, prox_fun_arg, prox_optim_arg) {
  if (beta_opt_type == 'LBFGS') {
    z_hat_new <- as.vector(crossprod(F_reparam, beta_new - s_hat_reparam))
    beta_res <- stats::optim(z_hat_new, fun_reparam_mean, grad_reparam_mean, F_reparam, s_hat_reparam, X, wts_new, gamma_curr, dual_step, method = 'BFGS')
    z_hat_new <- beta_res$par
    beta_new <- as.vector(F_reparam %*% z_hat_new + s_hat_reparam)
    converged <- beta_res$convergence == 0
  }
  else if (beta_opt_type == 'proximal_gradient') {
    prox_res <- do.call('proximal_gradient', c(list(fun_mean, grad_mean, beta_new, X = X, wts = wts_new, gamma = gamma_curr, dual_step = dual_step, prox_fun = prox_fun, prox_fun_arg = prox_fun_arg), prox_optim_arg), quote = T)
    beta_new <- prox_res$par
    converged <- prox_res$converged
  } else {
    stop("invalid beta optimization type")
  }
  RX_new <- RX_mean(X, beta_new)
  return(list(beta_new = beta_new, RX_new = RX_new, converged = converged))
}
