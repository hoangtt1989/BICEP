######canonical link functions
binomial_link <- function(mean_val) {
  1 / (1 + exp(-mean_val))
}
binomial_link_deriv <- function(mean_val) {
  exp(mean_val) / ((1 + exp(mean_val))^2)
}
gaussian_link_deriv <- function(mean_val) {
  rep(1, length(mean_val))
}
###function value for reparameterized beta
fun_reparam_glm_AL <- function(z_hat, F_reparam, s_hat_reparam, X, y, wts, gamma, dual_step, family) {
  link_fun <- switch(family,
                     gaussian = base::identity,
                     binomial = binomial_link,
                     poisson = exp)
  wts_X <- (wts * X)
  temp_term <- crossprod(wts_X, y - link_fun(X %*% (F_reparam %*% z_hat + s_hat_reparam)))
  as.numeric(as.vector(gamma) %*% as.vector(temp_term) + 0.5 * dual_step * sum((temp_term)^2))
}
###gradient for reparameterized beta
grad_reparam_glm_AL <- function(z_hat, F_reparam, s_hat_reparam, X, y, wts, gamma, dual_step, family) {
  link_fun <- switch(family,
                     gaussian = base::identity,
                     binomial = binomial_link,
                     poisson = exp)
  link_deriv <- switch(family,
                       gaussian = gaussian_link_deriv,
                       binomial = binomial_link_deriv,
                       poisson = exp)
  wts_X <- (wts * X)
  inner_link_term <- X %*% (F_reparam %*% z_hat + s_hat_reparam)
  link_deriv_term <- link_deriv(inner_link_term)
  X_F <- X %*% F_reparam
  ret <- - crossprod(X_F, (wts_X %*% gamma) * link_deriv_term) - dual_step * crossprod(X_F, (tcrossprod(wts_X) %*% (y - link_fun(inner_link_term))) * link_deriv_term)
  as.numeric(ret)
}
###
###function value for beta (for proximal gradient descent)
fun_glm_AL <- function(beta, X, y, wts, gamma, dual_step, family) {
  link_fun <- switch(family,
                     gaussian = base::identity,
                     binomial = binomial_link,
                     poisson = exp)
  wts_X <- (wts * X)
  temp_term <- crossprod(wts_X, y - link_fun(X %*% beta))
  as.numeric(as.vector(gamma) %*% as.vector(temp_term) + 0.5 * dual_step * sum((temp_term)^2))
}
###
###gradient for beta (for proximal gradient descent)
grad_glm_AL <- function(beta, X, y, wts, gamma, dual_step, family) {
  link_fun <- switch(family,
                     gaussian = base::identity,
                     binomial = binomial_link,
                     poisson = exp)
  link_deriv <- switch(family,
                       gaussian = gaussian_link_deriv,
                       binomial = binomial_link_deriv,
                       poisson = exp)
  wts_X <- (wts * X)
  inner_link_term <- X %*% beta
  link_deriv_term <- link_deriv(inner_link_term)
  ret <- - crossprod(X, (wts_X %*% gamma) * link_deriv_term) - dual_step * crossprod(X, (tcrossprod(wts_X) %*% (y - link_fun(inner_link_term))) * link_deriv_term)
  as.numeric(ret)
}
###
######


RX_canonical <- function(X, y, beta, link_fun) {
  R_vec <- as.vector(y - link_fun(X %*% beta))
  RX_mat <- R_vec * X
  return(RX_mat)
}
