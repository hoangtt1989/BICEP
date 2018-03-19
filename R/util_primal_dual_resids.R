primal_resid_calc <- function(Z_mat_new, wts_new) {
  sqrt( sum( crossprod(Z_mat_new, wts_new)^2 ) )
}

canonical_primal_resid_sc <- function(X, y, wts_new, beta_new, link_fun) {
  wts_X <- wts_new * X
  max(sqrt(sum((crossprod(wts_X, y))^2)), sqrt(sum((crossprod(wts_X, link_fun(X %*% beta_new)))^2)), 1)
}

mean_primal_resid_sc <- function(X, wts_new, beta_new) {
  p <- length(beta_new)
  n <- length(wts_new)
  ones_mat <- matrix(1, nrow = n, ncol = p)
  beta_new <- as.vector(beta_new)
  max( sqrt( sum( crossprod(X, wts_new)^2 ) ), sqrt( sum( ((beta_new * t(ones_mat)) %*% wts_new)^2 ) ), 1)
}

REL_wts_dual_resid <- function(delta_new, delta_curr, wts_new, Z_mat_curr, dual_step) {
  sqrt( sum( (-dual_step * (tcrossprod(Z_mat_curr) %*% (wts_new * (delta_new - delta_curr))) * (1 - delta_curr))^2 ) )
}

REL_delta_dual_resid <- function(Z_mat_new, Z_mat_curr, wts_new, delta_new, dual_step) {
  sqrt( sum( (-dual_step * (tcrossprod(Z_mat_curr, Z_mat_curr - Z_mat_new) %*% (wts_new * (1 - delta_new))) * wts_new )^2 ) )
}

REL_dual_resid_sc <- function(Z_mat_curr, gamma_new, wts_new, delta_curr) {
  Z_gam <- Z_mat_curr %*% gamma_new
  max(sqrt( sum( ( Z_gam * wts_new )^2 ) ), sqrt( sum( ( Z_gam * (1 - delta_curr) )^2 ) ) )
}


REL_dual_resids <- function(Z_mat_new, Z_mat_curr, delta_new, delta_curr, wts_new, gamma_new, dual_step) {
  # n <- length(wts_new)
  wts_dual_resid <- REL_wts_dual_resid(delta_new, delta_curr, wts_new, Z_mat_curr, dual_step)
  delta_dual_resid <- REL_delta_dual_resid(Z_mat_new, Z_mat_curr, wts_new, delta_new, dual_step)
  max_dual_resid <- max(wts_dual_resid, delta_dual_resid)
  dual_resid_scale <- max(REL_dual_resid_sc(Z_mat_curr, gamma_new, wts_new, delta_curr), 1)
  return(list(dual_resid = max_dual_resid, dual_resid_scale = dual_resid_scale))
}

CEL_wts_dual_resids <- function(Z_mat_new, Z_mat_curr, wts_new, gamma_new, dual_step) {
  # n <- length(wts_new)
  dual_resid <- sqrt(sum( (dual_step * Z_mat_curr %*% crossprod(Z_mat_new - Z_mat_curr, wts_new))^2 ))
  dual_resid_scale <- max(sqrt(sum((Z_mat_curr %*% gamma_new)^2)), 1)
  return(list(dual_resid = dual_resid, dual_resid_scale = dual_resid_scale))
}
