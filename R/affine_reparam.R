##################affine reparam helper
#' Construct the variables for the affine reparameterization
#' 
#' @param mean_init a vector of the MLE.
#' @param fix_index a vector containing the indices to be profiled.
#' @param alpha_add a vector of the same length as \code{fix_index} containing the adjustments away from the MLE.
#' @return A list with the new candidate value for \code{beta} and the matrices/vectors required for the reparameterization.
#' @export
reparam_helper <- function(mean_init, fix_index, alpha_add) {
  p <- length(mean_init)
  ##set up the reparam
  T_reparam <- diag(p)
  T_reparam <- T_reparam[fix_index, ]
  if (length(fix_index) == 1) {
    T_reparam <- t(as.matrix(T_reparam))
  }
  alpha_test <- mean_init[fix_index] + alpha_add
  s_hat_reparam <- rep(0, p)
  s_hat_reparam[fix_index] <- alpha_test
  F_reparam <- pracma::orth(pracma::nullspace(T_reparam))
  ##
  beta_test <- mean_init
  beta_test[fix_index] <- alpha_test
  return(list(beta_test = beta_test, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, alpha_test = alpha_test, T_reparam = T_reparam))
}
##################
