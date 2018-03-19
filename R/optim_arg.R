#' Control parameters for proximal gradient descent
#' 
#' @param tol_type the type of tolerance checking.
#' @param step_size initial step size.
#' @param accelerate a boolean to turn on acceleration.
#' @param line_search a boolean to turn on line search - it is NOT recommended to turn this off.
#' @param ls_beta by how much should the step-size increase or decrease in each line search iteration?
#' @param ls_eps you can miss the line search check by this much.
#' @param ls_maxit the maximum number of line search iterations.
#' @param maxit the maximum number of iterations.
#' @param eps the maximum value of the tolerance before declaring convergence.
#' @export
optim_prox_control <- function(tol_type = c('par', 'fval'), step_size = 1, accelerate = F, line_search = T, 
                               ls_beta = .75, ls_eps = 1e-10, ls_maxit = 50, maxit = 1000, eps = 1e-6) {
  return(list(tol_type = tol_type, step_size = step_size, accelerate = accelerate, line_search = line_search,
              ls_beta = ls_beta, ls_eps = ls_eps, ls_max_iter = ls_maxit, maxit = maxit, prox_eps = eps))
}

#' Control parameters for mirror descent
#' 
#' @param maxit the maximum number of iterations.
#' @param eps the maximum value of the tolerance before declaring convergence.
#' @param tol_type the type of tolerance checking.
#' @param vary_step for a given step size, a boolean indicating whether it should be decreased at every iteration.
#' @param line_search a boolean to turn on line search - it is NOT recommended to turn this off.
#' @param ls_beta by how much should the step-size increase or decrease in each line search iteration?
#' @param ls_eps you can miss the line search check by this much.
#' @param ls_maxit the maximum number of line search iterations.
#' @export
optim_mirror_control <- function(maxit = 1000, eps = 1e-8, tol_type = c('fval', 'wts', 'logelr'),
                                 vary_step = F, line_search = F, ls_beta = .75, ls_eps = 1e-10, ls_maxit = 100) {
  tol_type <- tol_type[1]
  return(list(maxit = maxit, mirror_eps = eps, tol_type = tol_type, vary_step = vary_step, line_search = line_search,
              ls_eps = ls_eps, ls_maxit = ls_maxit, ls_beta = ls_beta))
}

#' Control parameters for the proximal gradient descent update for outlier detection

#' @param maxit the maximum number of iterations.
#' @param eps the maximum value of the tolerance before declaring convergence.
#' @param accelerate a boolean to turn on acceleration.
#' @param line_search a boolean to turn on line search - it is NOT recommended to turn this off.
#' @param ls_beta by how much should the step-size increase or decrease in each line search iteration?
#' @param ls_eps you can miss the line search check by this much.
#' @param ls_maxit the maximum number of line search iterations.
#' @export
optim_outlier_control <- function(maxit = 2000, eps = 1e-5, accelerate = T, line_search = T,
                                  ls_beta = .75, ls_eps = 1e-12, ls_maxit = 100) {
  return(list(maxit = maxit, outlier_eps = eps, accelerate = accelerate, line_search = line_search,
              ls_eps = ls_eps, ls_maxit = ls_maxit, ls_beta = ls_beta))
}

#' Control parameters for the inner loop dual formulation
#' 
#' @param thresh convergence threshold.
#' @param itermax maximum number of iterations.
#' @param verbose print console output.
#' @export
optim_owen_inner_control <- function(thresh = 1e-16, itermax = 150, verbose = F) {
  return(list(thresh = thresh, itermax = itermax, verbose = verbose))
}

#' Control parameters for damped Newton
#' 
#' @param ls_beta by how much should the step-size increase or decrease in each line search iteration?
#' @param ls_eps you can miss the line search check by this much.
#' @param ls_c1 the fraction of the decrease in f predicted by linear extrapolation that we will accept.
#' @param ls_maxit the maximum number of line search iterations.
#' @param mod_hess the Hessian to be positive definite - recommended to leave this \code{FALSE}.
#' @export
optim_newton_control <- function(ls_beta = .7, ls_eps = 1e-12, ls_c1 = .25, ls_maxit = 100, mod_hess = F) {
  return(list(ls_beta = ls_beta, ls_eps = ls_eps, ls_c1 = ls_c1, ls_maxit = ls_maxit, mod_hess = mod_hess))
}
