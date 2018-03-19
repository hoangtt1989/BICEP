#' Proximal operator for affine reparameterization
#' 
#' @param inp a vector.
#' @param T_reparam a matrix.
#' @param alpha_reparam a vector or double.
#' @export
prox_affine <- function(inp, T_reparam, alpha_reparam) {
  inp - crossprod(T_reparam, solve(tcrossprod(T_reparam), T_reparam %*% inp - alpha_reparam))
}
#' Proximal operator for inequality testing
#' 
#' @param inp a vector.
#' @param test_idx the index for the values to be tested.
#' @param test_val the values for the test.
#' @param dir the direction of the inequality; greater than (\code{gt}) or less than (\code{lt}).
#' @export
prox_inequality <- function(inp, test_idx, test_val, dir = c('gt', 'lt')) {
  if (!all(test_idx <= length(inp) & test_idx >= 1)) {
    stop('test_idx must be contained in the observations')
  }
  dir <- dir[1]
  if (length(test_val) != length(test_idx)) {
    test_val <- rep(test_val[1], length(test_idx))
  }
  cond_check <- inp[test_idx] >= test_val[order(test_idx)]
  if (dir == 'lt') {
    cond_check <- !cond_check
  }
  new_vals <- test_val
  new_vals[cond_check] <- inp[test_idx[cond_check]]
  ret <- inp
  ret[test_idx] <- new_vals
  return(ret)
}
#' Proximal operator for interval testing
#' 
#' @param inp a vector.
#' @param test_idx index for the values to be tested.
#' @param left_val a vector of values for the left interval. If the dimension does not agree with \code{test_idx}, the test is repeated for all indices in \code{test_idx}.
#' @param right_val the values for the right interval. If the dimension does not agree with \code{test_idx}, the test is repeated for all indices in \code{test_idx}.
#' @export
prox_interval <- function(inp, test_idx, left_val, right_val) {
  if (!(test_idx %in% 1:length(inp))) {
    stop('test_idx must be contained in the observations')
  }
  if (!all(left_val <= right_val)) {
    stop('left_val must be <= right_val')
  }
  if (length(left_val) != length(right_val)) {
    stop('left_val must have same length as right_val')
  }
  if (length(left_val) != length(test_idx)) {
    left_val <- rep(left_val[1], length(test_idx))
    right_val <- rep(right_val[1], length(test_idx))
  }
  order_idx <- order(test_idx)
  left_check <- inp[test_idx] < left_val[order_idx]
  right_check <- inp[test_idx] > right_val[order_idx]
  ret <- inp
  ret[test_idx[left_check]] <- left_val[left_check]
  ret[test_idx[right_check]] <- right_val[right_check]
  return(ret)
}
