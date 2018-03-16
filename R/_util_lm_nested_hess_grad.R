##################functions for hessian/gradient calculations
lm_fval_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X), rev_sign = F) {
  n <- nrow(X)
  beta_new <- F_reparam %*% z_reparam + s_hat_reparam
  R_new <- as.vector(y - X %*% beta_new)
  score_new <- R_new * X
  ret <- -sum(mllog(1 + score_new %*% gamma, llog_eps, der = 0))
  if (rev_sign) {
    ret <- -ret
  }
  return(ret)
}
lm_hess_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X)) {
  ##L_z_gamma
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  score_mat <- score_eq * X
  inner_sum <- 1 + score_eq * (X %*% gamma)
  temp_plog1 <- -mllog(inner_sum, llog_eps, der = 1)[, 2]
  temp_plog2 <- -mllog(inner_sum, llog_eps, der = 2)[, 3]
  temp1 <- matrix(0, nrow = p, ncol = ncol(X))
  temp2 <- temp1
  for (i in 1:n) {
    temp1 <- temp1 + temp_plog1[i] * crossprod(F_reparam, tcrossprod(X[i, ]))
    temp2 <- temp2 + temp_plog2[i] * crossprod(F_reparam, X[i, ] * as.numeric(X[i, ] %*% gamma)) %*% t(X[i, ]) * (y[i] - as.numeric(t(X[i, ]) %*% (F_reparam %*% z_reparam + s_hat_reparam)))
  }
  L_z_gamma <- -temp1 - temp2
  ##
  L_gamma_z <- t(L_z_gamma)
  L_gamma_gamma <- crossprod(score_mat, as.vector(-mllog(1 + score_mat %*% gamma, llog_eps, der = 2)[, 3]) * score_mat)
  ##L_z_z
  temp_res <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
    temp_res <- temp_res + temp_plog2[i] * (crossprod(F_reparam, X[i, ]) * as.numeric(as.numeric(gamma) %*% X[i, ])^2) %*% crossprod(X[i, ], F_reparam)
  }
  L_z_z <- temp_res
  ##
  hess <- L_z_z - L_z_gamma %*% solve(L_gamma_gamma, L_gamma_z)
}
lm_part_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X), rev_sign = F) {
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  inner_sum <- 1 + score_eq * (X %*% gamma)
  temp_plog1 <- -mllog(inner_sum, llog_eps, der = 1)[, 2]
  temp_res <- as.vector(rep(0, p))
  for (i in 1:n) {
    temp_res <- temp_res + temp_plog1[i] * crossprod(F_reparam, X[i, ] * as.numeric(X[i, ] %*% gamma))
  }
  ret <- -temp_res
  if (rev_sign) {
    ret <- -ret
  }
  return(ret)
}
lm_part_z_gamma <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X)) {
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  inner_sum <- 1 + score_eq * (X %*% gamma)
  temp_plog1 <- -mllog(inner_sum, llog_eps, der = 1)[, 2]
  temp_plog2 <- -mllog(inner_sum, llog_eps, der = 2)[, 3]
  temp1 <- matrix(0, nrow = p, ncol = ncol(X))
  temp2 <- temp1
  for (i in 1:n) {
    temp1 <- temp1 + temp_plog1[i] * crossprod(F_reparam, tcrossprod(X[i, ]))
    temp2 <- temp2 + temp_plog2[i] * crossprod(F_reparam, X[i, ] * as.numeric(X[i, ] %*% gamma)) %*% t(X[i, ]) * (y[i] - as.numeric(t(X[i, ]) %*% (F_reparam %*% z_reparam + s_hat_reparam)))
  }
  -temp1 - temp2
}
lm_part_z_z <- function(z_reparam, X, F_reparam, s_hat_reparam, y, gamma, llog_eps = 1 / nrow(X)) {
  n <- nrow(X)
  p <- ncol(F_reparam)
  score_eq <- as.vector(y - X %*% (F_reparam %*% z_reparam + s_hat_reparam))
  temp_plog2 <- -mllog(1 + score_eq * (X %*% gamma), llog_eps, der = 2)[, 3]
  temp_res <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
    temp_res <- temp_res + temp_plog2[i] * (crossprod(F_reparam, X[i, ]) * as.numeric(gamma %*% X[i, ])^2) %*% crossprod(X[i, ], F_reparam)
  }
  temp_res
}
lm_fval_gamma <- function(gamma, score_mat, llog_eps = 1 / nrow(score_mat)) {
  sum(mllog(1 + score_mat %*% gamma, llog_eps, der = 0))
}
lm_part_gamma <- function(gamma, score_mat, llog_eps = 1 / nrow(score_mat)) {
  crossprod(score_mat, mllog(1 + score_mat %*% gamma, llog_eps, der = 1)[, 2])
}
##################
