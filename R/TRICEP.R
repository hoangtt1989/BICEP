
##################affine reparam helper
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
  s_hat_reparam <- pinv(T_reparam) * alpha_test
  F_reparam <- pracma::orth(pracma::nullspace(T_reparam))
  ##
  beta_test <- mean_init
  beta_test[fix_index] <- alpha_test
  return(list(beta_test = beta_test, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, alpha_test = alpha_test, T_reparam = T_reparam))
}
##################

##################wrapper function for mean
#' @export
TRICEP_mean <- function(beta_test, X, F_reparam = NULL, s_hat_reparam = NULL, prox_fun = NULL, 
                          detect_outlier = F, q = NULL, reuse_delta = F,
                          prox_fun_arg = list(), 
                          wts_init = NULL,
                          vary_penalty = c('RB', 'RB2', 'none'), 
                          RB_mu = 10, RB_tau = 2, RB2_ksi = 2, outer_eps = 1e-8, outer_rel_eps = 1e-4, dual_step = 2, 
                          outer_maxit = 1000, dual_type = 1, wts_beta_rep = 1,
                          outer_tol_type = c('primal_dual', 'fval', 'primal'), 
                          outlier_loop_arg = list(),
                          mirror_arg = list(line_search = T), 
                          prox_optim_arg = list(), verbose = F) {
  TRICEP_glm(beta_test, X, y = NULL, F_reparam = F_reparam, s_hat_reparam = s_hat_reparam, prox_fun = prox_fun,
               detect_outlier = detect_outlier, q = q, reuse_delta = reuse_delta, prox_fun_arg = prox_fun_arg,
               family = 'mean', wts_init = wts_init, vary_penalty = vary_penalty, RB_mu = RB_mu, RB_tau = RB_tau, RB2_ksi = RB2_ksi,
               outer_eps = outer_eps, outer_rel_eps = outer_rel_eps, dual_step = dual_step, outer_maxit = outer_maxit, dual_type = dual_type,
               wts_beta_rep = wts_beta_rep, outer_tol_type = outer_tol_type, outlier_loop_arg = outlier_loop_arg, mirror_arg = mirror_arg,
               prox_optim_arg = prox_optim_arg, verbose = verbose)
}
##################


##################function for CEL glm admm (Rcpp)
#' @export
TRICEP_glm <- function(beta_test, X, y, F_reparam = NULL, s_hat_reparam = NULL, prox_fun = NULL, 
                         detect_outlier = F, q = NULL, reuse_delta = F,
                         prox_fun_arg = list(), 
                         family = c('gaussian', 'binomial', 'poisson', 'mean'),
                         wts_init = NULL, beta_opt_type = c('closed_form', 'LBFGS'), 
                         vary_penalty = c('RB', 'RB2', 'none'), 
                         RB_mu = 10, RB_tau = 2, RB2_ksi = 2, outer_eps = 1e-8, outer_rel_eps = 1e-4, dual_step = 2, 
                         outer_maxit = 1000, dual_type = 1, wts_beta_rep = 1,
                         outer_tol_type = c('primal_dual', 'fval', 'primal'), 
                         outlier_loop_arg = list(),
                         mirror_arg = list(line_search = T), lbfgs_arg = list(gtol = .1, invisible = 1), 
                         prox_optim_arg = list(), verbose = F) {
  ##checking inputs
  if (all(is.null(F_reparam), is.null(s_hat_reparam)) && is.null(prox_fun)) {
    stop('The proximal operator (prox_fun) or BOTH F_reparam and s_hat_reparam (for affine reparameterization) must be supplied')
  }
  beta_opt_type <- beta_opt_type[1]
  outer_tol_type <- outer_tol_type[1]
  vary_penalty <- vary_penalty[1]
  family <- family[1]
  primal_dual_flag <- ifelse(vary_penalty %in% c('RB', 'RB2') || outer_tol_type %in% c('primal', 'primal_dual'), T, F)
  if (!(outer_tol_type %in% c('fval', 'primal', 'primal_dual'))) {
    stop('supply a valid outer tolerance type')
  }
  if (!(vary_penalty %in% c('RB', 'RB2', 'none'))) {
    stop('supply a valid penalty varying type')
  }
  canonical_flag <- ifelse(family %in% c('gaussian', 'binomial', 'poisson'), T, F)
  if (canonical_flag) {
    link_fun <- switch(family,
                       gaussian = base::identity,
                       binomial = binomial_link,
                       poisson = exp) 
  }
  if (family != 'gaussian' && beta_opt_type == 'closed_form') {
    beta_opt_type <- 'LBFGS'
    if (verbose) {
      message('Using LBFGS for beta optimization')
    }
  }
  if (!is.null(prox_fun)) {
    beta_opt_type <- 'proximal_gradient'
    if (verbose) {
      message('Using supplied proximal operator (prox_fun) for proximal gradient descent for beta optimization')
    }
  }
  ##
  ##initial parameters
  tiny <- .Machine$double.eps
  tiny2 <- tiny^(.25)
  prob_size <- dim(X)
  n <- prob_size[1]
  p <- prob_size[2]
  beta_curr <- beta_test
  beta_new <- beta_curr
  if (is.null(wts_init)) {
    wts_init <- rep(1, n)
  } else {
    if (length(wts_init) != n) {
      stop('wts_init must be a vector with n elements')
    }
  }
  if (is.null(q)) {
    q <- floor(.05 * n)
  }
  delta_init <- rep(0, n)
  delta_new <- delta_init
  delta_curr <- delta_new
  wts_new <- wts_init
  wts_curr <- wts_new
  gamma_curr <- rep(0, length(beta_test))
  outer_converged <- F
  ##
  ##primal_dual stopping
  if (outer_tol_type == 'primal_dual') {
    primal_eps <- outer_eps * sqrt(p)
    dual_eps <- outer_eps * sqrt(n)
  }
  ##
  ##initialize holders
  outer_fval <- rep(NA, outer_maxit)
  mirror_converged_vec <- rep(NA, outer_maxit)
  beta_converged_vec <- rep(NA, outer_maxit)
  delta_converged_vec <- rep(NA, outer_maxit)
  rho_vec <- rep(NA, outer_maxit)
  ##
  ####precomputing
  if (canonical_flag) {
    RX_curr <- RX_canonical(X, y, beta_curr, link_fun)
  } else if (family == 'mean') {
    RX_curr <- RX_mean(X, beta_curr)
  }
  RX_new <- RX_curr
  outer_fval_curr <- sum(log(wts_new))
  ####
  #################main loop
  start_time <- Sys.time()
  for (j in 1:outer_maxit) {
    if (verbose) {
      if (! j %% 50) {
        message(paste('Iteration', j))
      }
    }
    ##wts and beta updates can be repeated
    for (wts_beta_iter in 1:wts_beta_rep) {
      #######mirror descent
      A_max <- n * dual_step * max(abs(tcrossprod(RX_new))) ##RX %&% t(RX)
      mirror_step <- 1 / ((A_max + max(abs(RX_new %*% gamma_curr))) * max(abs(1 - delta_new)))
      if (!detect_outlier) {
        delta_new <- delta_init
      }
      mirror_res <- do.call('mirror_descent_REL', c(list(wts = wts_new, delta = delta_new, mirror_step = mirror_step, gamma = gamma_curr, Z_mat = RX_new, dual_step = dual_step), mirror_arg), quote = T)
      wts_new <- as.vector(mirror_res$wts)
      #######
      #######delta update
      if (detect_outlier) {
        if (reuse_delta) {
          delta_inp <- delta_new
        } else {
          delta_inp <- delta_init
        }
        outlier_res <- do.call('prox_grad_REL', c(list(wts_new = wts_new, Z_mat = RX_new, gamma_curr = gamma_curr, dual_step = dual_step, delta_inp = delta_inp, q = q), outlier_loop_arg), quote = T)
        delta_new <- as.vector(outlier_res$delta_opt)
        delta_nits <- outlier_res$iter
      }
      #######
      #######beta update
      if (detect_outlier) {
        wts_inp <- wts_new * (1 - delta_new)
      } else {
        wts_inp <- wts_new
      }
      if (canonical_flag) {
        beta_res <- glm_beta_opt(beta_new, wts_inp, gamma_curr, y, X, F_reparam, s_hat_reparam, dual_step, family, beta_opt_type, prox_fun, prox_fun_arg, lbfgs_arg, prox_optim_arg)
        beta_new <- beta_res$beta_new
        RX_new <- beta_res$RX_new
      } else if (family == 'mean') {
        beta_res <- mean_opt_wrapper(beta_new, wts_inp, gamma_curr, X, F_reparam, s_hat_reparam, dual_step, beta_opt_type, prox_fun, prox_fun_arg, prox_optim_arg)
        beta_new <- beta_res$beta_new
        RX_new <- beta_res$RX_new
      }
      #######
    }
    beta_converged_vec[j] <- beta_res$converged
    mirror_converged_vec[j] <- mirror_res$iter
    if (detect_outlier) {
      delta_converged_vec[j] <- outlier_res$iter
    }
    if (verbose) {
      if (!mirror_res$converged) {
        warning(paste('Mirror descent did not converge at iteration', j))
      }
      if (!beta_res$converged) {
        warning(paste('Beta did not converge at iteration', j))
      }
      if (detect_outlier && !outlier_res$converged) {
        warning(paste('Delta did not converge at iteration', j))
      }
    }
    gamma_new <- gamma_curr + dual_step * crossprod(RX_new, wts_inp) ## gamma update
    #######dual residual
    if (primal_dual_flag) {
      if (detect_outlier) {
        dual_resid_res <- REL_dual_resids(RX_new, RX_curr, delta_new, delta_curr, wts_new, gamma_new, dual_step)
      } else {
        dual_resid_res <- CEL_wts_dual_resids(RX_new, RX_curr, wts_new, gamma_new, dual_step)
      }
      dual_resid <- dual_resid_res$dual_resid
      dual_resid_scale <- dual_resid_res$dual_resid_scale
    }
    #######
    #######primal residual
    if (primal_dual_flag) {
      primal_resid <- primal_resid_calc(RX_new, wts_inp)
      if (canonical_flag) {
        primal_resid_scale <- canonical_primal_resid_sc(X, y, wts_new, beta_new, link_fun)
      } else if (family == 'mean') {
        primal_resid_scale <- mean_primal_resid_sc(X, wts_new, beta_new)
      } else {
        primal_resid_scale <- 1
      }
    }
    #######
    ##update R, gamma, beta and wts
    beta_curr <- as.vector(beta_new)
    RX_curr <- RX_new
    wts_curr <- as.vector(wts_new)
    delta_curr <- as.vector(delta_new)
    gamma_curr <- gamma_new
    ##
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
    } else if (vary_penalty == 'RB2') {
      dual_step <- RB2_vary(dual_step, primal_resid, primal_resid_scale, dual_resid, dual_resid_scale, RB2_ksi, RB_tau, RB2_tau)
    } else {
      dual_step <- dual_step
    }
    rho_vec[j] <- dual_step
    #######
    #######stopping
    if (outer_tol_type == 'fval') {
      outer_tol <- outer_fval_diff
    } else if (outer_tol_type == 'primal') {
      outer_tol <- primal_resid
    } else if (outer_tol_type == 'primal_dual') {
      outer_tol <- primal_dual_stopping(primal_resid, dual_resid, primal_eps, dual_eps, primal_resid_scale, dual_resid_scale, outer_rel_eps)
      outer_tol <- ifelse(outer_tol, 0, 1000)
    }
    if (all(outer_tol < outer_eps) && !any(is.nan(outer_tol) || is.na(outer_tol) || is.infinite(outer_tol))) {
      outer_converged <- T
      break
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
  mirror_converged_vec <- mirror_converged_vec[1:j]
  beta_converged_vec <- beta_converged_vec[1:j]
  logelr <- outer_fval_curr
  if (detect_outlier) {
    delta_converged_vec <- delta_converged_vec[1:j]
    delta_opt <- delta_curr
    nz_patt <- which(delta_opt != 0)
  } else {
    delta_converged_vec <- NULL
    delta_opt <- NULL
    nz_patt <- NULL
  }
  if (!primal_dual_flag) {
    primal_resid <- 'not calculated'
    dual_resid <- 'not calculated'
  }
  ####
  
  return(list(logelr = logelr, logelr_stat = -2 * logelr, wts = wts_curr / n, sum_wts = sum(wts_curr) / n,
              gamma = gamma_curr, delta_opt = delta_opt, outlier_idx = nz_patt, beta_opt = beta_curr, outer_converged = outer_converged, time = as.double(tot_time, 'secs'), 
              outer_fval = -outer_fval, primal_resid = primal_resid, dual_resid = dual_resid, rho = rho_vec,
              outer_tol = outer_tol, mirror_converged = mirror_converged_vec, delta_converged = delta_converged_vec, beta_converged = beta_converged_vec, outer_nits = j))
}
##################
