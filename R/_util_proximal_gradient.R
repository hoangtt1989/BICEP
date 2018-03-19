surrogate_prox <- function(feval, par_new, par_old, grad_val_old, step_size, ...) {
  as.numeric(feval(par_old, ...) + as.numeric(grad_val_old) %*% as.numeric(par_new - par_old) + 1 / (2 * step_size) * sum((par_new - par_old)^2))
}

prox_line_search <- function(feval, par_new, par_old, grad_val_old, step_size, prox_fun, prox_fun_arg = list(), ls_beta = .75, ls_eps = 1e-10, ls_max_iter = 50, ...) {
  fun_diff <- as.numeric(feval(par_new, ...) - surrogate_prox(feval, par_new, par_old, grad_val_old, step_size, ...))
  ls_iter <- 0
  while(fun_diff > 0 && abs(fun_diff) > ls_eps && ls_iter < ls_max_iter) {
    ls_iter <- ls_iter + 1
    step_size <- step_size * ls_beta
    update_term <- par_old - step_size * grad_val_old
    par_new <- do.call('prox_fun', c(list(update_term), prox_fun_arg), quote = T)
    fun_diff <- as.numeric(feval(par_new, ...) - surrogate_prox(feval, par_new, par_old, grad_val_old, step_size, ...))
  }
  return(list(par = par_new, step_size = step_size))
}

proximal_gradient <- function(feval, geval, par, ..., prox_fun, prox_fun_arg = list(), 
                              tol_type = c('par', 'fval'), step_size = 1, accelerate = T, line_search = T, 
                              ls_beta = .75, ls_eps = 1e-10, ls_max_iter = 50, maxit = 1000, prox_eps = 1e-6) {
  ##checking inputs
  tol_type <- tol_type[1]
  ##initializing
  par_curr <- par
  ##
  ## for acceleration
  if (accelerate) {
    par_m1 <- par_curr
    par_m2 <- par_curr
  }
  ##
  fun_curr <- feval(par, ...)
  fun_vals <- rep(0, maxit)
  tiny <- .Machine$double.eps
  tiny2 <- .Machine$double.eps^(.25)
  j <- 0
  converged <- F
  while(!converged && j < maxit) {
    j <- j + 1
    ####acceleration with previous two iterates
    if (accelerate) {
      nu <- par_m1 + (j - 2) / (j + 1) * (par_m1 - par_m2)
      grad_val <- geval(nu, ...)
      update_term <- nu - step_size * grad_val
      par_new <- do.call('prox_fun', c(list(update_term), prox_fun_arg), quote = T)
      if (line_search) {
        ls_res <- prox_line_search(feval, par_new, nu, grad_val, step_size, prox_fun, prox_fun_arg, ls_beta, ls_eps, ls_max_iter, ...)
        par_new <- ls_res$par
        step_size <- ls_res$step_size
      }
    } else {
      step_size <- 1
      grad_val <- geval(par_curr, ...)
      par_new <- do.call('prox_fun', c(list(update_term), prox_fun_arg), quote = T)
      if (line_search) {
        ls_res <- prox_line_search(feval, par_new, par_curr, grad_val, step_size, prox_fun, prox_fun_arg, ls_beta, ls_eps, ls_max_iter, ...)
      }
    }
    fun_new <- feval(par_new, ...)
    fun_sc <- ifelse(abs(fun_curr) < tiny, tiny2, abs(fun_curr))
    fun_tol <- max(abs(fun_new - fun_curr)) / fun_sc
    fun_curr <- fun_new
    fun_vals[j] <- fun_curr
    par_sc <- ifelse(max(abs(par_curr)) < tiny, tiny2, max(abs(par_curr)))
    par_tol <- max(abs(par_new - par_curr)) / par_sc
    ## update values
    par_curr <- par_new
    if (accelerate) {
      par_m2 <- par_m1
      par_m1 <- par_new
    }
    ##
    tol <- switch(tol_type,
                  par = par_tol,
                  fval = fun_tol)
    if (tol < prox_eps && !any(is.na(tol) || is.nan(tol) || is.infinite(tol))) {
      converged <- T
    }
  }
  fun_vals <- fun_vals[1:j]
  return(list(par = par_curr, prox_step = step_size, fval = fun_vals, iter = j, converged = converged))
}
