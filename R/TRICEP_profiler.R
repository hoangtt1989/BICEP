##################bootstrap for more general composite tests - using proximal operator
#' Bootstrap TRICEP/BICEP for general composite tests using a proximal operator
#' 
#' @param X a design matrix with observations in rows and variables in columns.
#' @param y a vector of the responses. Default is \code{NULL} in case \code{family = "mean"} (test a vector mean).
#' @param prox_fun a proximal operator. This should be an \code{R} object.
#' @param prox_fun_arg a list of arguments passed to \code{prox_fun}.
#' @param family either the GLM family or \code{"mean"} for a vector mean.
#' @param ... additional arguments passed to \code{\link{TRICEP_glm}}.
#' @param beta_init an initial guess for the \code{beta} vector.
#' @param B the number of boostrap replications.
#' @param prob the coverage probability.
#' @param parallel a boolean indicating whether each replicate should be parallelized. A backend for \code{\link[foreach]{foreach}} must be registered if this is \code{T}.
#' @param seed an integer for setting the seed.
#' @return A list with the bootstrapped test statistic, the p value, and the replicates.
#' @export
TRICEP_glm_boot <- function(X, y = NULL, prox_fun = prox_inequality, prox_fun_arg = list(test_idx = 1, test_val = 0), family = c('gaussian', 'binomial', 'poisson', 'mean'), ..., beta_init = NULL, B = 200, prob = .95, parallel = T, seed = 123) {
  ##checking inputs
  family <- family[1]
  ##
  n <- nrow(X)
  if (is.null(beta_init)) {
    beta_init <- switch(family,
                        gaussian = solve(crossprod(X), crossprod(X, y)),
                        binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                        poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                        mean = colMeans(X))
  }
  set.seed(seed)
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
  } else {
    `%fun%` <- foreach::`%do%`
  }
  i <- NULL
  init_res <- TRICEP_glm(beta_init, X, y, prox_fun = prox_fun, prox_fun_arg = prox_fun_arg, family = family, ...)
  start_time <- Sys.time()
  boot_res <- foreach::foreach(i = 1:B, .combine = cbind) %fun% {
    sample_ind <- sample(1:n, n, replace = T)
    X <- X[sample_ind, ]
    if (family != 'mean') {
      y <- y[sample_ind]
    }
    beta_init <- switch(family,
                        gaussian = solve(crossprod(X), crossprod(X, y)),
                        binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                        poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                        mean = colMeans(X))
    REL_res <- try(TRICEP_glm(beta_init, X, y, prox_fun = prox_fun, prox_fun_arg = prox_fun_arg, family = family, ...))
    stat <- ifelse(class(REL_res) == 'try-error', Inf, REL_res$logelr_stat)
    return(stat)
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  keep_idx <- which(!is.infinite(boot_res))
  keep_len <- length(keep_idx)
  boot_res <- boot_res[keep_idx]
  boot_pval <- length(which(boot_res >= init_res$logelr_stat)) / keep_len
  return(list(boot_stat = stats::quantile(boot_res, prob), boot_pval = boot_pval, boot_vec = as.numeric(boot_res), time = as.double(tot_time, 'secs')))
}
##################

##################bootstrap for coefficient - starting at MLE
#' Bootstrap TRICEP/BICEP for a coefficient, starting at the MLE.
#' 
#' @param fix_index a vector indicating the indices to be boostrapped.
#' @param X a design matrix with observations in rows and variables in columns.
#' @param y a vector of the responses. Default is \code{NULL} in case \code{family = "mean"} (test a vector mean).
#' @param family either the GLM family or \code{"mean"} for a vector mean.
#' @param B the number of boostrap replications.
#' @param prob the coverage probability.
#' @param mean_init the MLE for the coefficients.
#' @param algorithm a string indicating the algorithm to be used. \code{"nested"} is only compatible when \code{family = "gaussian"}.
#' @param parallel a boolean indicating whether each replicate should be parallelized. A backend for \code{\link[foreach]{foreach}} must be registered if this is \code{T}.
#' @param seed an integer for setting the seed.
#' @param ... additional arguments passed to \code{\link{TRICEP_glm}}.
#' @return A list with the bootstrapped test statistic and the replicates.
#' @export
TRICEP_glm_beta_boot <- function(fix_index, X, y = NULL, family = c('gaussian', 'binomial', 'poisson'), B = 200, prob = .95, mean_init = NULL, algorithm = c('TRICEP', 'nested'), parallel = T, seed = 123, ...) {
  ##checking inputs
  family <- family[1]
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('TRICEP', 'nested'))) {
    stop('algorithm must be TRICEP or nested')
  }
  if (algorithm != 'TRICEP' && family != 'gaussian') {
    stop('nested is currently only implemented for gaussian family')
  }
  ##
  ##initial parameters
  p <- ncol(X)
  if (is.null(mean_init)) {
    mean_init <- switch(family,
                        gaussian = solve(crossprod(X), crossprod(X, y)),
                        binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                        poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                        mean = colMeans(X))
  }
  ##
  ##set up the reparam
  reparam_res <- reparam_helper(mean_init, fix_index, 0)
  F_reparam <- reparam_res$F_reparam
  s_hat_reparam <- reparam_res$s_hat_reparam
  alpha_test <- reparam_res$alpha_test
  ##
  set.seed(seed)
  if (parallel) {
    `%fun%` <- doRNG::`%dorng%`
  } else {
    `%fun%` <- foreach::`%do%`
  }
  i <- NULL
  start_time <- Sys.time()
  boot_res <- foreach::foreach(i = 1:B, .combine = cbind) %fun% {
    res <- try(glm_boot_beta_interm_fun(fix_index, alpha_test, X, y, F_reparam, s_hat_reparam, algorithm, family, ...))
    stat <- ifelse(class(res) == 'try-error', Inf, res)
    return(stat)
  }
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  keep_idx <- which(!is.infinite(boot_res))
  keep_len <- length(keep_idx)
  boot_res <- boot_res[keep_idx]
  return(list(boot_stat = stats::quantile(boot_res, prob), boot_vec = as.numeric(boot_res), time = as.double(tot_time, 'secs')))
}

glm_boot_beta_interm_fun <- function(fix_index, alpha_test, X, y, F_reparam, s_hat_reparam, algorithm = c('TRICEP', 'nested'), family, ...) {
  algorithm <- algorithm[1]
  ##resample data
  n <- nrow(X)
  sample_ind <- sample(1:n, n, replace = T)
  X <- X[sample_ind, ]
  y <- y[sample_ind]
  ##
  ##MLE estimate for initial value
  beta_MLE <- switch(family,
                     gaussian = solve(crossprod(X), crossprod(X, y)),
                     binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                     poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                     mean = colMeans(X))
  beta_test <- beta_MLE
  beta_test[fix_index] <- alpha_test
  ##
  ##pass to optimization function
  res <- switch(algorithm,
                TRICEP = TRICEP_glm(beta_test, X, y, F_reparam, s_hat_reparam, family = family, ...),
                nested = CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, ...))
  ##
  ##extract the test statistic
  return(res$logelr_stat)
  ##
}
##################

##################compute EL for a fixed value of a component of the beta vector for glm
#' Compute CEL for a given value of a coefficient.
#' 
#' @param beta_fix_val the test values.
#' @param fix_index a vector indicating the indices to be tested.
#' @param X a design matrix with observations in rows and variables in columns.
#' @param y a vector of the responses. Default is \code{NULL} in case \code{family = "mean"} (test a vector mean).
#' @param family either the GLM family or \code{"mean"} for a vector mean.
#' @param mean_init the MLE for the coefficients.
#' @param algorithm a string indicating the algorithm to be used. \code{"nested"} is only compatible when \code{family = "gaussian"}.
#' @param ... additional arguments passed to \code{\link{TRICEP_glm}}.
#' @return A list with the output from \code{\link{TRICEP_glm}}.
#' @export
TRICEP_glm_beta_fixed <- function(beta_fix_val, fix_index, X, y = NULL, family = c('gaussian', 'binomial', 'poisson', 'mean'), mean_init = NULL, algorithm = c('TRICEP', 'nested'), ...) {
  ##checking inputs
  family <- family[1]
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('TRICEP', 'nested'))) {
    stop('algorithm must be TRICEP or nested')
  }
  if (algorithm != 'TRICEP' && family != 'gaussian') {
    stop('nested is currently only implemented for gaussian family')
  }
  ##
  ##initial parameters
  p <- ncol(X)
  if (is.null(mean_init)) {
    mean_init <- switch(family,
                        gaussian = solve(crossprod(X), crossprod(X, y)),
                        binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                        poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                        mean = colMeans(X))
  }
  alpha_add <- beta_fix_val - mean_init[fix_index] 
  TRICEP_glm_beta_add(alpha_add, fix_index, X, y, family, mean_init, algorithm, ...)
  ##
}
##################

##################compute EL for a fixed value of a component of the beta vector (alpha_add specifies how far away from MLE) for glm
#' Compute CEL by perturbing the MLE at a coefficient.
#' 
#' @param alpha_add a vector with the perturbations.
#' @param fix_index a vector indicating the indices to be tested.
#' @param X a design matrix with observations in rows and variables in columns.
#' @param y a vector of the responses. Default is \code{NULL} in case \code{family = "mean"} (test a vector mean).
#' @param family either the GLM family or \code{"mean"} for a vector mean.
#' @param mean_init the MLE for the coefficients.
#' @param algorithm a string indicating the algorithm to be used. \code{"nested"} is only compatible when \code{family = "gaussian"}.
#' @param ... additional arguments passed to \code{\link{TRICEP_glm}}.
#' @return A list with the output from \code{\link{TRICEP_glm}}.
#' @export
TRICEP_glm_beta_add <- function(alpha_add, fix_index, X, y = NULL, family = c('gaussian', 'binomial', 'poisson', 'mean'), mean_init = NULL, algorithm = c('TRICEP', 'nested'), ...) {
  ##checking inputs
  family <- family[1]
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('TRICEP', 'nested'))) {
    stop('algorithm must be TRICEP or nested')
  }
  if (algorithm != 'TRICEP' && family != 'gaussian') {
    stop('nested is currently only implemented for gaussian family')
  }
  ##
  ##initial parameters
  p <- ncol(X)
  if (is.null(mean_init)) {
    mean_init <- switch(family,
                        gaussian = solve(crossprod(X), crossprod(X, y)),
                        binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                        poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                        mean = colMeans(X))
  }
  ##
  ##set up the reparam
  reparam_res <- reparam_helper(mean_init, fix_index, alpha_add)
  beta_test <- reparam_res$beta_test
  F_reparam <- reparam_res$F_reparam
  s_hat_reparam <- reparam_res$s_hat_reparam
  ##run the chosen algorithm
  res <- switch(algorithm,
                TRICEP = TRICEP_glm(beta_test, X, y, F_reparam, s_hat_reparam, family = family, ...),
                nested = CEL_lm_nested(beta_test, X, y, F_reparam, s_hat_reparam, ...))
  ##
  return(c(list(mean_init = mean_init), res))
}
##################

##################find the EL confidence interval for a component of the beta vector for glm using TRICEP
#' Profile a CEL confidence interval for one component of the coefficient vector.
#' @param fix_index a vector indicating the indices to be tested.
#' @param X a design matrix with observations in rows and variables in columns.
#' @param y a vector of the responses. Default is \code{NULL} in case \code{family = "mean"} (test a vector mean).
#' @param family either the GLM family or \code{"mean"} for a vector mean.
#' @param conf_level the confidence level.
#' @param test_thresh either a string specifying the chi-square approximation or a double (such as a statistic calculated from the bootstrap).
#' @param upper_break the initial guess for the upper limit for root finding. Default uses a value based on the magnitude of the MLE.
#' @param lower_break the initial guess for the lower limit for root finding.
#' @param upper_divide divide the MLE by this much for \code{upper_break}.
#' @param upper_increase increase the initial guess for the upper limit by this amount.
#' @param algorithm a string indicating the algorithm to be used. \code{"nested"} is only compatible when \code{family = "gaussian"}.
#' @param ... additional arguments passed to \code{\link{TRICEP_glm}}.
#' @param verbose a boolean to allow console output.
#' @return A list with confidence interval attributes.
#' @export
TRICEP_glm_beta_profiler <- function(fix_index, X, y = NULL, family = c('gaussian', 'binomial', 'poisson', 'mean'), 
                                  conf_level = .05, test_thresh = 'chisq', upper_break = NULL, lower_break = 1e-6, 
                                  upper_divide = 2, upper_increase = 1e-1, algorithm = c('TRICEP', 'nested'), ..., verbose = F) {
  ##checking inputs
  family <- family[1]
  algorithm <- algorithm[1]
  if (!(algorithm %in% c('TRICEP', 'nested'))) {
    stop('algorithm must be TRICEP or nested')
  }
  if (algorithm != 'TRICEP' && family != 'gaussian') {
    stop('nested is currently only implemented for gaussian family')
  }
  ##
  ##get critical value
  if (test_thresh == 'chisq') {
    test_thresh <- stats::qchisq(1 - conf_level, df = 1)
  }
  ##
  ##get mean value
  beta_MLE <- switch(family,
                     gaussian = solve(crossprod(X), crossprod(X, y)),
                     binomial = stats::glm.fit(X, y, family = binomial(link = 'logit'))$coefficients,
                     poisson = stats::glm.fit(X, y, family = poisson(link = 'log'))$coefficients,
                     mean = colMeans(X))
  beta_fix <- beta_MLE[fix_index]
  ##
  ##make a function for root finding
  root_fun <- function(inp) {
    ret <- TRICEP_glm_beta_add(inp, fix_index, X, y, mean_init = beta_MLE, algorithm = algorithm, family = family, ...)
    ret$logelr_stat - test_thresh
  }
  ##
  #####make sure there is a sign change in the interval
  lower_init_fit <- root_fun(lower_break)
  if (lower_init_fit > 0) {
    stop('lower_break must produce a test statistic below the critical value - decrease it')
  }
  if (!is.null(upper_break)) {
    upper_init_fit <- root_fun(upper_break)
    if (upper_init_fit < 0) {
      stop('upper_break must produce a test statistic above the critical value - increase it')
    }
  } else { ## user did not supply upper_break - we will try to find one using the magnitude of beta
    upper_break <- abs(beta_fix) / upper_divide
    upper_init_fit <- root_fun(upper_break)
    while (upper_init_fit < 1e-2) {
      if (verbose) {
        message('increasing upper_break to produce a test statistic above the critical value')
      }
      upper_break <- upper_break + upper_increase
      upper_init_fit <- root_fun(upper_break)
    }
  } 
  #####
  ###use root finding to get the interval
  ##upper
  start_time <- Sys.time()
  upper_int <- stats::uniroot(root_fun, lower = lower_break, upper = upper_break, check.conv = F)
  end_time <- Sys.time()
  upper_time <- end_time - start_time
  ##
  ##lower
  ###make sure there is a sign change
  upper_break <- -upper_int$root
  while(root_fun(upper_break) < 0) {
    if (verbose) {
      message('decreasing upper_break to produce a test statistic above the critical value (for lower interval)')
    }
    upper_break <- upper_break - upper_increase
  }
  ###
  start_time <- Sys.time()
  lower_int <- stats::uniroot(root_fun, lower = upper_break, upper = -lower_break, check.conv = F)
  end_time <- Sys.time()
  lower_time <- end_time - start_time
  ##
  ###
  return(list(upper_int = beta_fix + upper_int$root, lower_int = beta_fix + lower_int$root, mean_val = beta_fix, time = as.double(upper_time + lower_time, unit = 'secs')))
}
##################
