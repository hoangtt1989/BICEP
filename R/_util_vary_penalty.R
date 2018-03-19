##functions for residual balancing
RB_vary <- function(dual_step, primal_resid, dual_resid, RB_mu, RB_tau) {
  if (primal_resid > RB_mu * dual_resid) {
    dual_step <- RB_tau * dual_step
  } else if (dual_resid > RB_mu * primal_resid) {
    dual_step <- dual_step / RB_tau
  }
  return(dual_step)
}

RB2_vary <- function(dual_step, primal_resid, primal_resid_scale, dual_resid, dual_resid_scale, RB2_ksi, RB_tau) {
  ###varying tau
  primal_dual_div <- primal_resid / dual_resid
  dual_primal_div <- 1 / primal_dual_div
  sqrt_primal_dual_div <- sqrt(primal_dual_div / RB2_ksi)
  sqrt_dual_primal_div <- sqrt(dual_primal_div * RB2_ksi)
  if (sqrt_primal_dual_div < RB_tau && sqrt_primal_dual_div >= 1) {
    RB2_tau <- sqrt_primal_dual_div
  } else if (sqrt_dual_primal_div < 1 && sqrt_dual_primal_div > (1 / RB_tau)) {
    RB2_tau <- sqrt_primal_dual_div
  } else {
    RB2_tau <- RB_tau
  }
  ###
  primal_resid_rel <- primal_resid / primal_resid_scale
  dual_resid_rel <- dual_resid / dual_resid_scale
  if (primal_resid_rel > RB2_ksi * RB_mu * dual_resid_rel) {
    dual_step <- RB2_tau * dual_step
  } else if (dual_resid_rel > RB_mu * primal_resid_rel / RB2_ksi) {
    dual_step <- dual_step / RB2_tau
  }
  return(dual_step)
}

primal_dual_stopping <- function(primal_resid, dual_resid, primal_eps, dual_eps, primal_resid_scale, dual_resid_scale, outer_rel_eps) {
  primal_eps_scale <- primal_eps + outer_rel_eps * primal_resid_scale
  dual_eps_scale <- dual_eps + outer_rel_eps * dual_resid_scale
  converged <- ifelse(primal_resid < primal_eps_scale && dual_resid < dual_eps_scale, T, F)
}
