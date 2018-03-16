#include <Rcpp.h>
#include <RcppEigen.h>
#include <queue>
#include <iostream>
#include <math.h>
using namespace Rcpp;
using namespace RcppEigen;

// [[Rcpp::depends(RcppEigen)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::VectorXi;
using namespace std;


int isnan(double x) { return x != x; }
int isinf(double x) { return !isnan(x) && isnan(x - x); }

MatrixXd AtA(const MatrixXd & A) { // transpose of A times itself
  int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint());
}

MatrixXd AAt (const MatrixXd & A) { // A times its transpose
  int m(A.rows());
  return MatrixXd(m, m).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A);
}


double bregman_entropy_cpp(const VectorXd & x_new, 
                           const VectorXd & x_curr) {
  return x_new.dot(x_new.array().log().matrix()) - x_curr.dot(x_curr.array().log().matrix()) - (x_new - x_curr).dot((1.0 + x_curr.array().log()).matrix());
}


// [[Rcpp::export]]
Eigen::VectorXd closed_beta_lm_AL_cpp(const Eigen::Map<Eigen::VectorXd> & wts, 
                               const Eigen::Map<Eigen::VectorXd> & gamma, 
                               const Eigen::Map<Eigen::VectorXd> & y, 
                               const Eigen::Map<Eigen::MatrixXd> & X, 
                               const Eigen::Map<Eigen::MatrixXd> & F_reparam, 
                               const Eigen::Map<Eigen::VectorXd> & s_hat_reparam,
                               const double & dual_step) {
  MatrixXd Xtr_W_X = X.adjoint() * wts.asDiagonal() * X;
  MatrixXd lhs = F_reparam.adjoint() * AtA(Xtr_W_X) * F_reparam;
  VectorXd rhs = F_reparam.adjoint() * Xtr_W_X * (X.adjoint() * (wts.asDiagonal() * y) + std::pow(dual_step, -1) * gamma - Xtr_W_X * s_hat_reparam);
  VectorXd soln = lhs.llt().solve(rhs);
  return soln;
}


//////////////////////////// REL ////////////////////////////

struct md_fun_par {
  VectorXd delta;
  VectorXd gamma;
  MatrixXd Z_mat;
  double dual_step;
};



double fun_wts_REL(const VectorXd & wts,
                   const md_fun_par & extra_par) { 
  VectorXd inner_term = extra_par.Z_mat.adjoint() * (wts.array() * (1.0 - extra_par.delta.array())).matrix();
  return -wts.array().log().sum() + inner_term.dot(extra_par.gamma) + .5 * extra_par.dual_step * (inner_term).squaredNorm();
}



VectorXd grad_wts_REL(const VectorXd & wts, 
                      const md_fun_par & extra_par) {
  ArrayXd delta_term = 1.0 - extra_par.delta.array();
  ArrayXd linear_term = (extra_par.Z_mat * extra_par.gamma).array();
  ArrayXd quad_term = (AAt(extra_par.Z_mat) * (wts.array() * delta_term).matrix()).array();
  return -wts.array().inverse().matrix() + ((linear_term + extra_par.dual_step * quad_term) * delta_term).matrix();
}

// [[Rcpp::export]]
List mirror_descent_REL(const Eigen::Map<Eigen::VectorXd> & wts,
                        const Eigen::Map<Eigen::VectorXd> & delta,
                        const Eigen::Map<Eigen::VectorXd> & gamma, 
                        const Eigen::Map<Eigen::MatrixXd> & Z_mat,
                        const double & dual_step,
                        double mirror_step = 1.0, 
                        const int & maxit = 1000, 
                        const double & mirror_eps = 1e-8,
                        const std::string & tol_type = "fval", 
                        const bool & vary_step = false,
                        const bool & line_search = false, 
                        const double & ls_eps = 1e-10, 
                        const int & ls_maxit = 1000,
                        const double & ls_beta = .75) {
  //// initializing
  md_fun_par extra_par = {delta, gamma, Z_mat, dual_step};
  VectorXd x_curr = wts;
  VectorXd y_curr = wts;
  VectorXd nu_curr = wts;
  VectorXd nu_new = nu_curr;
  VectorXd x_new = x_curr;
  VectorXd y_new = y_curr;
  int n = wts.size();
  int p = gamma.size();
  VectorXd wts_curr = wts;
  VectorXd wts_new = wts_curr;
  bool converged = false;
  int i = 0;
  double theta_curr;
  double theta_new;
  VectorXd grad_update(p);
  // line search
  double fval_y;
  double surr_curr;
  bool ls_converged = false;
  int iter_ls;
  double ls_diff;
  if (line_search) {
    mirror_step = 1.0;
  }
  //
  if (vary_step && !line_search) {
    mirror_step *= std::pow(2.0 * std::log(n), .5);
  }
  List ret;
  ////
  // initial values for tol
  double fval_new;
  double fval_curr;
  double fval_sc;
  double wts_sc;
  double logelr_new;
  double logelr_curr;
  double logelr_sc;
  double tiny = 1e-15;
  double tiny2 = std::pow(tiny, .25);
  if (!tol_type.compare("fval")) {
    fval_curr = fun_wts_REL(wts_curr, extra_par);
  } else if (!tol_type.compare("logelr")) {
    logelr_curr = wts_curr.array().log().sum();
  }
  double tol;
  //
  while (!converged && i < maxit) {
    i += 1;
    theta_curr = 2.0 / (i + 1.0); //automatically use acceleration
    theta_new = 2.0 / (i + 2.0);
    grad_update = grad_wts_REL(y_curr, extra_par);
    if (vary_step && !line_search) {
      mirror_step *= std::pow(i, -.5);
    }
    nu_new = (nu_curr.array() * exp((-mirror_step * grad_update / theta_curr).array())).matrix();
    nu_new = n * nu_new / nu_new.sum();
    x_new = (1.0 - theta_curr) * x_curr + theta_curr * nu_new;
    y_new = (1.0 - theta_new) * x_new + theta_new * nu_new;
    wts_new = x_new;
    if (line_search) {
      fval_new = fun_wts_REL(wts_new, extra_par);
      fval_y = fun_wts_REL(y_curr, extra_par);
      surr_curr = fval_y + grad_update.dot(wts_new - y_curr) + std::pow(mirror_step, -1.0) * bregman_entropy_cpp(wts_new, y_curr);
      ls_converged = fval_new < surr_curr;
      if (!ls_converged || isinf(ls_converged) || isnan(ls_converged)) {
        for (iter_ls = 0; iter_ls < maxit; ++iter_ls) {
          mirror_step *= ls_beta;
          nu_new = (nu_curr.array() * exp((-mirror_step * grad_update / theta_curr).array())).matrix();
          nu_new = n * nu_new / nu_new.sum();
          x_new = (1 - theta_curr) * x_curr + theta_curr * nu_new;
          y_new = (1 - theta_new) * x_new + theta_new * nu_new;
          wts_new = x_new;
          fval_new = fun_wts_REL(wts_new, extra_par);
          surr_curr = fval_y + grad_update.dot(wts_new - y_curr) + std::pow(mirror_step, -1.0) * bregman_entropy_cpp(wts_new, y_curr);
          ls_diff = surr_curr - fval_new;
          if ((!isinf(ls_diff) && !isnan(ls_diff)) && (ls_diff >= 0.0 || (std::abs(ls_diff) < ls_eps))) {
            ls_converged = true;
            break;
          }
        }
      }
    }
    //
    //update current values
    x_curr = wts_new;
    y_curr = y_new;
    nu_curr = nu_new;
    //
    //check tolerance
    if (!tol_type.compare("fval")) {
      if (!line_search) {
        fval_new = fun_wts_REL(wts_new, extra_par);
      }
      fval_sc = std::abs(fval_curr);
      if (fval_sc < tiny) {
        fval_sc = tiny2;
      }
      tol = std::abs(fval_new - fval_curr) / fval_sc;
      fval_curr = fval_new;
    } else if (!tol_type.compare("wts")) {
      wts_sc = wts_curr.array().abs().maxCoeff();
      if (wts_sc < tiny) {
        wts_sc = tiny2;
      }
      tol = (wts_new - wts_curr).array().abs().maxCoeff() / wts_sc;
    } else if(!tol_type.compare("logelr")) {
      logelr_new = wts_new.array().log().sum();
      logelr_sc = std::abs(logelr_new);
      if (logelr_sc < tiny) {
        logelr_sc = tiny2;
      }
      tol = std::abs(logelr_new - logelr_curr) / logelr_sc;
      logelr_curr = logelr_new;
    }
    //
    // update weights, exit if tol
    wts_curr = wts_new;
    if (tol < mirror_eps) {
      converged = true;
    }
    //
  }
  ret["wts"] = wts_curr;
  ret["iter"] = i;
  ret["tol"] = tol;
  ret["converged"] = converged;
  ret["ls_converged"] = ls_converged;
  return ret;
}

//////// functions for quantile thresholding ////////

vector<int> top_i_pq(const VectorXd & v, const int & n) {
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;
  
  for (int i = 0; i != v.size(); ++i) {
    if (pq.size() < n)
      pq.push(Elt(v(i), i));
    else {
      Elt elt = Elt(v(i), i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }
  
  result.reserve(pq.size());
  while (!pq.empty()) {
    result.push_back(pq.top().second);
    pq.pop();
  }
  
  return result;
}

VectorXd quantile_thresh_pq(const VectorXd & inp, const int & thresh_val) {
  int n = inp.size();
  VectorXd zs = VectorXd::Zero(n);
  if (thresh_val > 0) {
    VectorXd inp_abs = inp.array().abs().matrix();
    vector<int> order_idx = top_i_pq(inp_abs, thresh_val);
    int curr_idx;
    for(int i = 0; i < thresh_val; i++) {
      curr_idx = order_idx[i];
      zs(curr_idx) = inp(curr_idx);
    }
  } else {
    zs = inp;
  }
  return zs;
}

////////

VectorXd delta_prox_REL(const VectorXd & inp, const int & thresh_val) {
  int n = inp.size();
  VectorXd zs = VectorXd::Zero(n);
  if (thresh_val > 0) {
    vector<int> top_idx = top_i_pq(inp, thresh_val);
    int curr_idx;
    for(int i = 0; i < top_idx.size(); i++) {
      curr_idx = top_idx[i];
      zs(curr_idx) = std::min(std::max(inp(curr_idx), 0.0), 1.0);
    }
  } else {
    zs = inp;
  }
  return zs;
}

struct REL_fun_par {
  VectorXd wts;
  VectorXd gamma;
  MatrixXd Z_mat;
  double dual_step;
};

double fun_delta_REL_cpp(const VectorXd & delta,
                     const REL_fun_par & extra_par) {
  VectorXd inner_term = extra_par.Z_mat.adjoint() * (extra_par.wts.array() * (1.0 - delta.array())).matrix();
  return inner_term.dot(extra_par.gamma) + .5 * extra_par.dual_step * (inner_term).squaredNorm();
}


VectorXd grad_delta_REL_cpp(const VectorXd & delta,
                            const REL_fun_par & extra_par) {
  return - (extra_par.Z_mat * extra_par.gamma + extra_par.dual_step * AAt(extra_par.Z_mat) * (extra_par.wts.array() * (1.0 - delta.array())).matrix()).cwiseProduct(extra_par.wts);
}

double surrogate_REL_cpp(const VectorXd & delta_new,
                         const VectorXd & delta_curr,
                         const VectorXd & grad_curr,
                         const double & fval_curr,
                         const double & step_size) {
  return fval_curr + grad_curr.dot(delta_new - delta_curr) + 1.0 / (2.0 * step_size) * (delta_new - delta_curr).squaredNorm();
}


struct line_search_par {
  bool line_search;
  double ls_beta;
  double ls_eps;
  int ls_maxit;
};


void REL_line_search(VectorXd & delta_new,
                     const VectorXd & delta_curr,
                     const VectorXd & grad_curr,
                     const REL_fun_par & extra_par,
                     double & step_size,
                     const int & q,
                     const line_search_par ls_inp) {
  double fval_curr = fun_delta_REL_cpp(delta_curr, extra_par);
  double fun_diff = fun_delta_REL_cpp(delta_new, extra_par) - surrogate_REL_cpp(delta_new, delta_curr, grad_curr, fval_curr, step_size);
  int ls_iter = 0;
  while (fun_diff > 0 && abs(fun_diff) > ls_inp.ls_eps && ls_iter < ls_inp.ls_maxit) {
    ls_iter += 1;
    step_size *= ls_inp.ls_beta;
    delta_new = delta_prox_REL(delta_curr - step_size * grad_curr, q);
    fun_diff = fun_delta_REL_cpp(delta_new, extra_par) - surrogate_REL_cpp(delta_new, delta_curr, grad_curr, fval_curr, step_size);
  }
}

// [[Rcpp::export]]
List prox_grad_REL(const Eigen::Map<Eigen::VectorXd> & wts_new,
                   const Eigen::Map<Eigen::MatrixXd> & Z_mat,
                   const Eigen::Map<Eigen::VectorXd> & gamma_curr,
                   const double & dual_step,
                   const Eigen::Map<Eigen::VectorXd> & delta_inp,
                   const int & q,
                   double delta_step = 1.0,
                   const int & maxit = 2000,
                   const double & outlier_eps = 1e-5,
                   const bool & accelerate = true,
                   const bool & line_search = true,
                   const double & ls_eps = 1e-12,
                   const int & ls_maxit = 100,
                   const double & ls_beta = .75) {
  //initializing
  line_search_par ls_inp = {line_search, ls_beta, ls_eps, ls_maxit};
  REL_fun_par extra_par = {wts_new, gamma_curr, Z_mat, dual_step};
  VectorXd delta_new = delta_inp;
  VectorXd delta_curr = delta_inp;
  VectorXd delta_update = delta_inp;
  VectorXd delta_grad_val;
  VectorXd update_term;
  double delta_diff;
  VectorXd nu;
  VectorXd delta_m1 = delta_inp;
  VectorXd delta_m2 = delta_inp;
  bool converged = false;
  int i = 0;
  //
  while (!converged && i < maxit) {
    i += 1;
    if (accelerate) {
      delta_update = delta_m1 + (i - 2.0) / (i + 1.0) * (delta_m1 - delta_m2);
    } else {
      delta_update = delta_curr;
    }
    delta_grad_val = grad_delta_REL_cpp(delta_update, extra_par);
    delta_new = delta_prox_REL(delta_update - delta_step * delta_grad_val, q);
    if (ls_inp.line_search) {
      if (!accelerate) {
        delta_step = 1.0;
      }
      REL_line_search(delta_new, delta_update, delta_grad_val, extra_par, delta_step, q, ls_inp);
    }
    delta_diff = (delta_new - delta_curr).array().abs().maxCoeff();
    delta_curr = delta_new;
    if (accelerate) {
      delta_m2 = delta_m1;
      delta_m1 = delta_curr;
    }
    if (delta_diff < outlier_eps) {
      converged = true;
    }
  }
  List ret;
  ret["delta_opt"] = delta_curr;
  ret["iter"] = i;
  ret["tol"] = delta_diff;
  ret["converged"] = converged;
  return ret;
}


