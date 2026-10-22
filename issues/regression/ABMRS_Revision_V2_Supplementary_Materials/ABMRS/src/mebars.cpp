// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <algorithm>

#include "mebars.h"

double MEBARS::_birth(int i) {
  double p = c * std::min(1.0, std::pow((ns[i]-ks[i])/(ks[i]+1.0), 1.0-gamma));
  return p;
}

double MEBARS::_death(int i) {
  double p = c * std::min(1.0, std::pow(ks[i]/(ns[i]-ks[i]+1.0), 1.0-gamma));
  return p;
}

double MEBARS::_relocate(int i) {
  double p = 1.0 - _birth(i) - _death(i);
  return p;
}

void MEBARS::_knots() {
  // initial knots in each dimension
  for(int dim=0;dim<d;dim++) {
    Eigen::VectorXd knots_dim = Eigen::VectorXd::Zero(ns[dim]);
    double step = 1.0/(ns[dim]+1);
    knots_dim(0) = step;
    for(int i=1;i<ns[dim];i++) {
      knots_dim(i) = knots_dim(i-1) + step;
    }
    all_knots.push_back(knots_dim);
  }
}

MEBARS::MEBARS(const Eigen::MatrixXd & _x, const Eigen::VectorXd & _y,
               const Eigen::RowVectorXd & _xmin, const Eigen::RowVectorXd & _xmax, int _mpart, double _eps,
               Rcpp::NumericVector _para, Rcpp::IntegerVector _num, Rcpp::List _spline) {
  x = _x; y = _y;
  gamma = _para(0); c = _para(1);
  m = _y.size(); d = x.cols();
  mpart = _mpart;
  eps = _eps;
  // compute the number of knots
  ks = std::vector<int>(d);
  fix_ks = Rcpp::LogicalVector(d);
  for(int i=0;i<d;i++) {
    int _k = _num(i);
    fix_ks[i] = (_k>0) ? true : false;
    ks[i] = (_k>0) ? _k : 1;
  }

  ns = std::vector<int>(d);
  for(int i=0;i<d;i++) {
    int _n = _num(i+d);
    ns[i] = (_n>0) ? _n : int(m*_para(i+2));
  }

  degrees = Rcpp::as<std::vector<int>>(_spline["degrees"]);
  intercepts = Rcpp::as<Rcpp::LogicalVector>(_spline["intercepts"]);

  // transform x to t
  // xmin = x.colwise().minCoeff(); xmax = x.colwise().maxCoeff();
  xmin = _xmin; xmax = _xmax;
  t = (x.rowwise() - xmin).array().rowwise() / (xmax-xmin).array();
  tmin = Eigen::RowVectorXd::Zero(d);
  tmax = Eigen::RowVectorXd::Ones(d);

  _knots();
  _initial();
  // maximum likelihood estimation
  MLERegression mle = mle_regression(t, y, xi, degrees, intercepts, tmin, tmax, mpart, eps);
  beta_mle = mle.beta;
  beta_reduced_mle = mle.beta_reduced;
  valid_cols = mle.valid_cols;
  sigma_mle = mle.sigma;
  U_chol = mle.llt.matrixU(); // Store the upper triangular factor
  nv = beta_mle.size();
  // sample beta from posterior
  double shrink_factor = m / (m + 1.0);
  Eigen::VectorXd beta_mean = shrink_factor * beta_reduced_mle;
  Eigen::VectorXd z = Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(beta_reduced_mle.size()));
  Eigen::VectorXd beta_reduced = beta_mean + std::sqrt(shrink_factor) * sigma_mle * U_chol.triangularView<Eigen::Upper>().solve(z);
  // expand beta_reduced to full beta
  beta = Eigen::VectorXd::Zero(nv);
  for(int i=0;i<valid_cols.size();i++) {
    beta(valid_cols[i]) = beta_reduced(i);
  }
  
  xis = Rcpp::List::create();
  betas = Rcpp::List::create();
  sigmas = Rcpp::List::create();
}

void MEBARS::_initial() {
  // select initial knots in each dimension
  for(int dim=0;dim<d;dim++) {
    int k = ks[dim];
    int n = ns[dim];
    Eigen::VectorXi idx = Rcpp::as<Eigen::VectorXi>(Rcpp::sample(n, k)).array()-1;
    Eigen::VectorXi idx_rem = Eigen::VectorXi::LinSpaced(n, 0, n-1);
    Eigen::VectorXd xi_dim = Eigen::VectorXd::Zero(k);
    for(int i=0;i<k;i++) {
      xi_dim(i) = all_knots[dim](idx(i));
      idx_rem(idx(i)) = -1; // remark selected knots as -1
    }
    xi.push_back(xi_dim);

    Eigen::VectorXd remain_knots_dim = Eigen::VectorXd::Zero(n-k);
    int j = 0;
    for(int i=0;i<n;i++) {
      if(idx_rem(i)>=0) {
        remain_knots_dim(j++) = all_knots[dim](idx_rem(i));
      }
    }
    remain_knots.push_back(remain_knots_dim);
  }
}

// Metropolis-Hastings update
void MEBARS::_update() {
  // randomly select a dimension to update
  int dim = Rcpp::sample(d,1)[0] - 1;
  int k = ks[dim];
  int n = ns[dim];
  bool fix_k = fix_ks[dim];

  double type = Rcpp::runif(1, 0.0, 1.0)[0];
  double birth = _birth(dim);
  double death = _death(dim);
  // when k = 0, only birth is allowed
  if(k==0) {type = 0.0;}
  // when k is fixed, only relocate is allowed
  if(fix_k) {birth = 0.0; death = 0.0;}

  Eigen::VectorXd xi_new, remain_new;
  int k_new;

  // birth scheme
  if(type < birth) {
    k_new = k + 1;
    int idx = Rcpp::sample(n-k,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n-k_new);
    xi_new.head(k) = xi[dim]; xi_new(k) = remain_knots[dim](idx);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
    remain_new.head(idx) = remain_knots[dim].head(idx);
    remain_new.tail(n-k_new-idx) = remain_knots[dim].tail(n-k_new-idx);
  }

  // death scheme
  else if(type < (birth+death)) {
    k_new = k - 1;
    int idx = Rcpp::sample(k,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n-k_new);
    xi_new.head(idx) = xi[dim].head(idx);
    xi_new.tail(k_new-idx) = xi[dim].tail(k_new-idx);
    remain_new.head(n-k) = remain_knots[dim]; remain_new(n-k) = xi[dim](idx);
  }

  // relocate scheme
  else {
    k_new = k;
    int idx_0 = Rcpp::sample(k,1)[0] - 1;
    int idx_1 = Rcpp::sample(n-k,1)[0] - 1;
    xi_new = xi[dim]; remain_new = remain_knots[dim];
    xi_new(idx_0) = remain_knots[dim](idx_1);
    remain_new(idx_1) = xi[dim](idx_0);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
  }

  // compute MLE
  std::vector<Eigen::VectorXd> xi_temp = xi;
  xi_temp[dim] = xi_new;
  MLERegression mle_new = mle_regression(t, y, xi_temp, degrees, intercepts, tmin, tmax, mpart, eps);
  Eigen::VectorXd beta_mle_new = mle_new.beta;
  Eigen::VectorXd beta_reduced_mle_new = mle_new.beta_reduced;
  std::vector<Eigen::Index> valid_cols_new = mle_new.valid_cols;
  double sigma_mle_new = mle_new.sigma;
  Eigen::MatrixXd U_chol_new = mle_new.llt.matrixU();
  int nv_new = beta_mle_new.size();

  // compute marginal likelihood ratio
  double like_ratio = std::pow(m, (nv-nv_new)/2.0) * std::pow(sigma_mle/sigma_mle_new, m);

  // compute acceptance probability
  double acc_prob = (k > 0) ? like_ratio : birth * like_ratio;

  // decide whether to accept the new state
  double acc_criteria = Rcpp::runif(1, 0.0, 1.0)[0];
  if(acc_criteria < acc_prob) {
    ks[dim] = k_new; xi[dim] = xi_new; remain_knots[dim] = remain_new;
    beta_mle = beta_mle_new; sigma_mle = sigma_mle_new; U_chol = U_chol_new; nv = nv_new;
    beta_reduced_mle = beta_reduced_mle_new; valid_cols = valid_cols_new;
  }

  // sample beta from posterior
  double shrink_factor = m / (m + 1.0);
  Eigen::VectorXd beta_mean = shrink_factor * beta_reduced_mle;
  Eigen::VectorXd z = Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(beta_reduced_mle.size()));
  Eigen::VectorXd beta_reduced = beta_mean + std::sqrt(shrink_factor) * sigma_mle * U_chol.triangularView<Eigen::Upper>().solve(z);
  // expand beta_reduced to full beta
  beta = Eigen::VectorXd::Zero(nv);
  for(int i=0;i<valid_cols.size();i++) {
    beta(valid_cols[i]) = beta_reduced(i);
  }
}


void MEBARS::rjmcmc(int burns, int steps) {
  // burns-in period
  for(int i=0;i<burns;i++) {
    _update();
  }
  // posterior samples of knots
  for(int i=0;i<steps;i++) {
    _update();
    
    _xis.push_back(xi);
    std::vector<Eigen::VectorXd> xi_scaled(xi.size());
    for (int j = 0; j < xi.size(); j++) {
      xi_scaled[j] = xi[j].array() * (xmax[j] - xmin[j]) + xmin[j];
    }
    xis.push_back(xi_scaled);
    betas.push_back(beta);
    sigmas.push_back(sigma_mle);
  }
}

Eigen::MatrixXd MEBARS::predict(const Eigen::MatrixXd & x_new, bool nugget) {
  // transform x_new to t_new
  Eigen::MatrixXd t_new = (x_new.rowwise() - xmin).array().rowwise() / (xmax-xmin).array();
  // predict for new data using each posterior sample
  Eigen::MatrixXd predictions(t_new.rows(), xis.length());
  for(int i=0;i<xis.length();i++) {
    std::vector<Eigen::VectorXd> xi_sample = Rcpp::as<std::vector<Eigen::VectorXd>>(_xis[i]);
    Eigen::VectorXd beta_sample = Rcpp::as<Eigen::VectorXd>(betas[i]);
    Eigen::MatrixXd B = tensor_spline(t_new, xi_sample, degrees, intercepts, tmin, tmax);
    Eigen::VectorXd y_pred = B * beta_sample;
    if(nugget) {
      double sigma_sample = Rcpp::as<double>(sigmas[i]);
      y_pred += Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(t_new.rows(), 0.0, sigma_sample));
    }
    predictions.col(i) = y_pred;
  }
  return predictions;
}

Rcpp::List MEBARS::get_knots() {
  return xis;
}

Rcpp::List MEBARS::get_coefs() {
  return betas;
}

Rcpp::List MEBARS::get_resids() {
  return sigmas;
}

// expose Rcpp class
RCPP_MODULE(class_MEBARS) {
  using namespace Rcpp;

  class_<MEBARS>("MEBARS")
    .constructor<Eigen::MatrixXd, Eigen::VectorXd, Eigen::RowVectorXd, Eigen::RowVectorXd, int, double, NumericVector, IntegerVector, List>(
      "Construct MEBARS object for multivariate spline regression"
    )
    .method("rjmcmc", &MEBARS::rjmcmc, 
      "Run reversible jump MCMC algorithm")
    .method("predict", &MEBARS::predict, 
      "Predict posterior response values for new data")
    .method("knots", &MEBARS::get_knots, 
      "Get posterior samples of knots")
    .method("coefs", &MEBARS::get_coefs, 
      "Get posterior samples of regression coefficients")
    .method("resids", &MEBARS::get_resids, 
      "Get posterior samples of residual standard deviations")
  ;
}

RCPP_EXPOSED_CLASS(MEBARS)


