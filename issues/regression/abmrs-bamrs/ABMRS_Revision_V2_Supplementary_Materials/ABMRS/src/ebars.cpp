// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <algorithm>

#include "ebars.h"



double EBARS::_birth() {
  double p = c * std::min(1.0, std::pow((n-k)/(k+1.0), 1.0-gamma));
  return p;
}

double EBARS::_death() {
  double p = c * std::min(1.0, std::pow(k/(n-k+1.0), 1.0-gamma));
  return p;
}

double EBARS::_relocate() {
  double p = 1.0 - _birth() - _death();
  return p;
}

void EBARS::_knots() {
  knots = Eigen::VectorXd::Zero(n);
  double step = 1.0/(n+1);
  knots(0) = step;
  for(int i=1;i<n;i++) {
    knots(i) = knots(i-1) + step;
  }
}

EBARS::EBARS(const Eigen::VectorXd & _x, const Eigen::VectorXd & _y, double _xmin, double _xmax, int _mpart, double _eps,
             Rcpp::NumericVector _para, Rcpp::IntegerVector _num, Rcpp::List _spline) {
  x = _x; y = _y;
  gamma = _para(0); c = _para(1);
  mpart = _mpart; eps = _eps;
  m = _y.size();
  double _times = _para(2);
  // compute the number of knots
  int _k = _num(0);
  fix_k = (_k>0) ? true : false;
  k = (_k>0) ? _k : 1;

  int _n = _num(1);
  n = (_n>0) ? _n : int(m*_times);

  degree = int(_spline["degree"]);
  intercept = _spline["intercept"];

  // transform x to t
  // xmin = x.minCoeff(); xmax = x.maxCoeff();
  xmin = _xmin; xmax = _xmax;
  t = (x.array()-xmin)/(xmax-xmin);
  tmin = 0.0;
  tmax = 1.0;

  _knots();
  _initial();
  // maximum likelihood estimation
  MLERegression mle = mle_regression(t, y, xi, degree, intercept, tmin, tmax, mpart, eps);
  beta_mle = mle.beta;
  sigma_mle = mle.sigma;
  beta_reduced_mle = mle.beta_reduced;
  valid_cols = mle.valid_cols;
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

void EBARS::_initial() {
  Eigen::VectorXi idx = Rcpp::as<Eigen::VectorXi>(Rcpp::sample(n, k)).array()-1;
  Eigen::VectorXi idx_rem = Eigen::VectorXi::LinSpaced(n, 0, n-1);
  xi = Eigen::VectorXd::Zero(k);
  for(int i=0;i<k;i++) {
    xi(i) = knots(idx(i));
    idx_rem(idx(i)) = -1; // remark selected knots as -1
  }

  remain_knots = Eigen::VectorXd::Zero(n-k);
  int j = 0;
  for(int i=0;i<n;i++) {
    if(idx_rem(i)>=0) {
      remain_knots(j++) = knots(idx_rem(i));
    }
  }
}

// Metropolis-Hastings update
void EBARS::_update() {
  double type = Rcpp::runif(1, 0.0, 1.0)[0];
  double birth = _birth();
  double death = _death();
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
    xi_new.head(k) = xi; xi_new(k) = remain_knots(idx);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
    remain_new.head(idx) = remain_knots.head(idx);
    remain_new.tail(n-k_new-idx) = remain_knots.tail(n-k_new-idx);
  }

  // death scheme
  else if(type < (birth+death)) {
    k_new = k - 1;
    int idx = Rcpp::sample(k,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n-k_new);
    xi_new.head(idx) = xi.head(idx);
    xi_new.tail(k_new-idx) = xi.tail(k_new-idx);
    remain_new.head(n-k) = remain_knots; remain_new(n-k) = xi(idx);
  }

  // relocate scheme
  else {
    k_new = k;
    int idx_0 = Rcpp::sample(k,1)[0] - 1;
    int idx_1 = Rcpp::sample(n-k,1)[0] - 1;
    xi_new = xi; remain_new = remain_knots;
    xi_new(idx_0) = remain_knots(idx_1);
    remain_new(idx_1) = xi(idx_0);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
  }

  // compute MLE
  MLERegression mle_new = mle_regression(t, y, xi_new, degree, intercept, tmin, tmax, mpart, eps);
  Eigen::VectorXd beta_mle_new = mle_new.beta;
  double sigma_mle_new = mle_new.sigma;
  Eigen::MatrixXd U_chol_new = mle_new.llt.matrixU();
  Eigen::VectorXd beta_reduced_mle_new = mle_new.beta_reduced;
  std::vector<Eigen::Index> valid_cols_new = mle_new.valid_cols;
  int nv_new = beta_mle_new.size();

  // compute marginal likelihood ratio: m^((k-k')/2)*(RSS_k/RSS_k')^(m/2)
  // double like_ratio = std::pow(m, (k-k_new)/2.0) * std::pow(sigma_mle/sigma_mle_new, m);
  double like_ratio = std::pow(m, (nv-nv_new)/2.0) * std::pow(sigma_mle/sigma_mle_new, m);

  // compute acceptance probability
  double acc_prob = (k > 0) ? like_ratio : birth * like_ratio;

  // decide whether to accept the new state
  double acc_criteria = Rcpp::runif(1, 0.0, 1.0)[0];
  if(acc_criteria < acc_prob) {
    k = k_new; xi = xi_new; remain_knots = remain_new;
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


void EBARS::rjmcmc(int burns, int steps) {
  // burns-in period
  for(int i=0;i<burns;i++) {
    _update();
  }
  // posterior samples of knots
  for(int i=0;i<steps;i++) {
    _update();
    _xis.push_back(xi);
    xis.push_back(xi.array()*(xmax-xmin) + xmin);
    betas.push_back(beta);
    sigmas.push_back(sigma_mle);
  }
}

Eigen::MatrixXd EBARS::predict(const Eigen::VectorXd & x_new, bool nugget) {
  // transform x_new to t_new
  Eigen::VectorXd t_new = (x_new.array()-xmin)/(xmax-xmin);
  // predict for new data using each posterior sample
  Eigen::MatrixXd predictions(t_new.size(), xis.length());
  for(int i=0;i<xis.length();i++) {
    Eigen::VectorXd xi_sample = Rcpp::as<Eigen::VectorXd>(_xis[i]);
    Eigen::VectorXd beta_sample = Rcpp::as<Eigen::VectorXd>(betas[i]);
    Eigen::MatrixXd B = spline(t_new, xi_sample, degree, intercept, tmin, tmax);
    Eigen::VectorXd y_pred = B * beta_sample;
    if(nugget) {
      double sigma_sample = Rcpp::as<double>(sigmas[i]);
      y_pred += Rcpp::as<Eigen::VectorXd>(Rcpp::rnorm(t_new.size(), 0.0, sigma_sample));
    }
    predictions.col(i) = y_pred;
  }
  return predictions;
}

Rcpp::List EBARS::get_knots() {
  return xis;
}

Rcpp::List EBARS::get_coefs() {
  return betas;
}

Rcpp::List EBARS::get_resids() {
  return sigmas;
}

// expose Rcpp class
RCPP_MODULE(class_EBARS) {
  using namespace Rcpp;

  class_<EBARS>("EBARS")
    .constructor<Eigen::VectorXd, Eigen::VectorXd, double, double, int, double, NumericVector, IntegerVector, List>(
      "Construct EBARS object for univariate spline regression"
    )
    .method("rjmcmc", &EBARS::rjmcmc, 
      "Run reversible jump MCMC algorithm")
    .method("predict", &EBARS::predict, 
      "Predict posterior response values for new data")
    .method("knots", &EBARS::get_knots, 
      "Get posterior samples of knots")
    .method("coefs", &EBARS::get_coefs, 
      "Get posterior samples of regression coefficients")
    .method("resids", &EBARS::get_resids, 
      "Get posterior samples of residual standard deviations")
  ;
}

RCPP_EXPOSED_CLASS(EBARS)

