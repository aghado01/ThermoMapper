#ifndef MEBARS_H
#define MEBARS_H

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "utils.h"


//' @name ClassMEBARS
//' @title extended Bayesian adaptive regression multivariate spline
//' @description a cpp class
//'
//' @field new constructor, see `mebars`
//' @field rjmcmc run reversible jump MCMC algorithm
//' @field predict predict posterior response values for new data
//' @field knots get posterior samples of knots
//' @field coefs get posterior samples of regression coefficients
//' @field resids get posterior samples of residual standard deviations
class MEBARS {
private:
  std::vector<int> ks; // the number of underlying knots
  std::vector<int> ns; // the number of all potential knots
  int m; // the number of training points
  int d; // the dimension of predictors
  int nv; // the dimension of tensor product spline spaces
  std::vector<int> degrees; // the degree of polynomials
  int mpart; // the minimum number of non-zero points in a basis function
  double eps; // a small constant to define "non-zero"
  double gamma, c; // constants
  Eigen::RowVectorXd xmin, xmax; // range of x
  Eigen::RowVectorXd tmin, tmax; // range of t
  Rcpp::LogicalVector fix_ks; // whether k is given and fixed
  Rcpp::LogicalVector intercepts; // whether an intercept is included in the basis

  Eigen::MatrixXd x;
  Eigen::VectorXd y;
  Eigen::MatrixXd t;

  std::vector<Eigen::VectorXd> all_knots; // all knots
  std::vector<Eigen::VectorXd> xi; // selected knots
  std::vector<Eigen::VectorXd> remain_knots; // remaining knots

  Rcpp::List _xis; // posterior samples of normalized selected knots
  Rcpp::List xis; // posterior samples of selected knots
  Rcpp::List betas; // posterior samples of regression coefficients
  Rcpp::List sigmas; // posterior samples of residual standard deviation

  Eigen::VectorXd beta_mle; // estimated beta by MLE
  Eigen::VectorXd beta_reduced_mle; // estimated reduced beta by MLE
  std::vector<Eigen::Index> valid_cols; // valid columns in the design matrix for reduced model
  double sigma_mle; // estimated residual standard deviation by MLE
  Eigen::MatrixXd U_chol; // Upper triangular Cholesky factor of the design matrix' cross-product BtB
  Eigen::VectorXd beta; // posterior regression coefficients by Bayesian inference


  // compute probabilities of the next movement type
  double _birth(int i);
  double _death(int i);
  double _relocate(int i);

  void _knots(); // create equally-spaced knot locations
  void _initial(); // initialize

  void _update(); // one step movement in RJMCMC

public:
  // initialize MEBARS
  MEBARS(const Eigen::MatrixXd & _x, const Eigen::VectorXd & _y,
         const Eigen::RowVectorXd & _xmin, const Eigen::RowVectorXd & _xmax,
         int _mpart, double _eps,
         Rcpp::NumericVector _para,
         Rcpp::IntegerVector _num,
         Rcpp::List _spline);
  void rjmcmc(int burns, int steps); // reversible jump MCMC
  // predict response values on x_new
  Eigen::MatrixXd predict(const Eigen::MatrixXd & x_new, bool nugget); // if nugget is true, add noise to the predictions
  Rcpp::List get_knots(); // return posterior knot samples
  Rcpp::List get_coefs(); // return posterior regression coefficients samples
  Rcpp::List get_resids(); // return posterior residual standard deviation samples
};

#endif
