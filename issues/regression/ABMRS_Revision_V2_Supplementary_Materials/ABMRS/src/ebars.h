#ifndef EBARS_H
#define EBARS_H

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "utils.h"



//' @name ClassEBARS
//' @title extended Bayesian adaptive regression univariate spline
//' @description a cpp class
//'
//' @field new constructor, see `ebars`
//' @field rjmcmc run reversible jump MCMC algorithm
//' @field predict predict posterior response values for new data
//' @field knots get posterior samples of knots
//' @field coefs get posterior samples of regression coefficients
//' @field resids get posterior samples of residual standard deviations
class EBARS {
private:
  int k; // the number of underlying knots
  int n; // the number of all potential knots
  int m; // the number of training points
  int nv; // the dimension of spline space
  int degree; // the degree of polynomials
  int mpart; // the minimum number of non-zero points in a basis function
  double eps; // a small constant to define "non-zero"
  double gamma, c; // constants
  double xmin, xmax; // range of x
  double tmin, tmax; // range of t
  bool fix_k; // whether k is given and fixed
  bool intercept; // whether an intercept is included in the basis

  Eigen::VectorXd x;
  Eigen::VectorXd y;
  Eigen::VectorXd t;

  Eigen::VectorXd knots; // all knots
  Eigen::VectorXd xi; // selected knots
  Eigen::VectorXd remain_knots; // remaining knots

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
  double _birth();
  double _death();
  double _relocate();

  void _knots(); // create equally-spaced knot locations
  void _initial(); // initialize

  void _update(); // one step movement in RJMCMC

public:
  // initialize EBARS
  EBARS(const Eigen::VectorXd & _x, const Eigen::VectorXd & _y,
        double _xmin, double _xmax,
        int _mpart, double _eps,
        Rcpp::NumericVector _para,
        Rcpp::IntegerVector _num,
        Rcpp::List _spline);
  void rjmcmc(int burns, int steps); // reversible jump MCMC
  // predict response values on x_new
  Eigen::MatrixXd predict(const Eigen::VectorXd & x_new, bool nugget); // if nugget is true, add noise to the predictions
  Rcpp::List get_knots(); // return posterior knot samples
  Rcpp::List get_coefs(); // return posterior regression coefficients samples
  Rcpp::List get_resids(); // return posterior residual standard deviation samples
};

#endif
