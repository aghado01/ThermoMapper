#ifndef BINEBARS_H
#define BINEBARS_H

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include "utils.h"

/*
//' @name ClassBinEBARS
//' @title bivariate extended Bayesian adaptive regression spline
//' @description a class of [binebars()]
//'
//' @field new constructor, see `binebars`
//' @field rjmcmc reversible jump MCMC, a wrapper is `mcmc`
//' @field predict predict by surface spline regression with EBARS
//' @field knots return estimated knots
//' @field samples return posterior samples
class BinEBARS {
private:
  int k_1, k_2; // the number of underlying knots in the direction x_1 and x_2
  int n_1, n_2; // the number of all potential knots in the direction x_1 and x_2
  int m; // the number of training points
  int degree_1, degree_2; // the degree of polynomials
  double gamma, c; // constants
  Eigen::RowVector2d xmin, xmax; // range of x
  bool fix_k_1, fix_k_2; // whether k is given and fixed
  bool intercept_1, intercept_2; // whether an intercept is included in the basis

  Eigen::MatrixXd x;
  Eigen::VectorXd y;
  Eigen::MatrixXd t;

  Eigen::VectorXd knots_1, knots_2; // all knots
  Eigen::VectorXd xi_1, xi_2; // selected knots
  Eigen::VectorXd remain_knots_1, remain_knots_2; // remaining knots
  Rcpp::List xis_1, xis_2; // posterior samples

  Eigen::VectorXd beta; // estimated beta by MLE
  double sigma; // estimated sigma by MLE


  // compute probabilities of the next movement type
  double _birth_1(); double _birth_2();
  double _death_1(); double _death_2();
  double _relocate_1(); double _relocate_2();

  void _knots(); // create equally-spaced knot locations
  void _initial(); // initialize

  bool _jump_1(); // one possible jump in the direction x_1, may be rejected
  bool _jump_2(); // one possible jump in the direction x_2, may be rejected
  void _update(); // one step movement in RJMCMC

public:
  // initialize BinEBARS
  BinEBARS(const Eigen::MatrixXd & _x, const Eigen::VectorXd & _y,
           Rcpp::NumericVector _para = Rcpp::NumericVector::create(1.0,0.3,2.0,2.0),
           Rcpp::IntegerVector _num = Rcpp::IntegerVector::create(-1,-1,-1,-1),
           Rcpp::List _spline = Rcpp::List::create(Rcpp::Named("degree_1") = 3,
                                                   Rcpp::Named("degree_2") = 3,
                                                   Rcpp::Named("intercept_1") = false,
                                                   Rcpp::Named("intercept_2") = false));
  void rjmcmc(int burns = 200, int steps = 200); // reversible jump MCMC
  // predict response values on x_new
  Eigen::VectorXd predict(const Eigen::MatrixXd & x_new);
  Rcpp::List get_knots(); // return estimated knots
  Rcpp::List get_samples(); // return posterior samples
};
*/


#endif
