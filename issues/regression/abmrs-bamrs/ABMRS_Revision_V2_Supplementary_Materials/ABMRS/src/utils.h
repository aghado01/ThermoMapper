#ifndef UTILS_H
#define UTILS_H

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// Structure to hold MLE regression results, 
// including coefficients, sigma, and LLT decomposition
struct MLERegression {
  Eigen::VectorXd beta;
  double sigma;
  Eigen::LLT<Eigen::MatrixXd> llt;
  Eigen::VectorXd beta_reduced;
  std::vector<Eigen::Index> valid_cols;
  
  MLERegression(const Eigen::VectorXd& b, double s, const Eigen::LLT<Eigen::MatrixXd>& l,
                const Eigen::VectorXd& b_red, const std::vector<Eigen::Index>& vc)
    : beta(b), sigma(s), llt(l), beta_reduced(b_red), valid_cols(vc) {}
};



// ------ univariate spline regression ------

// estimate beta and sigma in the univariate spline regression by least square estimation
MLERegression mle_regression(const Eigen::VectorXd & x, const Eigen::VectorXd & y,
                             const Eigen::VectorXd & xi,
                             int degree = 3, bool intercept = false,
                             double xmin = 0.0, double xmax = 1.0,
                             int mpart = 0, double eps = 1e-8);

/*
Rcpp::List spline_regression(const Eigen::VectorXd & x, const Eigen::VectorXd & y,
                             const Eigen::VectorXd & xi,
                             int degree = 3, bool intercept = false);
*/

/*
// predict response values on x_new by spline regression
Eigen::VectorXd spline_predict(const Eigen::VectorXd & x_new, const Eigen::VectorXd & xi,
                               const Eigen::VectorXd & beta,
                               int degree = 3, bool intercept = false);
*/

// generate a univariate B-spline basis matrix
Eigen::MatrixXd spline(const Eigen::VectorXd & x, const Eigen::VectorXd & xi,
                       int degree = 3, bool intercept = false, double xmin = 0.0, double xmax = 1.0);



/*
// ------ bivariate spline regression ------

// estimate beta and sigma in the bivariate spline regression by least square estimation
MLERegression mle_regression(const Eigen::MatrixXd & x,
                             const Eigen::VectorXd & y,
                             const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                             int degree_1, int degree_2,
                             bool intercept_1, bool intercept_2);

                       
// estimate beta and sigma in the surface spline regression by least square estimation
Rcpp::List surface_spline_regression(const Eigen::MatrixXd & x,
                                     const Eigen::VectorXd & y,
                                     const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                                     int degree_1 = 3, int degree_2 = 3,
                                     bool intercept_1 = false, bool intercept_2 = false);
*/


/*
// estimate beta and sigma in the cube spline regression by least square estimation
Rcpp::List cube_spline_regression(const Eigen::MatrixXd & x,
                                  const Eigen::VectorXd & y,
                                  const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2, const Eigen::VectorXd & xi_3,
                                  int degree_1 = 3, int degree_2 = 3, int degree_3 = 3,
                                  bool intercept_1 = false, bool intercept_2 = false, bool intercept_3 = false);
*/

                                  
/*
//' create a bivariate tensor product spline basis matrix
//' @param x a numeric matrix, (m,2), each row indicates a predictor value (x_1,x_2).
//' @param xi_1 a numeric vector, indicates `k_1` knots in x_1.
//' @param xi_2 a numeric vector, indicates `k_2` knots in x_2
//' @param degree_1 int, the degree of polynomial in x_1, default value is `3`.
//' @param degree_2 int, the degree of polynomial in x_2, default value is `3`.
//' @param intercept_1 bool, whether the intercept is included in the basis in x_1,
//' default value is `FALSE`.
//' @param intercept_2 bool, whether the intercept is included in the basis in x_2,
//' default value is `FALSE`.
//'
//' @returns a bivariate tensor product B-spline basis matrix.
//'
//' @export
// [[Rcpp::export(bi_tensor_spline)]]
Eigen::MatrixXd bi_tensor_spline(const Eigen::MatrixXd & x,
                              const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                              int degree_1 = 3, int degree_2 = 3,
                              bool intercept_1 = false, bool intercept_2 = false);
*/


// ------ multivariate spline regression ------

// estimate beta and sigma in the multivariate spline regression by least square estimation
MLERegression mle_regression(const Eigen::MatrixXd & x,
                             const Eigen::VectorXd & y,
                             const std::vector<Eigen::VectorXd> & xis,
                             const std::vector<int> & degrees,
                             const Rcpp::LogicalVector & intercepts,
                             const Eigen::RowVectorXd & xmin,
                             const Eigen::RowVectorXd & xmax,
                             int mpart = 0, double eps = 1e-8);

// General tensor spline for arbitrary dimensions
Eigen::MatrixXd tensor_spline(const Eigen::MatrixXd & x,
                              const std::vector<Eigen::VectorXd> & xis,
                              std::vector<int> degrees,
                              Rcpp::LogicalVector intercepts,
                              const Eigen::RowVectorXd & xmin,
                              const Eigen::RowVectorXd & xmax);


/*
//' create a trivariate tensor product spline basis matrix
//' @param x a numeric matrix, (m,3), each row indicates a predictor value (x_1,x_2,x_3).
//' @param xi_1 a numeric vector, indicates `k_1` knots in x_1.
//' @param xi_2 a numeric vector, indicates `k_2` knots in x_2
//' @param xi_3 a numeric vector, indicates `k_3` knots in x_3
//' @param degree_1 int, the degree of polynomial in x_1, default value is `3`.
//' @param degree_2 int, the degree of polynomial in x_2, default value is `3`.
//' @param degree_3 int, the degree of polynomial in x_3, default value is `3`.
//' @param intercept_1 bool, whether the intercept is included in the basis in x_1,
//' default value is `FALSE`.
//' @param intercept_2 bool, whether the intercept is included in the basis in x_2,
//' default value is `FALSE`.
//' @param intercept_3 bool, whether the intercept is included in the basis in x_3,
//' default value is `FALSE`.
//'
//' @returns a trivariate tensor product B-spline basis matrix.
//'
//' @export
// [[Rcpp::export(tri_tensor_spline)]]
Eigen::MatrixXd tri_tensor_spline(const Eigen::MatrixXd & x,
                              const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2, const Eigen::VectorXd & xi_3,
                              int degree_1 = 3, int degree_2 = 3, int degree_3 = 3,
                              bool intercept_1 = false, bool intercept_2 = false, bool intercept_3 = false);
*/

/*
// predict response values on x_new by surface spline regression
Eigen::VectorXd surface_spline_predict(const Eigen::MatrixXd & x_new,
                                       const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                                       const Eigen::VectorXd & beta,
                                       int degree_1 = 3, int degree_2 = 3,
                                       bool intercept_1 = false, bool intercept_2 = false);
*/


/*
// predict response values on x_new by cube spline regression
Eigen::VectorXd cube_spline_predict(const Eigen::MatrixXd & x_new,
                                    const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2, const Eigen::VectorXd & xi_3,
                                    const Eigen::VectorXd & beta,
                                    int degree_1 = 3, int degree_2 = 3, int degree_3 = 3,
                                    bool intercept_1 = false, bool intercept_2 = false, bool intercept_3 = false);
*/


#endif
