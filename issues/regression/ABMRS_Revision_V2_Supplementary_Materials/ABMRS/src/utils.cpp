// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <algorithm>

#include "utils.h"

// ------ univariate spline regression ------

// MLE of univariate spline regression coefficients
MLERegression mle_regression(const Eigen::VectorXd & x, const Eigen::VectorXd & y,
                             const Eigen::VectorXd & xi, int degree, bool intercept,
                             double xmin, double xmax, int mpart, double eps) {
  // generate B-spline design matrix
  Eigen::MatrixXd B_full = spline(x, xi, degree, intercept, xmin, xmax);
  double m = B_full.rows();
  const Eigen::Index n_cols = B_full.cols();

  // Find columns with non-zero count > mpart
  std::vector<Eigen::Index> valid_cols;
  if (mpart > 0) {
    for (Eigen::Index j = 0; j < n_cols; ++j) {
      const Eigen::Index nnz = (B_full.col(j).array().abs() > eps).count();
      if (nnz > mpart) {
        valid_cols.push_back(j);
      }
    }
  } else {
    // If mpart <= 0, use all columns
    for (Eigen::Index j = 0; j < n_cols; ++j) {
      valid_cols.push_back(j);
    }
  }

  // Extract valid columns into B
  Eigen::MatrixXd B(B_full.rows(), valid_cols.size());
  for (size_t i = 0; i < valid_cols.size(); ++i) {
    B.col(i) = B_full.col(valid_cols[i]);
  }

  // Compute MLE on reduced basis
  Eigen::LLT<Eigen::MatrixXd> llt(B.transpose()*B + 1e-8*Eigen::MatrixXd::Identity(B.cols(), B.cols()));
  Eigen::VectorXd beta_reduced = llt.solve(B.transpose()*y);

  // Map reduced beta back to full size (zero for excluded columns)
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(n_cols);
  for (size_t i = 0; i < valid_cols.size(); ++i) {
    beta[valid_cols[i]] = beta_reduced[i];
  }

  double sigma = (y - B_full*beta).norm() / std::sqrt(m);
  return MLERegression(beta, sigma, llt, beta_reduced, valid_cols);
}

/*
Rcpp::List spline_regression(const Eigen::VectorXd & x, const Eigen::VectorXd & y,
                             const Eigen::VectorXd & xi, int degree, bool intercept) {
  // generate B-spline design matrix
  Rcpp::Environment splines = Rcpp::Environment::namespace_env("splines");
  Rcpp::Function bs = splines["bs"];

  Eigen::MatrixXd B = Rcpp::as<Eigen::MatrixXd>(bs(Rcpp::Named("x")=Rcpp::wrap(x),
                         Rcpp::Named("knots")=Rcpp::wrap(xi),
                         Rcpp::Named("degree")=degree,
                         Rcpp::Named("intercept")=intercept,
                         Rcpp::Named("Boundary.knots")=Rcpp::NumericVector::create(0.0,1.0)));
  Eigen::VectorXd beta = (B.transpose()*B).ldlt().solve(B.transpose()*y);
  double sigma = (y-B*beta).norm() / std::sqrt(x.size());
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("sigma")=sigma);
}
*/
  
/*
Eigen::VectorXd spline_predict(const Eigen::VectorXd & x_new, const Eigen::VectorXd & xi,
                               const Eigen::VectorXd & beta,
                               int degree, bool intercept) {
  // generate B-spline design matrix
  Rcpp::Environment splines = Rcpp::Environment::namespace_env("splines");
  Rcpp::Function bs = splines["bs"];

  Eigen::MatrixXd B = Rcpp::as<Eigen::MatrixXd>(bs(Rcpp::Named("x")=Rcpp::wrap(x_new),
                                                   Rcpp::Named("knots")=Rcpp::wrap(xi),
                                                   Rcpp::Named("degree")=degree,
                                                   Rcpp::Named("intercept")=intercept,
                                                   Rcpp::Named("Boundary.knots")=Rcpp::NumericVector::create(0.0,1.0)));
  Eigen::VectorXd y_new = B*beta;
  return y_new;
}
*/


Eigen::MatrixXd spline(const Eigen::VectorXd & x,
                       const Eigen::VectorXd & xi,
                       int degree, bool intercept,
                       double xmin, double xmax) {
  // generate B-spline design matrix
  Rcpp::Environment splines = Rcpp::Environment::namespace_env("splines");
  Rcpp::Function bs = splines["bs"];

  Eigen::MatrixXd B = Rcpp::as<Eigen::MatrixXd>(
    bs(
      Rcpp::Named("x")=Rcpp::wrap(x),
      Rcpp::Named("knots")=Rcpp::wrap(xi),
      Rcpp::Named("degree")=degree,
      Rcpp::Named("intercept")=intercept,
      Rcpp::Named("Boundary.knots")=Rcpp::NumericVector::create(xmin,xmax)
    )
  );
  return B;
}


/*
// ------ bivariate spline regression ------

// MLE of bivariate spline regression coefficients

MLERegression mle_regression(const Eigen::MatrixXd & x,
                             const Eigen::VectorXd & y,
                             const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                             int degree_1, int degree_2,
                             bool intercept_1, bool intercept_2) {
  // generate B-spline design matrix
  const std::vector<Eigen::VectorXd> xis{xi_1, xi_2};
  const std::vector<int> degrees{degree_1, degree_2};
  const Rcpp::LogicalVector intercepts{intercept_1, intercept_2};
  
  Eigen::MatrixXd B = tensor_spline(x, xis, degrees, intercepts);
  double m = B.rows();
  // compute MLE
  Eigen::LLT<Eigen::MatrixXd> llt(B.transpose()*B + 1e-8*Eigen::MatrixXd::Identity(B.cols(), B.cols()));
  Eigen::VectorXd beta = llt.solve(B.transpose()*y);
  double sigma = (y-B*beta).norm() / std::sqrt(m);
  return MLERegression(beta, sigma, llt);
}


Eigen::MatrixXd bi_tensor_spline(const Eigen::MatrixXd & x,
                              const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                              int degree_1, int degree_2,
                              bool intercept_1, bool intercept_2) {
  std::vector<Eigen::VectorXd> xis = {xi_1, xi_2};
  std::vector<int> degrees = {degree_1, degree_2};
  Rcpp::LogicalVector intercepts = {intercept_1, intercept_2};
  return tensor_spline(x, xis, degrees, intercepts);
}
*/


/*
Eigen::MatrixXd tri_tensor_spline(const Eigen::MatrixXd & x,
                              const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2, const Eigen::VectorXd & xi_3,
                              int degree_1, int degree_2, int degree_3,
                              bool intercept_1, bool intercept_2, bool intercept_3) {
  std::vector<Eigen::VectorXd> xis = {xi_1, xi_2, xi_3};
  std::vector<int> degrees = {degree_1, degree_2, degree_3};
  Rcpp::LogicalVector intercepts = {intercept_1, intercept_2, intercept_3};
  return tensor_spline(x, xis, degrees, intercepts);
}



Rcpp::List cube_spline_regression(const Eigen::MatrixXd & x,
                                     const Eigen::VectorXd & y,
                                     const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2, const Eigen::VectorXd & xi_3,
                                     int degree_1, int degree_2, int degree_3,
                                     bool intercept_1, bool intercept_2, bool intercept_3) {
  Eigen::MatrixXd B = tri_tensor_spline(x, xi_1, xi_2, xi_3, degree_1, degree_2, degree_3, intercept_1, intercept_2, intercept_3);
  Eigen::VectorXd beta = (B.transpose()*B).ldlt().solve(B.transpose()*y);
  double sigma = (y-B*beta).norm() / std::sqrt(x.rows());
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("sigma")=sigma);
}
*/


/*
Rcpp::List surface_spline_regression(const Eigen::MatrixXd & x,
                                     const Eigen::VectorXd & y,
                                     const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                                     int degree_1, int degree_2,
                                     bool intercept_1, bool intercept_2) {
  Eigen::MatrixXd B = bi_tensor_spline(x, xi_1, xi_2, degree_1, degree_2, intercept_1, intercept_2);
  Eigen::VectorXd beta = (B.transpose()*B).ldlt().solve(B.transpose()*y);
  double sigma = (y-B*beta).norm() / std::sqrt(x.rows());
  return Rcpp::List::create(Rcpp::Named("beta")=beta,
                            Rcpp::Named("sigma")=sigma);
}


Eigen::VectorXd surface_spline_predict(const Eigen::MatrixXd & x_new,
                                       const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2,
                                       const Eigen::VectorXd & beta,
                                       int degree_1, int degree_2,
                                       bool intercept_1, bool intercept_2) {
  Eigen::MatrixXd B = bi_tensor_spline(x_new, xi_1, xi_2, degree_1, degree_2, intercept_1, intercept_2);
  Eigen::VectorXd y_new = B*beta;
  return y_new;
}
*/


/*
Eigen::VectorXd cube_spline_predict(const Eigen::MatrixXd & x_new,
                                       const Eigen::VectorXd & xi_1, const Eigen::VectorXd & xi_2, const Eigen::VectorXd & xi_3,
                                       const Eigen::VectorXd & beta,
                                       int degree_1, int degree_2, int degree_3,
                                       bool intercept_1, bool intercept_2, bool intercept_3) {
  Eigen::MatrixXd B = tri_tensor_spline(x_new, xi_1, xi_2, xi_3, degree_1, degree_2, degree_3, intercept_1, intercept_2, intercept_3);
  Eigen::VectorXd y_new = B*beta;
  return y_new;
}
*/



// ----- multivariate spline regression ------

// MLE of multivariate spline regression coefficients
MLERegression mle_regression(const Eigen::MatrixXd & x,
                             const Eigen::VectorXd & y,
                             const std::vector<Eigen::VectorXd> & xis,
                             const std::vector<int> & degrees,
                             const Rcpp::LogicalVector & intercepts,
                             const Eigen::RowVectorXd & xmin,
                             const Eigen::RowVectorXd & xmax,
                             int mpart, double eps) {
  // generate B-spline design matrix
  Eigen::MatrixXd B_full = tensor_spline(x, xis, degrees, intercepts, xmin, xmax);
  double m = B_full.rows();
  const Eigen::Index n_cols = B_full.cols();

  // Find columns with non-zero count > mpart
  std::vector<Eigen::Index> valid_cols;
  if (mpart > 0) {
    for (Eigen::Index j = 0; j < n_cols; ++j) {
      const Eigen::Index nnz = (B_full.col(j).array().abs() > eps).count();
      if (nnz > mpart) {
        valid_cols.push_back(j);
      }
    }
  } else {
    // If mpart <= 0, use all columns
    for (Eigen::Index j = 0; j < n_cols; ++j) {
      valid_cols.push_back(j);
    }
  }

  // Extract valid columns into B
  Eigen::MatrixXd B(B_full.rows(), valid_cols.size());
  for (size_t i = 0; i < valid_cols.size(); ++i) {
    B.col(i) = B_full.col(valid_cols[i]);
  }

  // Compute MLE on reduced basis
  Eigen::LLT<Eigen::MatrixXd> llt(B.transpose()*B + 1e-8*Eigen::MatrixXd::Identity(B.cols(), B.cols()));
  Eigen::VectorXd beta_reduced = llt.solve(B.transpose()*y);

  // Map reduced beta back to full size (zero for excluded columns)
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(n_cols);
  for (size_t i = 0; i < valid_cols.size(); ++i) {
    beta[valid_cols[i]] = beta_reduced[i];
  }

  double sigma = (y - B_full*beta).norm() / std::sqrt(m);
  return MLERegression(beta, sigma, llt, beta_reduced, valid_cols);
}


//' create a general tensor product spline basis matrix for arbitrary dimensions
//' @param x a numeric matrix, (m,d), each row indicates a predictor value.
//' @param xis a list of numeric vectors, each element contains knots for one dimension.
//' @param degrees an integer vector, degrees for each dimension.
//' @param intercepts a logical vector, whether intercepts are included for each dimension.
//' @param xmin a numeric vector, minimum boundary for each input dimension.
//' @param xmax a numeric vector, maximum boundary for each input dimension.
//'
//' @returns a tensor product B-spline basis matrix.
//'
//' @export
//' @examples
//' # 2D example
//' x <- matrix(runif(100), ncol=2)
//' xmin <- c(0,0)
//' xmax <- c(1,1)
//' xis <- list(c(0.3, 0.6), c(0.4, 0.5))
//' B <- tensor_spline(x, xis, c(3,3), c(FALSE,FALSE), xmin, xmax)
//' 
//' # 3D example
//' x <- matrix(runif(150), ncol=3)
//' xmin <- c(0,0,0)
//' xmax <- c(1,1,1)
//' xis <- list(c(0.2, 0.5), c(0.3, 0.7), c(0.4))
//' B <- tensor_spline(x, xis, c(3,3,3), c(TRUE,TRUE,TRUE), xmin, xmax)
// [[Rcpp::export(tensor_spline)]]
Eigen::MatrixXd tensor_spline(const Eigen::MatrixXd & x,
                              const std::vector<Eigen::VectorXd> & xis,
                              std::vector<int> degrees,
                              Rcpp::LogicalVector intercepts,
                              const Eigen::RowVectorXd & xmin,
                              const Eigen::RowVectorXd & xmax) {
  int d = xis.size(); // number of dimensions
  int m = x.rows();   // number of observations

  // assert that input dimensions match
  if(x.cols() != d) {
    throw std::invalid_argument("Input dimension mismatch: x.cols() != number of xis");
  }
  if(degrees.size() != d) {
    throw std::invalid_argument("Input dimension mismatch: degrees.size() != number of xis");
  }
  if(intercepts.size() != d) {
    throw std::invalid_argument("Input dimension mismatch: intercepts.size() != number of xis");
  }
  if(xmin.size() != d || xmax.size() != d) {
    throw std::invalid_argument("Input dimension mismatch: xmin/xmax.size() != number of xis");
  }
  
  // Get R's splines::bs function
  Rcpp::Environment splines = Rcpp::Environment::namespace_env("splines");
  Rcpp::Function bs = splines["bs"];
  
  // Generate B-spline basis for each dimension
  std::vector<Eigen::MatrixXd> bases(d);
  int total_cols = 1;

  for(int j = 0; j < d; j++) {
    bases[j] = Rcpp::as<Eigen::MatrixXd>(
      bs(
        Rcpp::Named("x") = Rcpp::wrap(x.col(j)),
        Rcpp::Named("knots") = Rcpp::wrap(xis[j]),
        Rcpp::Named("degree") = degrees[j],
        Rcpp::Named("intercept") = intercepts[j],
        Rcpp::Named("Boundary.knots") = Rcpp::NumericVector::create(xmin[j], xmax[j])
      )
    );
    total_cols *= bases[j].cols();
  }
  
  // Build tensor product
  Eigen::MatrixXd B(m, total_cols);
  
  for(int i = 0; i < m; i++) {
    // Start with the first dimension
    Eigen::RowVectorXd row_product = bases[0].row(i);
    
    // Kronecker product with remaining dimensions
    for(int j = 1; j < d; j++) {
      row_product = Eigen::kroneckerProduct(row_product, bases[j].row(i)).eval();
    }
    
    B.row(i) = row_product;
  }
  
  return B;
}
