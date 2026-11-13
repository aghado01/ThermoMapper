# wrappers for exposed Rcpp classes


#' @title extended Bayesian adaptive regression univariate spline
#' @description
#' see [ClassEBARS] for more details about class EBARS
#'
#'
#' @param x a numeric vector, indicates the predictor values.
#' @param y a numeric vector, indicates the response values, the same length as x.
#' @param xmin double, the minimum value of predictor, default is `min(x)`.
#' @param xmax double, the maximum value of predictor, default is `max(x)`.
#' @param gamma double, the exponent in extended BIC, default value is `1.0`.
#' It should be constrained to the range `[0.5, 2.0]`.
#' @param c double, the constant affects the proposal distribution in MCMC,
#' default value is `0.3`. It has to be contained in `[0,1]`.
#' @param times double, decides the number of all potential knots if `n<0`.
#' In this case, the number is `length(x) * times`, default value is `3.0`.
#' @param k int, the number of knots is fixed if `k>0`. Otherwise `k` will
#' be estimated in the algorithm, default value is `-1`.
#' @param n int, the number of all potential knots, default value is `-1`.
#' @param degree int, the degree of polynomial spline, default value is `3`.
#' @param intercept bool, whether the intercept is included in the basis,
#' default value is `TRUE`.
#' @param npart int, minimum number of non-zero points in a basis function. 
#' Defaults to `20` or `0.1 * length(x)`, whichever is smaller.
#' @param eps double, a small constant to define "non-zero", default value is `1e-8`.
#' @returns An R6 object of class EBARS with the following methods:
#' \describe{
#'   \item{\code{rjmcmc(burns, steps)}}{Run reversible jump MCMC algorithm.
#'     \itemize{
#'       \item \code{burns}: Number of burn-in iterations, default value is `5000`.
#'       \item \code{steps}: Number of posterior sampling iterations, default value is `5000`.
#'     }}
#'   \item{\code{predict(x_new, nugget)}}{Predict response values for new data.
#'     \itemize{
#'       \item \code{x_new}: Numeric vector of new predictor values
#'       \item \code{nugget}: Boolean indicating whether to add noise to predictions, default is `FALSE`
#'       \item Returns: Matrix of predictions, each column is a posterior sample
#'     }}
#'   \item{\code{knots()}}{Get posterior samples of knot locations.
#'     \itemize{
#'       \item Returns: List of knot vectors from each MCMC iteration
#'     }}
#'   \item{\code{coefs()}}{Get posterior samples of regression coefficients.
#'     \itemize{
#'       \item Returns: List of coefficient vectors from each MCMC iteration
#'     }}
#'   \item{\code{resids()}}{Get posterior samples of residual standard deviations.
#'     \itemize{
#'       \item Returns: List of sigma values from each MCMC iteration
#'     }}
#' }
#' @export
#' @examples
#' library(ABMRS)
#' library(splines)
#' set.seed(1234)
#' knot = c(0.4,0.4,0.4,0.4,0.7)
#' beta = matrix(c(2,-5,5,2,-3,-1,2),ncol=1)
#' m = 200
#'
#' # generate train data set
#' x = c(0,sort(runif(m, 0, 1))[c(-1,-m)],1)
#' B = ns(x,knots=knot,intercept=TRUE,Boundary.knots=c(0,1))
#' y = B%*%beta
#' y_h = y + rnorm(m,0,0.05)
#'
#' # run EBARS
#' a = ebars(x,y_h,times = 2)
#' a$rjmcmc(burns=1000,steps=1000)
#'
#' # generate test data set
#' m_new = 50
#' x_new = runif(m_new,0,1)
#' y_new = predict(B,x_new)%*%beta
#' pred = a$predict(x_new, nugget = FALSE)
#' y_hat = rowMeans(pred)
#'
ebars <- function(x, y, xmin = NULL, xmax = NULL, gamma = 1.0, c = 0.3, times = 3,
                  k = -1, n = -1, degree = 3, intercept = TRUE, npart = NULL, eps = 1e-8) {
  R6EBARS$new(x, y, xmin, xmax, gamma, c, times, k, n, degree, intercept, npart, eps)
}

# R6 wrapper class for EBARS
R6EBARS <- R6::R6Class(
  "EBARS",
  private = list(
    obj = NULL
  ),
  public = list(
    initialize = function(x, y, xmin = NULL, xmax = NULL, gamma = 1.0, c = 0.3, times = 3, 
    k = -1, n = -1, degree = 3, intercept = TRUE, npart = NULL, eps = 1e-8) {
      if (is.null(xmin)) {
        xmin <- min(x)
      }
      if (is.null(xmax)) {
        xmax <- max(x)
      }
      if (is.null(npart)) {
        npart <- min(20, floor(0.1 * length(x)))
      }
      
      private$obj <- EBARS$new(x, y, xmin, xmax,  npart, eps, c(gamma, c, times), c(k, n), 
      list("degree" = degree, "intercept" = intercept))
    },
    rjmcmc = function(burns = 5000, steps = 5000) {
      private$obj$rjmcmc(burns, steps)
    },
    predict = function(x_new, nugget = FALSE) {
      private$obj$predict(x_new, nugget)
    },
    knots = function() {
      private$obj$knots()
    },
    coefs = function() {
      private$obj$coefs()
    },
    resids = function() {
      private$obj$resids()
    }
  )
)




#' @title extended Bayesian adaptive regression multivariate spline
#' @description
#' see [ClassMEBARS] for more details about class MEBARS
#'
#'
#' @param x a numeric matrix, indicates the predictor values, each column is a predictor.
#' @param y a numeric vector, indicates the response values.
#' @param xmin a double vector with the length `ncol(x)`, the minimum values of predictors,
#' default is `apply(x, 2, min)`.
#' @param xmax a double vector with the length `ncol(x)`, the maximum values of predictors,
#' default is `apply(x, 2, max)`.
#' @param gamma double, the exponent in extended BIC, default value is `1.0`.
#' It should be constrained to the range `[0.5, 2.0]`.
#' @param c double, the constant affects the proposal distribution in MCMC,
#' default value is `0.3`. It has to be contained in `[0,1]`.
#' @param times a double vector with the length `ncol(x)`, decides the number of all potential knots if `ns<0`.
#' In this case, the number is `nrow(x) * times`, default value is `c(3, 3, ...)`.
#' @param ks an integer vector with the length `ncol(x)`, the number of knots is fixed if `ks>0`. Otherwise `ks` will
#' be estimated in the algorithm, default value is `c(-1, -1, ...)`.
#' @param ns an integer vector with the length `ncol(x)`, the number of all potential knots, default value is `c(-1, -1, ...)`.
#' @param degrees an integer vector with the length `ncol(x)`, the degree of polynomial spline, default value is `c(3, 3, ...)`.
#' @param intercepts a logical vector with the length `ncol(x)`, whether the intercept is included in the basis,
#' default value is `c(TRUE, TRUE, ...)`.
#' @param npart int, minimum number of non-zero points in a basis function. 
#' Defaults to `20` or `0.1 * nrow(x)`, whichever is smaller.
#' @param eps double, a small constant to define "non-zero", default value is `1e-8`.
#' @returns An R6 object of class MEBARS with the following methods:
#' \describe{
#'   \item{\code{rjmcmc(burns, steps)}}{Run reversible jump MCMC algorithm.
#'     \itemize{
#'       \item \code{burns}: Number of burn-in iterations, standard value is `20000`
#'       \item \code{steps}: Number of posterior sampling iterations, standard value is `20000`
#'     }}
#'   \item{\code{predict(x_new, nugget)}}{Predict response values for new data.
#'     \itemize{
#'       \item \code{x_new}: Numeric matrix of new predictor values
#'       \item \code{nugget}: Boolean indicating whether to add noise to predictions, default is `FALSE`
#'       \item Returns: Matrix of predictions, each column is a posterior sample
#'     }}
#'   \item{\code{knots()}}{Get posterior samples of knot locations.
#'     \itemize{
#'       \item Returns: List of knot vectors from each MCMC iteration
#'     }}
#'   \item{\code{coefs()}}{Get posterior samples of regression coefficients.
#'     \itemize{
#'       \item Returns: List of coefficient vectors from each MCMC iteration
#'     }}
#'   \item{\code{resids()}}{Get posterior samples of residual standard deviations.
#'     \itemize{
#'       \item Returns: List of sigma values from each MCMC iteration
#'     }}
#' }
#' @export
#' @examples
#' library(ABMRS)
#' library(splines)
#' set.seed(1234)
#' # tensor spline
#' beta = matrix(rnorm(40,0,1),ncol=1)
#' fss <- function(x,y){
#'   xix = c(0.2,0.3)
#'   xiy = c(0.5,0.5,0.5,0.5,0.7)
#'   B = tensor_spline(cbind(x,y),list(xix,xiy), c(3,3), c(FALSE,FALSE), c(0,0), c(1,1))
#'   return(B%*%beta)
#' }
#' # parameters' configuration
#' m_train = 1000; m_test = 200
#' noise = 0.1
#' # generate train data set
#' x_1 = c(runif(m_train-2,0,1),0,1)
#' x_2 = c(runif(m_train-2,0,1),0,1)
#' y = fss(x_1,x_2)
#' y_h = y + rnorm(m_train,0,noise)
#' # generate test set
#' x_1_new = runif(m_test,0,1)
#' x_2_new = runif(m_test,0,1)
#' y_new = fss(x_1_new,x_2_new)
#' # run suface EBARS
#' time_start = Sys.time()
#' my_mebars = mebars(cbind(x_1,x_2), y_h)
#' my_mebars$rjmcmc(100,100)
#' time_end = Sys.time()
#' time_end - time_start
#' pred = my_mebars$predict(cbind(x_1_new,x_2_new), nugget = FALSE)
#' y_hat = rowMeans(pred)
#' mean((y_new-y_hat)^2)
#'
mebars <- function(x, y, xmin = NULL, xmax = NULL, gamma = 1.0, c = 0.3, times = NULL,
                   ks = NULL, ns = NULL, degrees = NULL, intercepts = NULL, npart = NULL, eps = 1e-8) {
  R6MEBARS$new(x, y, xmin, xmax, gamma, c, times, ks, ns, degrees, intercepts, npart, eps)
}

# R6 wrapper class for MEBARS
R6MEBARS <- R6::R6Class(
  "MEBARS",
  private = list(
    obj = NULL
  ),
  public = list(
    initialize = function(x, y, xmin = NULL, xmax = NULL, gamma = 1.0, c = 0.3, times = NULL, 
    ks = NULL, ns = NULL, degrees = NULL, intercepts = NULL, npart = NULL, eps = 1e-8) {
      d <- ncol(x)
      if (is.null(xmin)) {
        xmin <- apply(x, 2, min)
      }
      if (is.null(xmax)) {
        xmax <- apply(x, 2, max)
      }
      if (is.null(times)) {
        times <- rep(3, d)
      }
      if (is.null(ks)) {
        ks <- rep(-1, d)
      }
      if (is.null(ns)) {
        ns <- rep(-1, d)
      }
      if (is.null(degrees)) {
        degrees <- rep(3, d)
      }
      if (is.null(intercepts)) {
        intercepts <- rep(TRUE, d)
      }
      if (is.null(npart)) {
        npart <- min(20, floor(0.1 * nrow(x)))
      }
      
      private$obj <- MEBARS$new(x, y, xmin, xmax, npart, eps, c(gamma, c, times), c(ks, ns), 
      list("degrees" = degrees, "intercepts" = intercepts))
    },
    rjmcmc = function(burns = 20000, steps = 20000) {
      private$obj$rjmcmc(burns, steps)
    },
    predict = function(x_new, nugget = FALSE) {
      private$obj$predict(x_new, nugget)
    },
    knots = function() {
      private$obj$knots()
    },
    coefs = function() {
      private$obj$coefs()
    },
    resids = function() {
      private$obj$resids()
    }
  )
)