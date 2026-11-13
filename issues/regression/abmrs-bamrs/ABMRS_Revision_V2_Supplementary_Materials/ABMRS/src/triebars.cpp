// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <algorithm>

#include "triebars.h"


/*
double TriEBARS::_birth_1() {
  double p = c * std::min(1.0, std::pow((n_1-k_1)/(k_1+1.0), 1.0-gamma));
  return p;
}

double TriEBARS::_birth_2() {
  double p = c * std::min(1.0, std::pow((n_2-k_2)/(k_2+1.0), 1.0-gamma));
  return p;
}

double TriEBARS::_birth_3() {
  double p = c * std::min(1.0, std::pow((n_3-k_3)/(k_3+1.0), 1.0-gamma));
  return p;
}

double TriEBARS::_death_1() {
  double p = c * std::min(1.0, std::pow(k_1/(n_1-k_1+1.0), 1.0-gamma));
  return p;
}

double TriEBARS::_death_2() {
  double p = c * std::min(1.0, std::pow(k_2/(n_2-k_2+1.0), 1.0-gamma));
  return p;
}

double TriEBARS::_death_3() {
  double p = c * std::min(1.0, std::pow(k_3/(n_3-k_3+1.0), 1.0-gamma));
  return p;
}

double TriEBARS::_relocate_1() {
  double p = 1.0 - _birth_1() - _death_1();
  return p;
}

double TriEBARS::_relocate_2() {
  double p = 1.0 - _birth_2() - _death_2();
  return p;
}

double TriEBARS::_relocate_3() {
  double p = 1.0 - _birth_3() - _death_3();
  return p;
}

void TriEBARS::_knots() {
  // initial knots in x_1
  knots_1 = Eigen::VectorXd::Zero(n_1);
  double step_1 = 1.0/(n_1+1);
  knots_1(0) = step_1;
  for(int i=1;i<n_1;i++) {
    knots_1(i) = knots_1(i-1) + step_1;
  }

  // initial knots in x_2
  knots_2 = Eigen::VectorXd::Zero(n_2);
  double step_2 = 1.0/(n_2+1);
  knots_2(0) = step_2;
  for(int i=1;i<n_2;i++) {
    knots_2(i) = knots_2(i-1) + step_2;
  }

  // initial knots in x_3
  knots_3 = Eigen::VectorXd::Zero(n_3);
  double step_3 = 1.0/(n_3+1);
  knots_3(0) = step_3;
  for(int i=1;i<n_3;i++) {
    knots_3(i) = knots_3(i-1) + step_3;
  }
}

TriEBARS::TriEBARS(const Eigen::MatrixXd & _x, const Eigen::VectorXd & _y,
                   Rcpp::NumericVector _para, Rcpp::IntegerVector _num, Rcpp::List _spline) {
  x = _x; y = _y;
  gamma = _para(0); c = _para(1);
  m = _y.size();
  // compute the number of knots
  if(_num(0)>0) {
    k_1 = _num(0); fix_k_1 = true;
  } else {
    k_1 = 1; fix_k_1 = false;
  }

  if(_num(1)>0) {
    k_2 = _num(1); fix_k_2 = true;
  } else {
    k_2 = 1; fix_k_2 = false;
  }

  if(_num(2)>0) {
    k_3 = _num(2); fix_k_3 = true;
  } else {
    k_3 = 1; fix_k_3 = false;
  }

  n_1 = (_num(3)>0) ? _num(3) : int(m*_para(2));
  n_2 = (_num(4)>0) ? _num(4) : int(m*_para(3));
  n_3 = (_num(5)>0) ? _num(5) : int(m*_para(4));

  degree_1 = int(_spline["degree_1"]); degree_2 = int(_spline["degree_2"]); degree_3 = int(_spline["degree_3"]);
  intercept_1 = _spline["intercept_1"]; intercept_2 = _spline["intercept_2"]; intercept_3 = _spline["intercept_3"];

  // transform x to t
  xmin = x.colwise().minCoeff(); xmax = x.colwise().maxCoeff();
  t = (x.rowwise() - xmin).array().rowwise() / (xmax-xmin).array();

  _knots();
  _initial();
  // maximum likelihood estimation
  Rcpp::List pars = cube_spline_regression(t,y,xi_1,xi_2,xi_3,degree_1,degree_2,degree_3,intercept_1,intercept_2,intercept_3);
  beta = Rcpp::as<Eigen::VectorXd>(pars["beta"]);
  sigma = pars["sigma"];

  xis_1 = Rcpp::List::create();
  xis_2 = Rcpp::List::create();
  xis_3 = Rcpp::List::create();
}

void TriEBARS::_initial() {
  // initial x_1 knots
  Eigen::VectorXi idx_1 = Rcpp::as<Eigen::VectorXi>(Rcpp::sample(n_1, k_1)).array()-1;
  Eigen::VectorXi idx_rem_1 = Eigen::VectorXi::LinSpaced(n_1, 0, n_1-1);
  xi_1 = Eigen::VectorXd::Zero(k_1);
  for(int i=0;i<k_1;i++) {
    xi_1(i) = knots_1(idx_1(i));
    idx_rem_1(idx_1(i)) = -1; // remark selected knots as -1
  }

  remain_knots_1 = Eigen::VectorXd::Zero(n_1-k_1);
  int j_1 = 0;
  for(int i=0;i<n_1;i++) {
    if(idx_rem_1(i)>=0) {
      remain_knots_1(j_1++) = knots_1(idx_rem_1(i));
    }
  }

  // initial x_2 knots
  Eigen::VectorXi idx_2 = Rcpp::as<Eigen::VectorXi>(Rcpp::sample(n_2, k_2)).array()-1;
  Eigen::VectorXi idx_rem_2 = Eigen::VectorXi::LinSpaced(n_2, 0, n_2-1);
  xi_2 = Eigen::VectorXd::Zero(k_2);
  for(int i=0;i<k_2;i++) {
    xi_2(i) = knots_2(idx_2(i));
    idx_rem_2(idx_2(i)) = -1; // remark selected knots as -1
  }

  remain_knots_2 = Eigen::VectorXd::Zero(n_2-k_2);
  int j_2 = 0;
  for(int i=0;i<n_2;i++) {
    if(idx_rem_2(i)>=0) {
      remain_knots_2(j_2++) = knots_2(idx_rem_2(i));
    }
  }

  // initial x_3 knots
  Eigen::VectorXi idx_3 = Rcpp::as<Eigen::VectorXi>(Rcpp::sample(n_3, k_3)).array()-1;
  Eigen::VectorXi idx_rem_3 = Eigen::VectorXi::LinSpaced(n_3, 0, n_3-1);
  xi_3 = Eigen::VectorXd::Zero(k_3);
  for(int i=0;i<k_3;i++) {
    xi_3(i) = knots_3(idx_3(i));
    idx_rem_3(idx_3(i)) = -1; // remark selected knots as -1
  }

  remain_knots_3 = Eigen::VectorXd::Zero(n_3-k_3);
  int j_3 = 0;
  for(int i=0;i<n_3;i++) {
    if(idx_rem_3(i)>=0) {
      remain_knots_3(j_3++) = knots_3(idx_rem_3(i));
    }
  }
}

bool TriEBARS::_jump_1() {
  double type = Rcpp::runif(1, 0.0, 1.0)[0];
  double birth = _birth_1();
  double death = _death_1();
  if(k_1==0 && type>=birth) {return true;}
  if(fix_k_1) {birth = 0.0;death = 0.0;}

  Eigen::VectorXd xi_new, remain_new;
  int k_new;

  // birth scheme
  if(type < birth) {
    k_new = k_1 + 1;
    int idx = Rcpp::sample(n_1-k_1,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n_1-k_new);
    xi_new.head(k_1) = xi_1; xi_new(k_1) = remain_knots_1(idx);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
    remain_new.head(idx) = remain_knots_1.head(idx);
    remain_new.tail(n_1-k_new-idx) = remain_knots_1.tail(n_1-k_new-idx);
  }

  // death scheme
  else if(type < (birth+death)) {
    k_new = k_1 - 1;
    int idx = Rcpp::sample(k_1,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n_1-k_new);
    xi_new.head(idx) = xi_1.head(idx);
    xi_new.tail(k_new-idx) = xi_1.tail(k_new-idx);
    remain_new.head(n_1-k_1) = remain_knots_1; remain_new(n_1-k_1) = xi_1(idx);
  }

  // relocate scheme
  else {
    k_new = k_1;
    int idx_0 = Rcpp::sample(k_1,1)[0] - 1;
    int idx_1 = Rcpp::sample(n_1-k_1,1)[0] - 1;
    xi_new = xi_1; remain_new = remain_knots_1;
    xi_new(idx_0) = remain_knots_1(idx_1);
    remain_new(idx_1) = xi_1(idx_0);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
  }

  // compute MLE
  Rcpp::List pars_new = cube_spline_regression(t,y,xi_new,xi_2,xi_3,degree_1,degree_2,degree_3,intercept_1,intercept_2,intercept_3);
  Eigen::VectorXd beta_new = Rcpp::as<Eigen::VectorXd>(pars_new["beta"]);
  double sigma_new = pars_new["sigma"];

  // compute marginal likelihood ratio: m^(((k+p-1)(l+p-1)-(k'+p-1)(l'+p-1))/2)*(RSS/RSS')^(m/2)
  double like_ratio = std::pow(m, (k_1-k_new)*(k_2+degree_2+intercept_2)*(k_3+degree_3+intercept_3)/2.0) * std::pow(sigma/sigma_new, m);

  // compute accept probability
  double acc_prob = like_ratio;

  // decide whether to move
  double acc = Rcpp::runif(1, 0.0, 1.0)[0];
  if(acc < acc_prob) {
    k_1 = k_new; xi_1 = xi_new; remain_knots_1 = remain_new;
    beta = beta_new; sigma = sigma_new;
    return true;
  } else {
    return false;
  }
}

bool TriEBARS::_jump_2() {
  double type = Rcpp::runif(1, 0.0, 1.0)[0];
  double birth = _birth_2();
  double death = _death_2();
  if(k_2==0 && type>=birth) {return true;}
  if(fix_k_2) {birth = 0.0; death = 0.0;}

  Eigen::VectorXd xi_new, remain_new;
  int k_new;

  // birth scheme
  if(type < birth) {
    k_new = k_2 + 1;
    int idx = Rcpp::sample(n_2-k_2,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n_2-k_new);
    xi_new.head(k_2) = xi_2; xi_new(k_2) = remain_knots_2(idx);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
    remain_new.head(idx) = remain_knots_2.head(idx);
    remain_new.tail(n_2-k_new-idx) = remain_knots_2.tail(n_2-k_new-idx);
  }

  // death scheme
  else if(type < (birth+death)) {
    k_new = k_2 - 1;
    int idx = Rcpp::sample(k_2,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n_2-k_new);
    xi_new.head(idx) = xi_2.head(idx);
    xi_new.tail(k_new-idx) = xi_2.tail(k_new-idx);
    remain_new.head(n_2-k_2) = remain_knots_2; remain_new(n_2-k_2) = xi_2(idx);
  }

  // relocate scheme
  else {
    k_new = k_2;
    int idx_0 = Rcpp::sample(k_2,1)[0] - 1;
    int idx_2 = Rcpp::sample(n_2-k_2,1)[0] - 1;
    xi_new = xi_2; remain_new = remain_knots_2;
    xi_new(idx_0) = remain_knots_2(idx_2);
    remain_new(idx_2) = xi_2(idx_0);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
  }

  // compute MLE
  Rcpp::List pars_new = cube_spline_regression(t,y,xi_1,xi_new,xi_3,degree_1,degree_2,degree_3,intercept_1,intercept_2,intercept_3);
  Eigen::VectorXd beta_new = Rcpp::as<Eigen::VectorXd>(pars_new["beta"]);
  double sigma_new = pars_new["sigma"];

  // compute marginal likelihood ratio: m^(((k+p-1)(l+p-1)-(k'+p-1)(l'+p-1))/2)*(RSS/RSS')^(m/2)
  double like_ratio = std::pow(m, (k_1+degree_1+intercept_1)*(k_3+degree_3+intercept_3)*(k_2-k_new)/2.0) * std::pow(sigma/sigma_new, m);

  // compute accept probability
  double acc_prob = like_ratio;

  // decide whether to move
  double acc = Rcpp::runif(1, 0.0, 1.0)[0];
  if(acc < acc_prob) {
    k_2 = k_new; xi_2 = xi_new; remain_knots_2 = remain_new;
    beta = beta_new; sigma = sigma_new;
    return true;
  } else {
    return false;
  }
}


bool TriEBARS::_jump_3() {
  double type = Rcpp::runif(1, 0.0, 1.0)[0];
  double birth = _birth_3();
  double death = _death_3();
  if(k_3==0 && type>=birth) {return true;}
  if(fix_k_3) {birth = 0.0; death = 0.0;}

  Eigen::VectorXd xi_new, remain_new;
  int k_new;

  // birth scheme
  if(type < birth) {
    k_new = k_3 + 1;
    int idx = Rcpp::sample(n_3-k_3,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n_3-k_new);
    xi_new.head(k_3) = xi_3; xi_new(k_3) = remain_knots_3(idx);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
    remain_new.head(idx) = remain_knots_3.head(idx);
    remain_new.tail(n_3-k_new-idx) = remain_knots_3.tail(n_3-k_new-idx);
  }

  // death scheme
  else if(type < (birth+death)) {
    k_new = k_3 - 1;
    int idx = Rcpp::sample(k_3,1)[0] - 1;
    xi_new = Eigen::VectorXd::Zero(k_new);
    remain_new = Eigen::VectorXd::Zero(n_3-k_new);
    xi_new.head(idx) = xi_3.head(idx);
    xi_new.tail(k_new-idx) = xi_3.tail(k_new-idx);
    remain_new.head(n_3-k_3) = remain_knots_3; remain_new(n_3-k_3) = xi_3(idx);
  }

  // relocate scheme
  else {
    k_new = k_3;
    int idx_0 = Rcpp::sample(k_3,1)[0] - 1;
    int idx_3 = Rcpp::sample(n_3-k_3,1)[0] - 1;
    xi_new = xi_3; remain_new = remain_knots_3;
    xi_new(idx_0) = remain_knots_3(idx_3);
    remain_new(idx_3) = xi_3(idx_0);
    std::sort(xi_new.data(),xi_new.data()+xi_new.size());
  }

  // compute MLE
  Rcpp::List pars_new = cube_spline_regression(t,y,xi_1,xi_2,xi_new,degree_1,degree_2,degree_3,intercept_1,intercept_2,intercept_3);
  Eigen::VectorXd beta_new = Rcpp::as<Eigen::VectorXd>(pars_new["beta"]);
  double sigma_new = pars_new["sigma"];

  // compute marginal likelihood ratio: m^(((k+p-1)(l+p-1)-(k'+p-1)(l'+p-1))/2)*(RSS/RSS')^(m/2)
  double like_ratio = std::pow(m, (k_1+degree_1+intercept_1)*(k_2+degree_2+intercept_2)*(k_3-k_new)/2.0) * std::pow(sigma/sigma_new, m);

  // compute accept probability
  double acc_prob = like_ratio;

  // decide whether to move
  double acc = Rcpp::runif(1, 0.0, 1.0)[0];
  if(acc < acc_prob) {
    k_3 = k_new; xi_3 = xi_new; remain_knots_3 = remain_new;
    beta = beta_new; sigma = sigma_new;
    return true;
  } else {
    return false;
  }
}

void TriEBARS::_update() {
  bool state = false;
  int max_iter = 100;
  for(int i=0;i<max_iter;i++) {
    // randomly update x_1, x_2 or x_3
    double direc = Rcpp::runif(1,0.0,1.0)[0]; // which direction to jump
    if(direc < 0.333) {
      state = _jump_1();
    } else if(direc < 0.666) {
      state = _jump_2();
    } else {
      state = _jump_3();
    }
    if(state) {
      break;
    }
  }
}

void TriEBARS::rjmcmc(int burns, int steps) {
  for(int i=0;i<burns;i++) {
    _update();
  }
  for(int i=0;i<steps;i++) {
    _update();
    xis_1.push_back(xi_1.array()*(xmax(0)-xmin(0)) + xmin(0));
    xis_2.push_back(xi_2.array()*(xmax(1)-xmin(1)) + xmin(1));
    xis_3.push_back(xi_3.array()*(xmax(2)-xmin(2)) + xmin(2));
  }
}

Rcpp::List TriEBARS::get_knots() {
  return Rcpp::List::create(Rcpp::Named("xi_1")=(xi_1.array()*(xmax(0)-xmin(0)) + xmin(0)),
                            Rcpp::Named("xi_2")=(xi_2.array()*(xmax(1)-xmin(1)) + xmin(1)),
                            Rcpp::Named("xi_3")=(xi_3.array()*(xmax(2)-xmin(2)) + xmin(2)));
}

Rcpp::List TriEBARS::get_samples() {
  return Rcpp::List::create(Rcpp::Named("xi_1")=xis_1, Rcpp::Named("xi_2")=xis_2, Rcpp::Named("xi_3")=xis_3);
}

Eigen::VectorXd TriEBARS::predict(const Eigen::MatrixXd & x_new) {
  // transform x_new to t_new
  Eigen::MatrixXd t_new = (x_new.rowwise() - xmin).array().rowwise() / (xmax-xmin).array();
  // predict
  Eigen::VectorXd y_new = cube_spline_predict(t_new, xi_1, xi_2, xi_3, beta, degree_1, degree_2, degree_3, intercept_1, intercept_2, intercept_3);
  return y_new;
}


// expose Rcpp class
RCPP_MODULE(class_TriEBARS) {
  using namespace Rcpp;

  class_<TriEBARS>("TriEBARS")

    .constructor<Eigen::MatrixXd,Eigen::VectorXd,NumericVector,IntegerVector,List>("constructor")

    .method("rjmcmc", &TriEBARS::rjmcmc, "reversible jump MCMC")
    .method("predict", &TriEBARS::predict, "predict by cube spline regression with EBARS")
    .method("knots", &TriEBARS::get_knots, "return estimated knots")
    .method("samples", &TriEBARS::get_samples, "return posterior samples")
    ;
}

RCPP_EXPOSED_CLASS(TriEBARS)
*/



