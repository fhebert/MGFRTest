//[[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;
using Rcpp::as;

double logLnbeta(Eigen::MatrixXd X, Eigen::VectorXd Y, Eigen::VectorXd beta){
  Eigen::MatrixXd Xbeta = X*beta;
  double res = (Y.array()*Xbeta.array()-(1+Xbeta.array().exp()).log()).sum();
  return(res);
}

//[[Rcpp::export]]
Eigen::VectorXd GLMProbs(Eigen::MatrixXd X, Eigen::VectorXd Y){
  int n = X.rows();
  int p = X.cols();
  Eigen::VectorXd beta_init = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd pbeta = Eigen::VectorXd::Zero(n);
  pbeta = pbeta.array()+0.5;
  Eigen::MatrixXd PP = Eigen::MatrixXd::Zero(n,p);
  PP = PP.array()+0.25;
  Eigen::MatrixXd XP = X.array()*PP.array();
  Eigen::MatrixXd A = (XP.transpose()*X).inverse();
  Eigen::VectorXd new_beta = beta_init+A*(X.transpose()*(Y-pbeta));
  double d0 = -logLnbeta(X,Y,beta_init);
  double d1 = -logLnbeta(X,Y,new_beta);
  double d = std::abs(d1-d0)/(std::abs(d1)+0.1);
  while(d>0.00000001){
    beta_init = new_beta;
    pbeta = (1+(-X*beta_init).array().exp()).inverse();
    PP -= PP;
    PP = PP.array().colwise()+pbeta.array()*(1-pbeta.array());
    XP = X.array()*PP.array();
    A.noalias() = (XP.transpose()*X).inverse();
    new_beta.noalias() += A*(X.transpose()*(Y-pbeta));
    d0 = -logLnbeta(X,Y,beta_init);
    d1 = -logLnbeta(X,Y,new_beta);
    d = std::abs(d1-d0)/(std::abs(d1)+0.1);
  }
  pbeta = (1+(-X*new_beta).array().exp()).inverse();
  return(pbeta);
}

//[[Rcpp::export]]
Eigen::MatrixXd PermGLMProbs(Eigen::MatrixXd X, Eigen::MatrixXd Y0){
  int n = X.rows();
  int N = Y0.cols();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n,N);
  for(int i=0;i<N;i++){
    res.col(i) = GLMProbs(X,Y0.col(i));
  }
  return(res);
}