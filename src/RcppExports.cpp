// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// GLMProbs
Eigen::VectorXd GLMProbs(Eigen::MatrixXd X, Eigen::VectorXd Y);
RcppExport SEXP _MGFRTest_GLMProbs(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(GLMProbs(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// PermGLMProbs
Eigen::MatrixXd PermGLMProbs(Eigen::MatrixXd X, Eigen::MatrixXd Y0);
RcppExport SEXP _MGFRTest_PermGLMProbs(SEXP XSEXP, SEXP Y0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Y0(Y0SEXP);
    rcpp_result_gen = Rcpp::wrap(PermGLMProbs(X, Y0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MGFRTest_GLMProbs", (DL_FUNC) &_MGFRTest_GLMProbs, 2},
    {"_MGFRTest_PermGLMProbs", (DL_FUNC) &_MGFRTest_PermGLMProbs, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MGFRTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}