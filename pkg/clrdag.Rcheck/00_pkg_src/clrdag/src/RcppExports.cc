#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// DC_ADMM
List DC_ADMM(mat X, mat A, mat Lambda, umat D, double tau, double mu, double rho, double tol_abs, double tol_rel, int dc_max_iter, int admm_max_iter, int auto_init, int trace_obj);
RcppExport SEXP _clrdag_DC_ADMM(SEXP XSEXP, SEXP ASEXP, SEXP LambdaSEXP, SEXP DSEXP, SEXP tauSEXP, SEXP muSEXP, SEXP rhoSEXP, SEXP tol_absSEXP, SEXP tol_relSEXP, SEXP dc_max_iterSEXP, SEXP admm_max_iterSEXP, SEXP auto_initSEXP, SEXP trace_objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< umat >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type tol_abs(tol_absSEXP);
    Rcpp::traits::input_parameter< double >::type tol_rel(tol_relSEXP);
    Rcpp::traits::input_parameter< int >::type dc_max_iter(dc_max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type admm_max_iter(admm_max_iterSEXP);
    Rcpp::traits::input_parameter< int >::type auto_init(auto_initSEXP);
    Rcpp::traits::input_parameter< int >::type trace_obj(trace_objSEXP);
    rcpp_result_gen = Rcpp::wrap(DC_ADMM(X, A, Lambda, D, tau, mu, rho, tol_abs, tol_rel, dc_max_iter, admm_max_iter, auto_init, trace_obj));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_clrdag_DC_ADMM", (DL_FUNC) &_clrdag_DC_ADMM, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_clrdag(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
