// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/MalariaIBM.h"
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _MalariaIBM_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// test__ctor
MalariaIBM::test test__ctor(int N_new);
RcppExport SEXP _MalariaIBM_test__ctor(SEXP N_newSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N_new(N_newSEXP);
    rcpp_result_gen = Rcpp::wrap(test__ctor(N_new));
    return rcpp_result_gen;
END_RCPP
}
// test__get_N
int test__get_N(MalariaIBM::RcppR6::RcppR6<MalariaIBM::test> obj_);
RcppExport SEXP _MalariaIBM_test__get_N(SEXP obj_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MalariaIBM::RcppR6::RcppR6<MalariaIBM::test> >::type obj_(obj_SEXP);
    rcpp_result_gen = Rcpp::wrap(test__get_N(obj_));
    return rcpp_result_gen;
END_RCPP
}
// test__set_N
void test__set_N(MalariaIBM::RcppR6::RcppR6<MalariaIBM::test> obj_, int N_new);
RcppExport SEXP _MalariaIBM_test__set_N(SEXP obj_SEXP, SEXP N_newSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MalariaIBM::RcppR6::RcppR6<MalariaIBM::test> >::type obj_(obj_SEXP);
    Rcpp::traits::input_parameter< int >::type N_new(N_newSEXP);
    test__set_N(obj_, N_new);
    return R_NilValue;
END_RCPP
}
// test__get_memLoc
void test__get_memLoc(MalariaIBM::RcppR6::RcppR6<MalariaIBM::test> obj_);
RcppExport SEXP _MalariaIBM_test__get_memLoc(SEXP obj_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MalariaIBM::RcppR6::RcppR6<MalariaIBM::test> >::type obj_(obj_SEXP);
    test__get_memLoc(obj_);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MalariaIBM_rcpp_hello_world", (DL_FUNC) &_MalariaIBM_rcpp_hello_world, 0},
    {"_MalariaIBM_test__ctor", (DL_FUNC) &_MalariaIBM_test__ctor, 1},
    {"_MalariaIBM_test__get_N", (DL_FUNC) &_MalariaIBM_test__get_N, 1},
    {"_MalariaIBM_test__set_N", (DL_FUNC) &_MalariaIBM_test__set_N, 2},
    {"_MalariaIBM_test__get_memLoc", (DL_FUNC) &_MalariaIBM_test__get_memLoc, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MalariaIBM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
