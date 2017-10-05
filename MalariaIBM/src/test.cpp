#include <MalariaIBM/test.h>
#include <Rcpp.h>

namespace MalariaIBM {

  // inline definition of constructor to accept default argument values
  test::test(const int &N_new){
    N = N_new;
    Rcpp::Rcout << "i'm an object of class test being born at memory location: " << this << std::endl;
  };

  test::~test(){
    Rcpp::Rcout << "i'm an object of class test getting killed at memory location: " << this << std::endl;
  };

  int test::get_N(){
    return(N);
  };

  // return number of emerging adults
  void test::set_N(const int &N_new){
    N = N_new;
  };

  void test::get_memLoc(){
    Rcpp::Rcout << "i'm an object of class test at memory location: " << this << std::endl;
  };

}
