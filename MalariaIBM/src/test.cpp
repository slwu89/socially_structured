#include <MalariaIBM/test.h>
#include <Rcpp.h>

namespace MalariaIBM {

  test::test(const int &N_new){
    N = N_new;
    std::cout << "i'm an object of class test being born at memory location: " << this << std::endl;
  };

  test::~test(){
    std::cout << "i'm an object of class test getting killed at memory location: " << this << std::endl;
  };

  int test::get_N(){
    return(N);
  };

  void test::set_N(const int &N_new){
    N = N_new;
  };

  void test::get_memLoc(){
    std::cout << "i'm an object of class test at memory location: " << this << std::endl;
  };

}
