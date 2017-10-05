// test class

#ifndef _MALARIAIBM_TEST_H_
#define _MALARIAIBM_TEST_H_

#include <Rcpp.h>

namespace MalariaIBM {

class test {
// public members
public:
  test(const int &N_new);
  ~test();

  int get_N();
  void set_N(const int &N_new);

  void get_memLoc();

// private members
private:

  int N;

};



}

#endif
