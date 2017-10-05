/*
  I'm a human!
    write me
    write me
    write me
    write me
    write me
*/

#ifndef _MALARIAIBM_HUMAN_H_
#define _MALARIAIBM_HUMAN_H_

#include <Rcpp.h>
#include <iostream>

namespace MalariaIBM {

class human {
// public members
public:
                              human(const int &id_new);
                              ~human();

  int                         get_id();
  std::string                 get_state();
  void                        set_state(const std::string &state_new);

  void                        get_memLoc();

// private members
private:

  int                         id;
  std::string                 state;      // life state

};



}

#endif
