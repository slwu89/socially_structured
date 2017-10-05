/*
  I'm a human!
    write me
    write me
    write me
    write me
    write me
*/

#ifndef _MALARIAIBM_HOUSE_H_
#define _MALARIAIBM_HOUSE_H_

#include <set>
#include <iostream>
#include <Rcpp.h>
#include "human.h"

namespace MalariaIBM {

class house {
// public members
public:
                              house(const int &id_new, const int &block_id_new);
                              ~house();

  int                         get_id();
  int                         get_block_id();

  void                        add_human(human* h);

  void                        get_memLoc();

// private members
private:

  int                         id;         // numeric id of house
  int                         block_id;   // numeric id of block
  std::set<human*>            humans;     // list of people here now


};



}

#endif
