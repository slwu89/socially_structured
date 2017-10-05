#include "house.h"
#include <Rcpp.h>

namespace MalariaIBM {

  house::house(const int &id_new, const int &block_id_new){
    id = id_new;
    block_id = block_id_new;
    std::cout << "i'm an object of class house being born at memory location: " << this << std::endl;
  };

  house::~house(){
    humans.clear();
    std::cout << "i'm an object of class house getting killed at memory location: " << this << std::endl;
  };

  int house::get_id(){
    return(id);
  };

  int house::get_block_id(){
    return(block_id);
  };

  void house::add_human(human* h){
    humans.insert(h);
  };

  void house::get_memLoc(){
    std::cout << "i'm an object of class house at memory location: " << this << std::endl;
  };

}
