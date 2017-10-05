#include "human.h"
#include <Rcpp.h>

namespace MalariaIBM {

  human::human(const int &id_new){
    id = id_new;
    std::cout << "i'm an object of class human being born at memory location: " << this << std::endl;
  };

  human::~human(){
    std::cout << "i'm an object of class human getting killed at memory location: " << this << std::endl;
  };

  int human::get_id(){
    return(id);
  };

  std::string human::get_state(){
    return(state);
  };

  void human::set_state(const std::string &state_new){
    state = state_new;
  };

  void human::get_memLoc(){
    std::cout << "i'm an object of class human at memory location: " << this << std::endl;
  };

}
