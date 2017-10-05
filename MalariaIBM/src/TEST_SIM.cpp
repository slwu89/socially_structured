/*
  DEBUG

    get the bugs out!

    write me
    write me
    write me
    write me
*/

#include <Rcpp.h>
#include <iostream>
#include "human.h"
#include "house.h"

using namespace MalariaIBM;

//' make a house
//'
//' get the bugs out! get the people in!
//'
//' @export
// [[Rcpp::export]]
void TEST_HOUSE(int nPeople){
  house house1(1,1);
  std::cout << "made a house, id: " << house1.get_id() << " " << std::endl;
  house1.get_memLoc();

  for(size_t i=1; i<=nPeople; i++){
    human* h = new human((int)(i));
    house1.add_human(h);
  }
};
