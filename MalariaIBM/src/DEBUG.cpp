/*
  DEBUG

    get the bugs out!

    write me
    write me
    write me
    write me
*/

#include <Rcpp.h>
#include <MalariaIBM/test.h>

// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

//' debug openmp!
//'
//' get the bugs out!
//'
//' @export
// [[Rcpp::export]]
void DEBUG_OPENMP(){

  int id;
  double wtime;

  std::cout << "  Number of processors available = " << omp_get_num_procs() << std::endl;
  std::cout << "  Number of threads =              " << omp_get_max_threads() << std::endl;

  wtime = omp_get_wtime();

  #pragma omp parallel private(id)
  {
      id = omp_get_thread_num();
      std::cout << "  This is process " << id << std::endl;
  }

    wtime = omp_get_wtime() - wtime;

    std::cout << "  Normal end of execution." << std::endl;
    std::cout << "  Elapsed wall clock time = " << wtime << std::endl;
};


//' debug openmp with classes
//'
//' get the bugs out!
//'
//' @export
// [[Rcpp::export]]
void DEBUG_CLASS_OPENMP(const int &N){
  std::vector<MalariaIBM::test> tests;
  for(size_t i=0; i<N; i++){
    tests.push_back(MalariaIBM::test((int)(i)));
  }

  int id;
  size_t i;
  #pragma omp parallel for private(i,id)
  for(i=0; i<N; i++){
    std::cout << "running on thread: " << omp_get_thread_num() << std::endl;
    tests[i].get_memLoc();
  }

};

//' debug openmp with classes
//'
//' get the bugs out!
//' check speed in serial
//'
//' @export
// [[Rcpp::export]]
void DEBUG_CLASS_SERIAL(const int &N){
  std::vector<MalariaIBM::test> tests;
  for(size_t i=0; i<N; i++){
    tests.push_back(MalariaIBM::test((int)(i)));
  }

  size_t i;
  for(i=0; i<N; i++){
    std::cout << "i: " << i << std::endl;
    tests[i].get_memLoc();
  }

};
