/*
  DEBUG

    get the bugs out!

    write me
    write me
    write me
    write me
*/

#include <Rcpp.h>

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

  std::cout << "  Number of processors available = " << omp_get_num_procs ( ) << std::endl;
  std::cout << "  Number of threads =              " << omp_get_max_threads ( ) << std::endl;

  wtime = omp_get_wtime ( );

  #pragma omp parallel private(id)
  {
      id = omp_get_thread_num ( );
      std::cout << "  This is process " << id << std::endl;
  }

    wtime = omp_get_wtime ( ) - wtime;

    std::cout << "  Normal end of execution." << std::endl;
    std::cout << "  Elapsed wall clock time = " << wtime << std::endl;
};
