// The script is from https://gist.github.com/brandonwillard/710a7bab8526db394bad
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iterator> 
#include <algorithm> 
#include <iostream> 

using namespace Rcpp;

struct SProd : public RcppParallel::Worker {
    // source matrix
    const arma::sp_fmat& spA;
    const arma::sp_fmat& spB;
    const std::set<int>& anrows;
    const std::set<int>& bnrows;
    //std::set<int>::iterator anrows_iter;
    //std::set<int>::iterator bnrows_iter;
    // destination matrix
    arma::sp_fmat& output;
    // initialize with source and destination
    SProd(const arma::sp_fmat& spA, 
          const arma::sp_fmat& spB,  
          const std::set<int>& anrows,
          const std::set<int>& bnrows,
          arma::sp_fmat& output) 
    : spA(spA), spB(spB), anrows(anrows), bnrows(bnrows), output(output) 
      //, anrows_iter(std::unique(spA.row_indices, spA.row_indices + spA.n_nonzero))
      //, bnrows_iter(std::unique(spB.row_indices, spB.row_indices + spB.n_nonzero))
    {
    }
    // take the square root of the range of elements requested
    void operator()(const std::size_t begin, const std::size_t end) {
//      std::cout << "begin=" << begin << ", end=" << end << std::endl << std::flush;
      std::set<int>::const_iterator b_start = bnrows.begin();
      std::advance(b_start, begin);
      std::set<int>::const_iterator b_end = bnrows.begin();
      std::advance(b_end, end+1);
      for (std::set<int>::const_iterator b_it = b_start; b_it != b_end; b_it++) {
        //std::cout << "b_row=" << *b_it << std::endl << std::flush;
//        if (Progress::check_abort())
//          break; //return(R_NilValue);
        for (std::set<int>::const_iterator a_it = anrows.begin(); a_it != anrows.end(); a_it++) {
//          if (Progress::check_abort())
//            break; //return(R_NilValue);
          const int k = (*b_it) * spA.n_rows + (*a_it);
          //std::cout << "k=" << k << std::endl << std::flush;
          output.row(k) = spA.row(*a_it) % spB.row(*b_it);
          //p.increment();
        }
      }
    }
  };


// [[Rcpp::export]]
int sparse_row_idx_mult_v2(arma::sp_fmat& spA){
	int a = spA.n_rows;
	return(a);
}


// [[Rcpp::export]]
arma::fmat sparse_row_idx_mult(arma::sp_fmat& spA, arma::sp_fmat & spB) {
//    Progress p(0, false);

    int a = spA.n_rows;	
    std::cout << "spA" << a << std::endl;	

//    std::cout << "sparse_row_idx_mult OK1" << std::endl;
//    std::cout << "spA.n_rows " << spA.n_rows << std::endl;
//    std::cout << "spA.n_nonzero " << spA.n_nonzero << std::endl;
//    std::cout << "sparse_row_idx_mult OK1" << std::endl;
//    arma::uvec rowind = spA.row_indices;
    std::cout << "sparse_row_idx_mult OK000" << std::endl;
    std::set<int> anrows(spA.row_indices, spA.row_indices + spA.n_nonzero);
    std::cout << "sparse_row_idx_mult OK2" << std::endl;
    std::set<int> bnrows(spB.row_indices, spB.row_indices + spB.n_nonzero);
    std::cout << "sparse_row_idx_mult OK3" << std::endl;
    const int bnrows_len = bnrows.size();
    const int res_rows = (spB.n_rows) * (spA.n_rows);
    arma::sp_fmat C(res_rows, spA.n_cols);
    std::cout << "B_nonzero_rows=" << bnrows_len << std::endl;
    std::cout << "res_rows=" << res_rows << std::endl;
//    SProd prodWorker(spA, spB, anrows, bnrows, p, C);
    SProd prodWorker(spA, spB, anrows, bnrows, C);
    RcppParallel::parallelFor(0, bnrows_len-1, prodWorker, 1000);
//    Rcpp::S4 Cout(Rcpp::wrap(C));
    arma::fmat Cout(C); 	
    return(Cout);
}
