#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;


// [[Rcpp::export]]
double innerProduct(NumericVector x, NumericVector y) {
   return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}


struct InnerProduct : public Worker
{   
   // source vectors
   const RVector<double> x;
   const RVector<double> y;
   
   // product that I have accumulated
   double product;
   
   // constructors
   InnerProduct(const NumericVector x, const NumericVector y) 
      : x(x), y(y), product(0) {}
   InnerProduct(const InnerProduct& innerProduct, Split) 
      : x(innerProduct.x), y(innerProduct.y), product(0) {}
   
   // process just the elements of the range I've been asked to
   void operator()(std::size_t begin, std::size_t end) {
      product += std::inner_product(x.begin() + begin, 
                                    x.begin() + end, 
                                    y.begin() + begin, 
                                    0.0);
   }
   
   // join my value with that of another InnerProduct
   void join(const InnerProduct& rhs) { 
     product += rhs.product; 
   }
};

// [[Rcpp::export]]
double parallelInnerProduct(NumericVector x, NumericVector y) {
   
   // declare the InnerProduct instance that takes a pointer to the vector data
   InnerProduct innerProduct(x, y);
   
   // call paralleReduce to start the work
   parallelReduce(0, x.length(), innerProduct);
   
   // return the computed product
   return innerProduct.product;
}


