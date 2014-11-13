#include <Rcpp.h>
#include "mnmcem.h"

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


double min0(double a, double b) {
  if(a < b)
    return a;
  return b;
}
  
double logAcceptRW(NumericMatrix current, NumericMatrix proposed, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ) {
  return min0(0, loglikelihoodCpp(proposed, kY, kBeta, kLambda, kX, kZ) - loglikelihoodCpp(current, kY, kBeta, kLambda, kX, kZ));
}


// [[Rcpp::export]]
NumericMatrix uSamplerRWCpp(NumericMatrix kU, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ, int B, double sd) {
  RNGScope scope;
  
  int nObs = kY.size();
  int kP = kX.ncol();
  int kR = kZ.ncol();
  int kC = kBeta.ncol() + 1;

  NumericMatrix sample(B * kR, kC - 1);
  NumericMatrix current(kR, kC - 1);
  NumericMatrix proposed(kR, kC - 1);
  
  for (int i = 0; i < kR; i++)
    for (int j = 0; j < kC - 1; j++)
      current(i, j) = sample(i, j) = kU(i, j); // rnorm(1);
  
  for (int i = 1; i < B; i++) {
    for(int j = 0; j < kC - 1; j++)
      proposed(_, j) = rnorm(kR, 0, sd);
    
    for (int ii = 0; ii < kR; ii++)
      for (int jj = 0; jj < kC - 1; jj++)
        proposed(ii, jj) = proposed(ii, jj) + current(ii, jj);
    
    // std::cout<<logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ)<<'\n';
    if (log(R::runif(0, 1)) < logAcceptRW(current, proposed, kY, kBeta, kLambda, kX, kZ)) {
      for (int ii = 0; ii < kR; ii++)
        for (int jj = 0; jj < kC - 1; jj++)
          current(ii, jj) = proposed(ii, jj);
    }
    
    
    for (int j = 0; j < kR; j++)
      for (int k = 0; k < kC - 1; k++)
        sample(kR * i + j, k) = current(j, k);
  }
  
  return sample;
  // return proposed;
}