#include <Rcpp.h>
#include "mnmcem.h"

//#include <RcppArmadillo.h>
using namespace Rcpp;
//// [[Rcpp::depends("RcppArmadillo")]]

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector lambdaMaxCpp(NumericMatrix kU) {
  int kC = kU.ncol() + 1;
  int kUobs = kU.nrow();
  NumericVector kLambda(kC - 1);
  
   
  for (int h = 0; h < kC - 1; h++) {
    double value = 0;
    for (int i = 0; i < kUobs; i++) {
      value = value + kU(i, h) * kU(i, h);
    }
    kLambda(h) = value/kUobs;
  }
  return kLambda;
}
