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
double qFunction2(NumericMatrix kU, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ) {
  int kR = kZ.ncol();
  int kK = kU.nrow()/kR;
  double value = 0;
  for (int m = 0; m < kK; m++) {
    value = value + loglikelihoodCpp(kU, kY, kBeta, kLambda, kX, kZ);
  }
  return(value);
}