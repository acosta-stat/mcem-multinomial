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
double qFunctionCpp(NumericMatrix kU, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ) {
  double value = 0;
  // double counter = 0;
  int nObs =  kY.size();
  int kP = kX.ncol();
  int kR = kZ.ncol();
  int kC = kBeta.ncol() + 1;
  int kK = kU.nrow()/kR;
  
  for (int i = 0; i < nObs; i++) {
    // First triple sum, we can just sum on the non-zero yijh directly if y belongs in the Cth class then this sum is 0
    for (int m = 0; m < kK; m++) {
      double wijh0 = 0;
      if (kY(i) < kC) {
        for (int j = 0; j < kP; j++) {
          wijh0 = wijh0 + kX(i, j) * kBeta(j, kY(i) - 1);
        }
        for (int j = 0; j < kR; j++) {
          wijh0 = wijh0 + kZ(i, j) * kU(m * kR +j, kY(i) - 1);
        }
        value = value + wijh0/kK;
      }
    }
     
    // Sum with the log
    for (int m = 0; m < kK; m++) {
      double wijh1 = 0;
      for (int h = 0; h < kC - 1; h++) {
        double wijh0 = 0;
        for (int j = 0; j < kP; j++) {
          wijh0 = wijh0 + kX(i, j) * kBeta(j, h);
        }
        for (int j = 0; j < kR; j++) {
          wijh0 = wijh0 + kZ(i, j) * kU(m * kR + j, h);
        }
        wijh1 = wijh1 + exp(wijh0);
      }
      value = value - log(1 + wijh1)/kK;
    }
    
    // Last sum (normal u)
    for (int m = 0; m < kK; m++) {
      double tmp0 = 0;
      for (int j = 0; j < kC - 1; j++) {
        for (int h = 0; h < kR; h++) {
          tmp0 = tmp0 - (0.5 * log(kLambda(j)) - 0.5 / kLambda(j) * kU(m * kR +h, j) * kU(m * kR +h, j))/kK;
        }
      }
    }
  }
  return value;
}