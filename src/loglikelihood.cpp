#include <Rcpp.h>
#include "mnmcem.h"

//#include <RcppArmadillo.h>
using namespace Rcpp;
//// [[Rcpp::depends("RcppArmadillo")]]

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

int indicator(int a, int b) {
  if(a == b) {
    return 1;
  } else {
    return 0;
  }
}

// [[Rcpp::export]]
double loglikelihoodCpp(NumericMatrix kU, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ) {
  double value = 0;
  // double counter = 0;
  int nObs =  kY.size();
  int kP = kX.ncol();
  int kR = kZ.ncol();
  int kC = kBeta.ncol() + 1;

  for (int i = 0; i < nObs; i++) {
    // First triple sum, we can just sum on the correct yijh directly
    double wijh0 = 0;
    if (kY(i) < kC) {
      for (int j = 0; j < kP; j++) {
        wijh0 = wijh0 + kX(i, j) * kBeta(j, kY(i) - 1);
      }
      for (int j = 0; j < kR; j++) {
        wijh0 = wijh0 + kZ(i, j) * kU(j, kY(i) - 1);
      }
      value = value + wijh0;
    }

    // Sum with the log
    NumericVector sumex(kC - 1);
    
    for (int h = 0; h < kC - 1; h++) {
      wijh0 = 0;
      for (int j = 0; j < kP; j++) {
        wijh0 = wijh0 + kX(i, j) * kBeta(j, h);
      }
      for (int j = 0; j < kR; j++) {
        wijh0 = wijh0 + kZ(i, j) * kU(j, h);
      }
      sumex(h) = wijh0;
      //wijh1 = wijh1 + exp(wijh0);
    }
    double mex = max(sumex);
    for (int ii = 0; ii < kC - 1; ii++)
      sumex(ii) = sumex(ii) - mex;
    sumex = exp(sumex);
    double wijh1 = exp(mex) * sum(sumex);
    value = value - log(1 + wijh1);
    //std::cout<<log(1 + wijh1)<<'\n';
  }
   
  // Last sum (on the U-s)
  double tmp0 = 0;
  for (int j = 0; j < kC - 1; j++) {
    for (int h = 0; h < kR; h++) {
      tmp0 = tmp0 - 0.5 * log(kLambda(j)) - 0.5 / kLambda(j) * kU(h, j) * kU(h, j);
    }
  }

  return value + tmp0;

  // std::vector<double> wijh(kY.size());
  /*
  NumericVector xx(6);
  for (int i = 0; i < 6; i++)
    xx(i) = xx(i) + R::rnorm(0, 1);
  xx = exp(xx);
  double mxx = max(xx);
  std::cout<<mxx<<'\n';
  double sxx = sum(xx);
  for (int i = 0; i < 6; i++)
    std::cout<<xx(i)<<'\n';
  std::cout<<sxx<<'\n';
  */
  
}