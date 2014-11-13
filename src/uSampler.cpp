#include <Rcpp.h>
#include "mnmcem.h"

using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

double min(double a, double b) {
  if(a < b)
    return a;
  return b;
}

bool inC0(NumericMatrix kU, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ, double kC0) {
  if(logRatioCpp(kU, kY, kBeta, kLambda, kX, kZ) < kC0)
    return TRUE;
  return FALSE;
}
  
double logAccept(NumericMatrix current, NumericMatrix proposed, double kC0, bool kYinC0, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ) {
  if (kYinC0)
    return kC0 - logRatioCpp(current, kY, kBeta, kLambda, kX, kZ);
  return min(0, kC0 - logRatioCpp(current, kY, kBeta, kLambda, kX, kZ) + logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ));
}

// [[Rcpp::export]]
NumericMatrix uSamplerCpp(NumericMatrix kU, IntegerVector kY, NumericMatrix kBeta, NumericVector kLambda, 
NumericMatrix kX, NumericMatrix kZ, int B, double kC0) {
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
    // std::cout<<inC0(current, kY, kBeta, kLambda, kX, kZ, kC0)<<'\n';
    
    int tmp = 0;
    int ctr = 0;
    while (tmp == 0) {
      for(int j = 0; j < kC - 1; j++)
        proposed(_, j) = rnorm(kR, 0, sqrt(kLambda(j)));
      if(log(R::runif(0, 1)) < logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ) - kC0)
        tmp = 1;
      ctr = ctr + 1;
      if(ctr == 10000) {
        std::cout<<"Lowering bound.\n";
        kC0 = kC0 - 0.1 * abs(kC0);
        ctr = 0;
      }
    }
    
    // std::cout<<logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ)<<'\n';
    // std::cout<<inC0(current, kY, kBeta, kLambda, kX, kZ, kC0)<<'\n';
    if (inC0(current, kY, kBeta, kLambda, kX, kZ, kC0)) {
      for (int ii = 0; ii < kR; ii++)
        for (int jj = 0; jj < kC - 1; jj++)
          current(ii, jj) = proposed(ii, jj);
    } else {
      if (log(R::runif(0, 1)) < logAccept(current, proposed, kC0, inC0(proposed, kY, kBeta, kLambda, kX, kZ, kC0), kY, kBeta, kLambda, kX, kZ)) {
          for (int ii = 0; ii < kR; ii++)
            for (int jj = 0; jj < kC - 1; jj++)
              current(ii, jj) = proposed(ii, jj);
      }
    }
    
    for (int j = 0; j < kR; j++)
      for (int k = 0; k < kC - 1; k++)
        sample(kR * i + j, k) = current(j, k);
    
    /*
    if(inC0(current, kY, kBeta, kLambda, kX, kZ, kC0)) {
      int tmp = 0;
      while (tmp == 0) {
        for(int j = 0; j < kC - 1; j++)
          proposed(_, j) = rnorm(kR, 0, sqrt(kLambda(j)));
        if(log(R::runif(0, 1)) < logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ) - kC0) {
          for (int ii = 0; ii < kR; ii++)
            for (int jj = 0; jj < kC - 1; jj++)
              current(ii,jj) = proposed(ii,jj);
          tmp = 1;
        }
      }
    } else {
      for(int j = 0; j < kC - 1; j++)
        proposed(_, j) = rnorm(kR, 0, sqrt(kLambda(j)));
        if(log(R::runif(0, 1)) < logAccept(current, proposed, kC0, inC0(proposed, kY, kBeta, kLambda, kX, kZ, kC0), kY, kBeta, kLambda, kX, kZ)) 
          for (int ii = 0; ii < kR; ii++)
            for (int jj = 0; jj < kC - 1; jj++)
              current(ii,jj) = proposed(ii,jj);
    }
    for (int j = 0; j < kR; j++)
      for (int k = 0; k < kC - 1; k++)
        sample(kR * i + j, k) = current(j, k);
  */
  }
  
  // std::cout<<log(R::runif(0, 1))<<'\n';
  // double t = R::rnorm(10,0.5);
  // std::cout<<t<<'\n';
  // sample(_,1) = rnorm(B * kR);
  // logAccept(current, proposed, kC0, kYinC0, kY, kBeta, kLambda, kX, kZ)
  return sample;
  // return proposed;
}