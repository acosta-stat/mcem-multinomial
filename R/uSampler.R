uSampler <- function(kU, kY, kBeta, kLambda, kX, kZ, B, r = 0.01) {
  nObs =  length(kY)
  kP = ncol(kX)
  kR = ncol(kZ)
  kC = ncol(kBeta) + 1
  
  c0 <- sum(table(kY) * log(table(kY)/length(kY))) + r * abs(sum(table(kY) * log(table(kY)/length(kY))))
  # print(c0)
  
  # Some functions
  in.C0 <- function(x) {
    return(logRatioCpp(x, kY, kBeta, kLambda, kX, kZ) < c0)
  }
  
  logAccept <- function(current, proposed, yinC0) {
    if (yinC0) {
      return(c0 - logRatioCpp(current, kY, kBeta, kLambda, kX, kZ))
    } else {
      return(min(1, c0 - logRatioCpp(current, kY, kBeta, kLambda, kX, kZ) + logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ)))
    }
  }
  
  usample <- matrix(0, B * kR, kC - 1)
  current <- kU
  usample[1:kR, ] <- current
  for (i in 2:B) {
    # print(in.C0(current))
    if (in.C0(current)) {
      tmp <- 0
      while (tmp == 0) {
        proposed <- matrix(rnorm(kR * (kC - 1), 0, sqrt(kLambda)), kR, kC - 1, byrow = TRUE)
        # print(logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ) - c0)
        if(log(runif(1)) < logRatioCpp(proposed, kY, kBeta, kLambda, kX, kZ) - c0) {
          current <- proposed
          tmp <- 1
        }
      }
    } else {
      proposed <- matrix(rnorm(kR * (kC - 1), 0, sqrt(kLambda)), kR, kC - 1, byrow = TRUE)
      if(log(runif(1)) < logAccept(current, proposed, in.C0(proposed))) {
        current <- proposed
      }
    }
    
    # print((kR * (i - 1) + 1):(kR * i))
    # print(in.C0(proposed))
    usample[(kR * (i - 1) + 1):(kR * i), ] <- current
  }
  return(usample)
}