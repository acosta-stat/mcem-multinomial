mcem.its <- function(kY, kX, kZ, it = 10, inc = 0.05, tk = 1000, sd0 = 0.05) {
  nObs <- length(kY)
  kC <- length(unique(kY))
  kP <- ncol(kX)
  kR <- ncol(kZ)
  kBeta.init <- matrix(0, kP, kC - 1)
  kLambda.init <- rep(1, kC - 1)
  kBeta <- matrix(0, kP * it, kC - 1)
  kLambda <- matrix(0, it, kC - 1)
  kBeta[1:kP, ] < kBeta.init
  kLambda[1, ] <- kLambda.init
  uStart <- matrix(0, kR, kC - 1)
  betaToMaxR <- function(kBeta, kU, kY, kLambda, kX, kZ) {
    -betaToMaxCpp(kU, kY, matrix(kBeta, ncol(kX), ncol(kU)), kLambda, kX, kZ)
  }
  for (i in 2:it) {
    uSample0 <- uSamplerRWCpp(as.matrix(uStart), kY, as.matrix(kBeta[(kP * (i - 2) + 1):(kP * (i - 1)), ]), kLambda[i - 1, ], as.matrix(kX), kZ, tk, sd0)
    kLambda[i, ] <- lambdaMaxCpp(uSample0)
    tmp0 <- optim(par = as.vector(kBeta[(kP * (i - 2) + 1):(kP * (i - 1)), ]), fn = betaToMaxR, kU = uSample0, kY = kY, kLambda = kLambda[i, ], kX = kX, kZ = kZ, method = 'BF' )
    if (tmp0$convergence > 0) {
      print(tmp0)
      stop("Maximization went wrong at iteration number ", i, ".")
    }
    kBeta[(kP * (i - 1) + 1):(kP * i), ] <- tmp0$par
    uStart <- uSample0[(tk - kR + 1):tk, ]
    tk <- tk + ceiling(inc * tk)
  }
  return(list(Beta = kBeta, Lambda = kLambda))
}
