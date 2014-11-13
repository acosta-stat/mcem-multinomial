mcem.its <- function(kY, kX, kZ, it = 10, inc = 0.05) {
  nObs <- length(kY)
  kC <- length(unique(kY))
  kP <- ncol(kX)
  kR <- ncol(kZ)
  kBeta.init <- matrix(0, kP, kC - 1)
  kLambda.init <- rep(1, kC - 1)
}