#
#  Ergodic Dynamics of Ising Model 
#  R simulation functions, data generation, analysis and other utilities
#  (c) 2013, 2014 by Dr.Mehmet Suzen
#  GPLv3 or higher
#

require("isingLenzMC"); #

#
# Using Analytic solution for ensemble magnetisation 
# using transfer matrix
#
#  Two state ising model: general case
#  http://micro.stanford.edu/~caiwei/me334/
#  Ch 7 page 35
# 
#Z(h, B, J, N) := exp(N*B*J) * ((cosh(B*h) + sqrt(sinh(B*h)^2+exp(-4*B*J)))^N + (cosh(B*h) - sqrt(sinh(B*h)^2+exp(-4*B*J)))^N);

# now the magnetisation (mulitply with 1/B) compute with maxima
# very long!
#diff(log(Z(h, B, J, N)), h);
#
#(%i4) diff(log(Z(h, B, J, N)), h);
#           B cosh(h B) sinh(h B)
#(%o4) ((---------------------------- + B sinh(h B))
#               - 4 B J       2
#        sqrt(%e        + sinh (h B))
#         - 4 B J       2                   N - 1
# (sqrt(%e        + sinh (h B)) + cosh(h B))      N
#                     B cosh(h B) sinh(h B)
# + (B sinh(h B) - ----------------------------)
#                         - 4 B J       2
#                  sqrt(%e        + sinh (h B))
#                     - 4 B J       2       N - 1
# (cosh(h B) - sqrt(%e        + sinh (h B)))      N)
#          - 4 B J       2                   N
#/((sqrt(%e        + sinh (h B)) + cosh(h B))
#                       - 4 B J       2       N
# + (cosh(h B) - sqrt(%e        + sinh (h B))) )

getEnsembleMagnetisation <- function(n, ikBT, J, H) {
# Derivative of Equation 35 of Wu notes at Standford
# Maxima output above
 if(n < 700) { # This is ad-hoc assumption is that we recover n -> inf solution above this
   cHB   <- cosh(H*ikBT)
   sHB   <- sinh(H*ikBT)
   e4B   <- exp(-4*ikBT*J)
   sq4B  <- sqrt(e4B+sHB^2)
   isq4B <- ikBT*cHB*sHB/sq4B
   # Per site ensemble averaged magnetisation
   eMag  <- isq4B + ikBT*sHB
   eMag  <- n * eMag * (sq4B+cHB)^(n-1)
   eMag  <- eMag + n*(ikBT*sHB - isq4B) * (cHB - sq4B)^(n-1)
   eMag  <- eMag * ((sq4B+cHB)^(n) + (cHB-sq4B)^(n))^(-1)
   if(is.nan(eMag)) return(0.9934344)
   return(eMag/ikBT/n)
  } else {
   return(0.9934344)
 }
}

#
# This is to perform 1D Ising simulation and Getting TM metric
# Thirumalai - Mountain Fluctuation Metric
#
magnetisationMetricTM1D <- function(ikBT, J, H, n, nstep, transP) {
  ensembleM <- getEnsembleMagnetisation(n, ikBT, J, H) # ensemble average magnetisation per-site
  x         <- genConfig1D(n)
  mySim     <- isPerform1D(ikBT, x, J, H, nstep, ensembleM, transP) 
  time      <- 1:mySim$naccept
  # Note that, Metropolis variable names! could be glauber as well depending transP value (1 or 2)
  magnetisationMetric <- mySim$omegaM[time]
  metropolisTime      <- mySim$times[time]
  out                 <- list(magnetisationMetric=magnetisationMetric, metropolisTime=metropolisTime);
  out
}


#
# Running a single sim (time vs. magnetisation metric and 
# generate data file with the parametrisation as fileName
#
runSim <- function(ikBT, J, H, n, nstep, transP) {
    if(transP == 1) dynamicsName <- 'metropolis';
    if(transP == 2) dynamicsName <- 'glauber';
    dataFileName <- paste(dynamicsName, "TmMetric-n=", sprintf('%.3e', n),
                          "ikBT=",  sprintf('%.3f', ikBT),
                          "H=",  sprintf('%.3f', H), ".dat",
                          sep="")
    out          <- magnetisationMetricTM1D(ikBT, J, H, n, nstep, transP)
    # Note that, Metropolis variable names! could be glauber as well depending transP value (1 or 2)
    ll           <- 2*length(out$metropolisTime) 
    timeMag      <- matrix(1:ll, ncol= 2, nrow=length(out$metropolisTime))
    timeMag[,1]  <- out$metropolisTime
    timeMag[,2]  <- out$magnetisationMetric
    write.table(timeMag, file=dataFileName, row.names=FALSE, col.names=FALSE) 
}

#
# Read the data file: coloumns: c("metropolisTime", "magnetisationTM")
readDat <- function(filePath) {
  matrixMC           <- as.matrix(read.csv(filePath, sep=" ", header=0))
  colnames(matrixMC) <- NULL
  return(matrixMC) 
}

#
# Obtain power law exponent alpha     x ~ t^alpha 
# with  Using formulation given by Newman-2005-Contemp. Phys. Vol 46
#   Appendix B
#
getPowerLaw <- function(x) {
  n     <- length(x)
  xMin  <- min(x)
  sumLn <- 0
  for(i in 1:n) {
    sumLn <- sumLn + log(x[i]/xMin)
  }
  alpha <- 1 + n/sumLn
  alpha
}
#
#

