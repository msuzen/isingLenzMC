# 
#  Ising Model MC functions
#    R functions: That maps the C code for measuring the performance
#  (c) 2013 by Dr.Mehmet Suzen
#  GPLv3 or higher
#


genConfig1D_R <- function(n) {
 mapply(genUniform, 1:n)
}

genUniform <- function(n=1) {
  x <- 0.0;
  a<-runif(1);
  if(a>=0.5) {
    x <-1.0
  } else {
    x <--1.0
  }
  x
}

flipConfig1D_R <- function(x) {
  rn    <-  round(runif(1)*length(x))
  x[rn] <- -1*x[rn]
  x 
}

lattice1DenergyNN_R <- function(x) {
  n      <- length(x);
  energy <- sum(x[2:n] * x[1:(n-1)])+x[n]*x[1]
  energy
}

totalEnergy1D_R <- function(x, J, H) {
  energy <- J*lattice1DenergyNN_R(x) + H*sumVec_R(x)
  energy
}

transitionProbability1D_R <- function(ikBT, x, xFlip, J, H) {
  energyOrg  <- totalEnergy1D_R(x, J, H)
  energyFlip <- totalEnergy1D_R(xFlip, J, H)
  DeltaE     <- energyOrg-energyFlip
  prob       <- exp(-ikBT*DeltaE)/(1 + exp(-ikBT*DeltaE))
  prob
}


sumVec_R <- function(x) {
  sum(x);
}
