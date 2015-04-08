# 
#  Ising Model Dynamics MC
#    R wrapper functions to C library functions
#  (c) 2013 by Dr.Mehmet Suzen
#  GPLv3 or higher
#

# Generate Random Configuration in 1D
#    Arguments  
#    n     : number of sites 
genConfig1D <- function(n) {
  x   <- rep(0.0,n)
  out <- .C("genConfig1D", n=as.integer(n), x=as.double(x));
  out$x
}

# Random single flip configuration in 1D
#   Arguments
#   x  : sites vector
#   Return random flip single site
flipConfig1D <- function(x) {
  out <- .C("flipConfig1D", n=as.integer(length(x)), x=as.double(x));
  out$x
}

# Random many flip configuration in 1D
#   Arguments
#   x  : sites
#   upperF : Number of repeats for random flip
flipConfig1Dmany <- function(x, upperF) {
  out <- .C("flipConfig1Dmany", n=as.integer(length(x)), x=as.double(x), upperF=as.integer(upperF));
  out$x
}

#  Nearest-Neighbour energy in periodic boundary conditions in 1D
#    Arguments  
#    ll     : number of sites 
#    vec    : site spin configuration 
#    energy : energy to be returned 
lattice1DenergyNN <- function(x) {
  energy <-0
  if(!is.numeric(x))
    stop("argument x must be numeric");
  out <- .C("lattice1DenergyNN", n=as.integer(length(x)), vec=as.double(x), energy=as.double(energy));
  out$energy
}

#
#  Sum given vector
#    Arguments
#    x    : site spin configuration 
sumVec <- function(x) {
  sum <-0
  if(!is.numeric(x))
    stop("argument x must be numeric");
  out <- .C("sumVec", n=as.integer(length(x)), vec=as.double(x), sum=as.double(sum));
  out$sum
}

# Total Energy in 1D
#   x : site spin configuration
#   J : interaction strength
#   H : external field
totalEnergy1D <- function(x, J, H) {
  energy <- 0
  if(!is.numeric(x)) 
    stop("argument x must be numeric");
  out <- .C("totalEnergy1D", n=as.integer(length(x)), vec=as.double(x), J=as.double(J), H=as.double(H), totalEnergy=as.double(energy));
  out$totalEnergy
}

# Compute transition probability
#  Arguments
#   ikBT : 1/kB*T (boltzmann factor)
#   x    : 1D Spin sites on the lattice
#   xflip: 1D Spin sites on the lattice: after a flip
#   J    : interaction strength
#   H    : external field
#   probSel : which transition probability to use 1 Metropolis 2 Glauber
transitionProbability1D <- function(ikBT, x, xflip, J, H, probSel) {
   prob <- 0.0;
  if(!is.numeric(x)) 
    stop("argument x must be numeric");
  if(!is.numeric(xflip)) 
    stop("argument xflip must be numeric");
  out <- .C("transitionProbability1D", ikBT=as.double(ikBT), n=as.integer(length(x)), vec=as.double(x), vecFlip=as.double(xflip), 
                                       J=as.double(J), H=as.double(H), prob=as.double(prob), probSel=as.integer(probSel));
  out$prob
}


# Take One MC Step : importance sampling!
#   Author : msuzen
#   Arguments
#   ikBT    : 1 over Temperature times Boltzmann constant
#   x    : 1D Spin sites on the lattice
#   J    : interaction strength
#   H    : external field
#   rout : Pair list, flip states (vec) and if step is accepted (accept)
#   probSel : which transition probability to use 1 Metropolis 2 Glauber
isStep1D <- function(ikBT, x, J, H, probSel) {
   prob <- 0.0;
   accept <- 0.0;
  if(!is.numeric(x)) 
    stop("argument x must be numeric");
    out <- .C("isStep1D", ikBT=as.double(ikBT), n=as.integer(length(x)), vec=as.double(x),  
                          J=as.double(J), H=as.double(H), prob=as.double(prob), accept=as.integer(accept), 
                          probSel=as.integer(probSel));
  rout <- list(vec=out$vec, accept=out$accept)
  rout
}

#  Perform importance sampling MC on 1D
#
#   Author : msuzen
#   Arguments
#   ikBT    : 1 over Temperature times Boltzmann constant
#   x    : 1D Spin sites on the lattice
#   J    : interaction strength
#   H    : external field
#   nstep     : number of MC steps requested
#   ensembleM : Value of the theoretical magnetization (could be thermodynamic limit value)
# Output pair list:
#   omegaM    : Fluctuating metric vector for Magnetisation (length of naccept)
#   naccept   : number of MC steps accepted
#   nreject   : number of MC steps rejected
#   probSel : which transition probability to use 1 Metropolis 2 Glauber
#
isPerform1D <- function(ikBT, x, J, H, nstep, ensembleM, probSel) {
  if(!is.numeric(x)) 
    stop("argument x must be numeric");
    omegaM  <- rep(0.0, nstep); # Let R allocate memory; assume all steps will be accepted
                                # index 1 to naccept will be reported in the vector
    times   <- rep(1.0, nstep); # Let R allocate memory: this is the so called time to accepted flip!
    naccept <- 0;
    nreject <- 0;
    out     <- .C("isPerform1D", ikBT=as.double(ikBT), n=as.integer(length(x)), vec=as.double(x),  
                             J=as.double(J), H=as.double(H), ensembleM=as.double(ensembleM), 
                             omegaM=as.double(omegaM), nstep=as.integer(nstep), naccept=as.integer(naccept),
                             nreject=as.integer(nreject), times=as.integer(times), probSel=as.integer(probSel));
  rout <- list(omegaM=out$omegaM, nreject=out$nreject, naccept=out$naccept, times=out$times);
  rout
}
