\name{isPerform1D}
\alias{isPerform1D}
\title{Perform metropolis MC on 1D Ising model}
\usage{
  isPerform1D(ikBT, x, J, H, nstep, ensembleM, probSel)
}
\arguments{
  \item{ikBT}{1/kB*T (Boltzmann factor)}
  \item{x}{1D Spin sites on the lattice.}
  \item{J}{Interaction strength}
  \item{H}{External field}
  \item{nstep}{Number of MC steps requested}
  \item{ensembleM}{Value of the theoretical magnetization (could be thermodynamic limit value)}
  \item{probSel}{Which transition probability to use. 1 for Metropolis 2 for Glauber}
}
\value{
  Returns a pair list containing values for omegaM, Fluctuating metric vector for Magnetisation (length of naccept),
  naccept, number of MC steps accepted and  nreject, number of MC steps rejected and times as accepted time steps.
  Times corresponds to times where flips occur, this is so-called transition times ('metropolis time' or 'single flip time') 
  to judge the timings between two accepted steps.
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, and an other flip sites, perform Metropolis Monte Carlo
  applying periodic boundary conditions, i.e., cyclic. 
  This function calls the C function 'isPerform1D'.
}
\examples{
  n        <- 10 # 10 spin sites
  mySites  <- genConfig1D(n) # Generate sites
  output   <- isPerform1D(1.0, mySites, 1.0, 0.0, 10, 0.5, 1) # Metropolis
  output   <- isPerform1D(1.0, mySites, 1.0, 0.0, 10, 0.5, 2) # Glauber
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
