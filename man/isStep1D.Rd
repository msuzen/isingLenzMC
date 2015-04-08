\name{isStep1D}
\alias{isStep1D}
\title{Carry one step Metropolis Monte Carlo on 1D ising model}
\usage{
  isStep1D(ikBT, x, J, H, probSel)
}
\arguments{
  \item{ikBT}{1/kB*T (Boltzmann factor)}
  \item{x}{1D Spin sites on the lattice.}
  \item{J}{Interaction strength}
  \item{H}{External field}
  \item{probSel}{Which transition probability to use. 1 for Metropolis 2 for Glauber}
}
\value{
  A pair list, flip states (vec) and if step is accepted (accept).
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively and the usual thermodynamic parameters ikBt, J and H.
  Perform 1 step metropolis Monte Carlo, applying periodic boundary conditions, 
  i.e., cyclic. This function calls the C function 'isStep1D'. Importance sampling 
  is applied.
}
\examples{
  n             <- 10 # 10 spin sites
  mySites       <- genConfig1D(n) # Generate sites
  # only short-range part
  isStep1D(1.0, mySites, 1.0, 0.0, 1) # Metropolis
  isStep1D(1.0, mySites, 1.0, 0.0, 2) # Glauber
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
