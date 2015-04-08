\name{transitionProbability1D}
\alias{transitionProbability1D}
\title{Compute transition probability using Boltzmann distribution.}
\usage{
   transitionProbability1D(ikBT, x, xflip, J, H, probSel)
}
\arguments{
  \item{ikBT}{1/kB*T (Boltzmann factor)}
  \item{x}{1D Spin sites on the lattice.}
  \item{xflip}{1D Spin sites on the lattice: after a flip.}
  \item{J}{Interaction strength}
  \item{H}{External field}
  \item{probSel}{Which transition probability to use. 1 for Metropolis 2 for Glauber}
}
\value{
  Returns transition probability.
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, and an other flip sites, return the transition probability, 
  applying periodic boundary conditions, i.e., cyclic. 
  This function calls the C function 'transitionProbability1D'.
}
\examples{
  n             <- 10             # 10 spin sites
  mySites       <- genConfig1D(n) # Generate sites
  mySitesNew    <- flipConfig1D(mySites)
  # only short-range part
  transitionProbability1D(1.0, mySites, mySitesNew, 1.0, 0.0, 1) # Metropolis
  transitionProbability1D(1.0, mySites, mySitesNew, 1.0, 0.0, 2) # Glauber
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
