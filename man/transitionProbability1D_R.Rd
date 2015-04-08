\name{transitionProbability1D_R}
\alias{transitionProbability1D_R}
\title{Compute transition probability using Boltzmann distribution.}
\usage{
  transitionProbability1D_R(ikBT, x, xFlip, J, H)
}
\arguments{
  \item{ikBT}{1/kB*T (Boltzmann factor)}
  \item{x}{1D Spin sites on the lattice.}
  \item{xFlip}{1D Spin sites on the lattice: after a flip.}
  \item{J}{Interaction strength}
  \item{H}{External field}
}
\value{
  Returns transition probability.
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, and an other flip sites, return the transition probability, 
  applying periodic boundary conditions, i.e., cyclic. 
  This function is pure R implementation.
}
\examples{
  n             <- 10 # 10 spin sites
  mySites       <- genConfig1D_R(n) # Generate sites
  mySitesNew    <- flipConfig1D_R(mySites)
  # only short-range part
  transitionProbability1D_R(1.0, mySites, mySitesNew, 1.0, 0.0)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
