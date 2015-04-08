\name{totalEnergy1D_R}
\alias{totalEnergy1D_R}
\title{Total energy in periodic boundary conditions in 1D}
\usage{
  totalEnergy1D_R(x, J, H)
}
\arguments{
  \item{x}{1D Spin sites on the lattice.}
  \item{J}{The strength of interaction.}
  \item{H}{The value of the external field.}
}
\value{
  Return the total energy.
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, return total energy, applying periodic boundary
  conditions, i.e., cyclic. This function is pure R implementation.
}
\examples{
  n             <- 10 # 10 spin sites
  mySites       <- genConfig1D_R(n) # Generate sites
  # only short-range part
  myTotalEnergy <- totalEnergy1D_R(mySites, 1.0, 0.0)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
