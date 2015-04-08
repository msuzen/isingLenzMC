\name{lattice1DenergyNN_R}
\alias{lattice1DenergyNN_R}
\title{Nearest-Neighbour energy in periodic boundary conditions in 1D}
\usage{
  lattice1DenergyNN_R(x)
}
\arguments{
  \item{x}{1D Spin sites on the lattice}
}
\value{
  Returns the nearest neighbour energy.
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, return nearest neighbour energy, applying periodic boundary
  conditions, i.e., cyclic. This function is a pure R implementation.
}
\examples{
  n          <- 10 # 10 spin sites
  mySites    <- genConfig1D_R(n) # Generate sites
  nnEnergy <- lattice1DenergyNN(mySites)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
