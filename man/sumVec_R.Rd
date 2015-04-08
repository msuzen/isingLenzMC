\name{sumVec_R}
\alias{sumVec_R}
\title{Sum given vector}
\usage{
  sumVec_R(x)
}
\arguments{
  \item{x}{1D Spin sites on the lattice}
}
\value{
  Returns the sum, corresponding the long-range part.
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, return the sum.  This function calls the C function 'sumVec'.
}
\examples{
  n          <- 10             # 10 spin sites
  mySites    <- genConfig1D_R(n) # Generate sites
  sumVecs    <- sumVec_R(mySites)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
