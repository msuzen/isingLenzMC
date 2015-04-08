\name{flipConfig1D_R}
\alias{flipConfig1D_R}
\title{Given Flip a site randomly}
\usage{
  flipConfig1D_R(x)
}
\arguments{
  \item{x}{1D Spin sites on the lattice}
}
\value{
  Returns vector that contains 1s or -1s. 
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, flip any of the site randomly. The function uses default 
  RNG (Marsienne-Twister) unless changed by the user, within R, to generate 
  a vector that contains 1s or -1s. This function is a pure R implementation
}
\examples{
  n          <- 10 # 10 spin sites
  mySites    <- genConfig1D_R(n) # Generate sites
  # now flip 
  mySitesNew <- flipConfig1D_R(mySites)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
