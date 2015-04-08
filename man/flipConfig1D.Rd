\name{flipConfig1D}
\alias{flipConfig1D}
\title{Given Flip a site randomly}
\usage{
  flipConfig1D(x)
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
  a vector that contains 1s or -1s. This function calls 'flipConfig1D' C function.
}
\examples{
  n          <- 10 # 10 spin sites
  mySites    <- genConfig1D(n) # Generate sites
  # now flip 
  mySitesNew <- flipConfig1D(mySites)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
