\name{flipConfig1Dmany}
\alias{flipConfig1Dmany}
\title{Flip a single site randomly many times}
\usage{
  flipConfig1Dmany(x, upperF)
}
\arguments{
  \item{x}{1D spin sites on the lattice.}
  \item{upperF}{The number of times}
}
\value{
  Returns vector that contains 1s or -1s. 
}
\description{
  Given a vector of flip sites, 1s or -1s, representing up and down spins 
  respectively, flip any of the site randomly, repeat it many times. 
  The function uses default RNG (Marsienne-Twister) unless changed by the user, 
  within R, to generate a vector that contains 1s or -1s. This function calls 
  'flipConfig1Dmany' C function.
}
\examples{
  n          <- 10 # 10 spin sites
  mySites    <- genConfig1D(n) # Generate sites
  # now flip 100 times
  mySitesNew <- flipConfig1Dmany(mySites, 100)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
