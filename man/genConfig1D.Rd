\name{genConfig1D}
\alias{genConfig1D}
\title{Generate one dimensional spin sites randomly}
\usage{
  genConfig1D(n)
}
\arguments{
  \item{n}{The number of spin sites on the lattice.}
}
\value{
  Returns vector that contains 1s or -1s. 
}
\description{
  The function uses default RNG (Marsienne-Twister) unless changed by the 
  user, within R, to generate a vector that contains 1 or -1. This reflects
  spin sites. This function calls 'genConfig1D' C function.
}
\examples{
  n <- 10 # 10 spin sites
  genConfig1D(n)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
