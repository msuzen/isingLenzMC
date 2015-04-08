\name{transferMatrix}
\alias{transferMatrix}
\title{Compute theoretical transfer matrix}
\usage{
  transferMatrix(ikBt, J, H)
}
\arguments{
  \item{ikBt}{1/kB*T (Boltzmann factor)}
  \item{J}{Interaction strength}
  \item{H}{External field}
}
\value{
  Returns transfer matrix and its eigenvalues in a pair list.
}
\description{
  Compute transfer matrix
}
\examples{
  transferMatrix(1.0, 1.0, 0)
}
\author{
  Mehmet Suzen <mehmet.suzen@physics.org>
}
