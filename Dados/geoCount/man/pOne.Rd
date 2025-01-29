\name{pOne}
\alias{pOne}

\title{Calculate One-side P-value}

\usage{
pOne(d.obs, d.base)
}

\arguments{
  \item{d.obs}{a value (or a vector) which indicates the distance for observed data}
  \item{d.base}{a vector (or a matrix) which indicates the samples of baseline distances}
}

\description{
  This function calculates one-side p-value(s) for observed distance(s) with respect to the samples of baseline distances.}


\value{
 A p-value (or a vector of p-values).
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# For single value
pOne(2, rnorm(10000)) 
# Visualize with plot_pRPS()
plot_pRPS(2, rnorm(10000), nm="d")
# For vector
pOne(1:3, matrix(rnorm(30000),,3))
}
}

\seealso{
\code{\link{e2dist}}, \code{\link{baseline.dist}}, \code{\link{pRPS}}, \code{\link{plot_pRPS}}.
}

\keyword{Transformed Residual Checking}
