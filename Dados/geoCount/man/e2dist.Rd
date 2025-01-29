\name{e2dist}
\alias{e2dist}

\title{Calculate Distances between Transformed Residuals and Standard Normal}

\usage{
e2dist(e.tran)
}

\arguments{
  \item{e.tran}{a vector which indicates the transformed residuals}
}

\description{
  This function calcualtes three types of distance between the empirical distribution of transformed residuals and standard normal.}

\details{
  \code{HellingerDist} and \code{KolmogorovDist} functions in \code{{distrEx}} are used to compute the distances. See \code{?HellingerDist} and \code{?KolmogorovDist} for details about how the distances are computed.
} 
\value{
 A vector with length 3 containing "Discrete Hellinger", "Smooth Hellinger" and "Kolmogorov" distances.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
require(distrEx)
e2dist(rnorm(200))       }
}

\seealso{
\code{\link{tranR}}, \code{\link{baseline.dist}}, \code{\link{pOne}}.
}

\keyword{Transformed Residual Checking}
