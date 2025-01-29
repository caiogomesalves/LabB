\name{baseline.dist}
\alias{baseline.dist}

\title{Generate Distance Samples to Build Baseline Distribution}

\usage{
baseline.dist(n, iter)
}

\arguments{
  \item{n}{the number of residuals}
  \item{iter}{the number of distance samples to generate}
}

\description{
  This function generates the samples of distance to build the baseline distribution for standard normal.}

\details{
  \code{HellingerDist} and \code{KolmogorovDist} functions in \code{{distrEx}} are used to compute the distances. See \code{?HellingerDist} and \code{?KolmogorovDist} for details about how the distances are computed.
} 
\value{
 A \eqn{iter \times 3} matrix for three types of distance: "Discrete Hellinger", "Smooth Hellinger" and "Kolmogorov". 
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# Time-consuming! Run once with large "iter" and 
# save the results for future use
d.base <- baseline.dist(50, iter=100) }
}

\seealso{
\code{\link{d.base}}, \code{\link{plot_baseline}}, \code{\link{pOne}}.
}

\keyword{Transformed Residual Checking}
