\name{pRPS}
\alias{pRPS}

\title{Calculate P-value and RPS}

\usage{
pRPS(T.obs, T.rep)
}

\arguments{
  \item{T.obs}{a value which indicates the observed statistic}
  \item{T.rep}{a vector which indicates the reference statistics}
}

\description{
  This function calculates p-value and relative predictive surprise (RPS) by comparing observed and reference statistics.}

\details{
  \code{density} function is used to compute the empirical density of reference statistics.
} 
\value{
 A vector of p-value and RPS.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
pRPS(2, rnorm(1000)) }
}

\seealso{
\code{\link{BMCT}}, \code{\link{plot_pRPS}}.
}

\keyword{Bayesian Model Checking}
