\name{plot_pRPS}
\alias{plot_pRPS}

\title{Plot Observed vs. Reference Diagnostic Statistics}

\usage{
plot_pRPS(T.obs, T.rep, nm = "x")
}

\arguments{
  \item{T.obs}{a value which indicates the observed statistic}
  \item{T.rep}{a vector which indicates the reference statistics}
  \item{nm}{the name of the diagnostic statistics}
}

\description{
  This function plots the observed statistic vs. the empirical density of reference statistics.}

\details{
  \code{density} function is used to compute the empirical density of reference statistics.
} 
\value{
 A plot.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
plot_pRPS(1, rnorm(1000)) }
}

\seealso{
\code{\link{BMCT}}, \code{\link{pRPS}}.
}

\keyword{Bayesian Model Checking}
