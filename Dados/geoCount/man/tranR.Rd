\name{tranR}
\alias{tranR}

\title{Calculate Transformed Residuals for Observed Data}

\usage{
tranR(Y.obs, Y.rep, discrete = FALSE)
}

\arguments{
  \item{Y.obs}{a vector which indicates the observed data set}
  \item{Y.rep}{a matrix which indicates the reference data sets}
  \item{discrete}{a logical value which indicates if the distribution of response variable is discrete}
}

\description{
  This function approximates transformed residuals for the observed data by using reference data.}


\value{
 A vector of transformed residuals.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1)
tranR(Y.obs, Y.rep) }
}

\seealso{
\code{\link{cdfU}}, \code{\link{plot_etran}}, \code{\link{e2dist}}.
}

\keyword{Transformed Residual Checking}
