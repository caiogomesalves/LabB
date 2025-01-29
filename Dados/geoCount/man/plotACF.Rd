\name{plotACF}
\alias{plotACF}

\title{Auto-correlation Plot for Latent Variables}

\usage{
plotACF(S.mcmc, lags = NULL, ...)
}

\arguments{
  \item{S.mcmc}{a matrix (or data.frame) with each row containing the posterior samples of one latent variable}
  \item{lags}{the maximum number of lags; the default "NULL" will result in \eqn{10*\log_{10}(N/m)} where N is the number of observations and m the number of series}
  \item{...}{more plotting parameters}
}

\description{
  This function plots auto-correlation curves for latent variables.}

\details{
  This function uses \code{acf} function to compute the estimates of auto-correlation.
} 
\value{
 No return value. A plot of auto-correlation curves.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
plotACF(res.m$S.posterior)
}
}


\keyword{MCMC}
