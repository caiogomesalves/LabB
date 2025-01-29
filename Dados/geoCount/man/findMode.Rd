\name{findMode}
\alias{findMode}

\title{Estimate Mode of the Posterior Samples}

\usage{
findMode(x, ...)
}

\arguments{
  \item{x}{a vector of posterior samples}
  \item{...}{other parameters used when estimating then empirical density function; see \code{?density}}
}

\description{
  This function estimates the mode of empirical density function for a given posterior samples.}

\details{
  This function uses \code{density} function to estimate the empirical density function.
} 
\value{
 The value of mode.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
findMode(rnorm(1000))
}
}


\keyword{MCMC}
