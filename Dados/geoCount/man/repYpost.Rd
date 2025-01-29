\name{repYpost}
\alias{repYpost}

\title{Generate Replicated Data with Posterior Samples of Latent Variables}

\usage{
repYpost(res.m, L, Y.family="Poisson")
}

\arguments{
  \item{res.m}{a list with elements containing the posterior samples of latent variables and parameters for observed locations}
  \item{L}{a vector of length n; it indicates the time duration during which the Poisson counts are accumulated, or the total number of trials for Binomial response}
  \item{Y.family}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}
}

\description{
  This function generates replicated data sets based on posterior samples of latent variables.}


\value{
 A matrix of replicated data sets.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep.post <- repYpost(res.m, L) }
}

\seealso{
\code{\link{repYeb}}
}

\keyword{Data Simulation}
