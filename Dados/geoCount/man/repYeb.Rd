\name{repYeb}
\alias{repYeb}

\title{Generate Replicated Data with Estimated Parameters}

\usage{
repYeb(N.sim, loc, L, X = NULL, rho.family = "rhoPowerExp", 
       Y.family="Poisson", res.m = NULL, est = "mode", 
       beta = NULL, sigma = NULL, phi = NULL, k = 1)
}

\arguments{
  \item{N.sim}{the number of replicated data sets to be simulated}
  \item{loc}{a \eqn{n \times 2} matrix which indicates the coordinates of observed locations}
  \item{L}{a vector of length n; it indicates the time duration during which the Poisson counts are accumulated, or the total number of trials for Binomial response}
  \item{X}{a \eqn{n \times p} covariate matrix; the default value "NULL" indicates no covariate}
  \item{rho.family}{take the value of \code{"rhoPowerExp"} or \code{"rhoMatern"} which indicates the powered exponential or Matern correlation function is used}
  \item{Y.family}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}
  \item{res.m}{a list with elements containing the posterior samples of latent variables and parameters for observed locations}
  \item{est}{take the value of \code{"mode"} which indicates the mode of posterior samples will be used as the parameter estimate; otherwise, the mean will be used}
  \item{beta}{a value which indicates the estimation for \eqn{\beta}; ignored if "res.m" is given}
  \item{sigma}{a value which indicates the estimation for \eqn{\sigma}; ignored if "res.m" is given}
  \item{phi}{a value which indicates the estimation for \eqn{\phi}; ignored if "res.m" is given}
  \item{k}{a value which indicates the estimation for \eqn{\kappa}; ignored if "res.m" is given and contains the posterior samples for \eqn{\kappa}}
}

\description{
  This function generates replicated data sets based on estimated parameters (given or from posterior samples).}


\value{
 A \eqn{n \times N.sim} matrix of replicated data sets.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# Estimate parameters from posterior samples
Yrep.eb <- repYeb(N.sim=2000, loc, L, res.m, est = "mode")
# Pre-determined parameters (also an efficient way to simulate massive data sets)
Yrep.eb <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1, 
                  k = 1)
}
}

\seealso{
\code{\link{repYpost}}, \code{\link{simData}}.
}

\keyword{Data Simulation}
