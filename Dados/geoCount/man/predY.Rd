\name{predY}
\alias{predY}

\title{Predict for Unsampled Locations}

\usage{
predY(res.m, loc, locp, X = NULL, Xp = NULL, Lp = 0, k = 1, 
      rho.family = "rhoPowerExp", Y.family = "Poisson", 
      parallel = NULL, n.cores = getOption("cores"), 
      cluster.type = "SOCK")
}

\arguments{
  \item{res.m}{a list with elements containing the posterior samples of latent variables and parameters for observed locations}
  \item{loc}{a matrix which indicates the coordinates of the observed locations}
  \item{locp}{a matrix which indicates the coordinates of the predicting locations}
  \item{X}{the covariate matrix for observed locations}
  \item{Xp}{the covariate matrix for predicting locations}
  \item{Lp}{a vector which indicates the time duration during which the Poisson counts are accumulated or the total number of trials for Binomial response; if 0 is found in the vector, 1 will be used to replace all the values in the vector}
  \item{k}{a value for fixed \eqn{\kappa}; ignored if there are posterior samples for \eqn{\kappa} in "res.m"}
  \item{rho.family}{take the value of \code{"rhoPowerExp"} or \code{"rhoMatern"} which indicates the powered exponential or Matern correlation function is used}
  \item{Y.family}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}
  \item{parallel}{the default input \code{NULL} indicates no parallel computing will be applied; any input value indicates parallel computing with the help of \code{{snowfall}} }
  \item{n.cores}{the number of CPUs that will be used for parallel computing; used only if \code{parallel} isn't \code{NULL}}
  \item{cluster.type}{type of cluster to be used for parallel computing; can be "SOCK", "MPI", "PVM", or "NWS"; used only if \code{parallel="snowfall"}}
}

\description{
  This function generates posterior predictive samples of latent and response variables for predicting locations.}

\details{
  This function performs parallel computing with the help of \code{{snowfall}} package.
}

\value{
 A list with elements:
 \item{latent.predict}{a matrix containing the posterior predictive samples for latent variables}
 \item{Y.predict}{a matrix containing the posterior predictive samples for response variables}
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Ypred <- predY(res.m, loc, locp, X=loc, Xp=locp, k=1, 
               rho.family = "rhoPowerExp", Y.family = "Poisson")
# require(snowfall)
# Ypred <- predY(res.m, loc, locp, X=loc, Xp=locp, 
#                parallel="snowfall", n.cores = 4)
Ypred.avg <- rowMeans(Ypred$Y); EYpred.avg <- rowMeans(exp(Ypred$latent))
}
}

\seealso{
   \code{\link{runMCMC}}.
}

\keyword{MCMC}
