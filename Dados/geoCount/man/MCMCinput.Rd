\name{MCMCinput}
\alias{MCMCinput}

\title{Settings for the MCMC Algorithm}

\usage{
MCMCinput(run = 200, run.S = 1, rho.family = "rhoPowerExp", 
        Y.family = "Poisson", 
        priorSigma = "Halft", parSigma = c(1, 1),
        ifkappa = 0, 
        scales = c(0.5, 1.65^2 + 0.8, 0.8, 0.7, 0.15), 
        phi.bound = c(0.005, 1), 
        initials = list(c(1), 1.5, 0.2, 1))
}

\arguments{
  \item{run}{the number of iterations}
  \item{run.S}{the number of internal iterations for latent variables}
  \item{rho.family}{take the value of \code{"rhoPowerExp"}, \code{"rhoMatern"}, or \code{"rhoSph"} which indicates the powered exponential, Matern, or Spherical correlation function is used}
  \item{Y.family}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}
  \item{priorSigma}{the prior distribution for \eqn{\sigma}, the options include "Halft" (positive-truncated t distribution), "InvGamma" (inverse gamma distribution), and "Reciprocal" (reciprocal distribution)}
  \item{parSigma}{the parameters for the prior distribution of \eqn{\sigma}: when priorSigma = "Halft" the first parameter is scale and the second is degree of freedom; when priorSigma = "InvGamma" the first parameter is shape and the second is scale; when priorSigma = "Reciprocal" both parameters are ignored}
  \item{ifkappa}{take zero or non-zero value which indicates whether \eqn{\kappa} should be sampled}
  \item{scales}{a vector which indicates the tuning parameters for \eqn{(S, \beta, \sigma,\phi,\kappa)} respectively}
  \item{phi.bound}{the upper and lower bound for \eqn{\phi} }
  \item{initials}{a list which indicates the initial values for \eqn{(\beta, \sigma,\phi,\kappa)} respectively}
}

\description{
  This function sets up the parameters and initial values used for the MCMC algorithms.
}

\details{
  During each iteration of Gibbs sampling process, the group of latent variables is updated "run.S" times to improve accuracy and reduce autocorrelations.
}

\value{
 A list of setting parameters.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  input <- MCMCinput( run = 10000, run.S = 10, 
          rho.family = "rhoPowerExp", 
          Y.family = "Poisson", 
          priorSigma = "Halft", parSigma = c(1, 1),
          ifkappa=0,
          scales=c(0.5, 1.5, 0.9, 0.6, 0.5), 
          phi.bound=c(0.005, 1), 
          initials=list(c(-1, 2, 1), 1, 0.1, 1) )
  res <- runMCMC(Y, L=0, loc=loc, X=loc, MCMCinput = input )
}
}

\seealso{
   \code{\link{runMCMC}}, \code{\link{runMCMC.sf}}.
}

\keyword{MCMC}
