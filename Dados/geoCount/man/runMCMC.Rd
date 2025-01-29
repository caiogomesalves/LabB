\name{runMCMC}
\alias{runMCMC}

\title{Perform Robust MCMC Algorithms for GLSM}

\usage{
runMCMC(Y, L = 0, loc, X = NULL, run = 200, run.S = 1, 
        rho.family = "rhoPowerExp", Y.family = "Poisson", 
        priorSigma = "Halft", parSigma = c(1, 1), ifkappa = 0, 
        scales = c(0.5, 1.65^2 + 0.8, 0.8, 0.7, 0.15), 
        phi.bound = c(0.005, 1), 
        initials = list(c(1), 1.5, 0.2, 1), 
        MCMCinput = NULL, partial = FALSE, famT = 1)
}

\arguments{
  \item{Y}{a vector of length n which indicates the response variables}
  \item{L}{a vector of length n; it indicates the time duration during which the Poisson counts are accumulated, or the total number of trials for Binomial response; if 0 is found in the vector, 1 will be used to replace all the values in the vector}
  \item{loc}{a \eqn{n \times 2} matrix which indicates the coordinates of locations}
  \item{X}{a \eqn{n \times p} covariate matrix; the default value "NULL" indicates no covariate}
  \item{run}{the number of iterations}
  \item{run.S}{the number of internal iterations for latern variables}
  \item{rho.family}{take the value of \code{"rhoPowerExp"}, \code{"rhoMatern"}, or \code{"rhoSph"} which indicates the powered exponential, Matern, or Spherical correlation function is used}
  \item{Y.family}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}
  \item{priorSigma}{the prior distribution for \eqn{\sigma}, the options include "Halft" (positive-truncated t distribution), "InvGamma" (inverse gamma distribution), and "Reciprocal" (reciprocal distribution)}
  \item{parSigma}{the parameters for the prior distribution of \eqn{\sigma}: when priorSigma = "Halft" the first parameter is scale and the second is degree of freedom; when priorSigma = "InvGamma" the first parameter is shape and the second is scale; when priorSigma = "Reciprocal" both parameters are ignored}
  \item{ifkappa}{take zero or non-zero value which indicates whether \eqn{\kappa} should be sampled}
  \item{scales}{a vector which indicates the tuning parameters for \eqn{(S, \beta, \sigma,\phi,\kappa)} respectively}
  \item{phi.bound}{the upper and lower bound for \eqn{\phi} }
  \item{initials}{a list which indicates the initial values for \eqn{(\beta, \sigma,\phi,\kappa)} respectively}
  \item{MCMCinput}{a list of alternative settings; usually the result from \code{MCMCinput} function}
  \item{partial}{a logical input which indicates whether partial posterior sampling should be used; only works for Y.family = "Poisson"}
  \item{famT}{take the value of 1, 2, or 3 which indicates the type of partial posterior sampling: 1 means "mean" diagnostic statistic is used, 2 means "maximum", and 3 means "minimum"; ignored if \code{partial=FALSE}}
}

\description{
  This function performs robust MCMC algorithms for generalized linear spatial models and generates posterior samples for latent variables and hyper-parameters.
}

\details{
  Group updating scheme, Langevin algorithms, and Data-based parameterization are applied to improve the robustness and efficiency of MCMC algorithms. The flat priors are used to guarantee an appropriate posterior. See my dissertation for more details.
  
  During each iteration of Gibbs sampling process, the group of latent variables is updated "run.S" times to improve accuracy and reduce autocorrelations.
}

\value{
 A list with elements:
  \item{S.posterior}{a \eqn{n \times run} matrix containing the posterior samples for latent variables}
  \item{m.posterior}{a \eqn{(p+1) \times run} matrix (in case of p covariate variables) or a vector with length "run" (no covariate case), containing the posterior samples for \eqn{\beta}}
  \item{s.posterior}{a vector with length "run" containing the posterior samples for \eqn{\sigma}}
  \item{a.posterior}{a vector with length "run" containing the posterior samples for \eqn{\phi}}
  \item{k.posterior}{a vector with length "run" containing the posterior samples for \eqn{\kappa} in the case that "ifkappa" is set to non-zero value}
  \item{AccRate}{a vector which indicates the acceptance rates}
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  data(Weed)
  Y <- Weed[,3]
  loc <- unifLoc(Weed[,1:2])
  L <- rep(1, length(Y))
  input.Weed <- MCMCinput( run=1000, run.S=10, rho.family="rhoPowerExp", 
          Y.family = "Poisson", 
          priorSigma = "Halft", parSigma = c(1, 1),
          ifkappa=0,
          scales=c(0.5, 3.5, 0.9, 0.6, 0.5), 
          phi.bound=c(0.005, 1), 
          initials=list(c(-1), 1, 0.1, 1) )
  res <- runMCMC(Y, L=L, loc=loc, X=NULL, MCMCinput = input.Weed )          
}
}

\seealso{
   \code{\link{MCMCinput}}, \code{\link{runMCMC.sf}}.
}

\keyword{MCMC}
