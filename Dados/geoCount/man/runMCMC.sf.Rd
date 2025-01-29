\name{runMCMC.sf}
\alias{runMCMC.sf}

\title{Perform Robust MCMC Algorithms for GLSM in Parallel}

\usage{
runMCMC.sf(Y, L=0, loc, X=NULL, 
    run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", 
    priorSigma = "Halft", parSigma = c(1, 1),
    ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1), 
    MCMCinput=NULL, partial = FALSE, famT=1,
    n.chn = 2, n.cores = getOption("cores"), cluster.type="SOCK")
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
  \item{MCMCinput}{a list of alternative settings}
  \item{partial}{a logical input which indicats whether partial posterior sampling should be used; only works for Y.family = "Poisson"}
  \item{famT}{take the value of 1, 2, or 3 which indicates the type of partial posterior sampling: 1 means "mean" diagnostic statistic is used, 2 means "maximum", and 3 means "minimum"; ignored if partial=FALSE}
  \item{n.chn}{the number of Markov chain sets that will be generated in parallel}
  \item{n.cores}{the number of CPUs that will be used to generate parallel Markov chains}
  \item{cluster.type}{type of cluster to be used for parallel computing; can be "SOCK", "MPI", "PVM", or "NWS"}
}

\description{
  This function performs robust MCMC algorithms in a parallel way for generalized linear spatial models and generates posterior samples for latent variables and hyper-parameters.   
}

\details{
  Essentially, this function runs \code{\link{runMCMC}} function simultaneously on different CPUs (if there are more than one CPU available) with different initial values. In the case the number of available CPUs is less than "n.chn", Markov chains will be put in a queue.
  
  This function performs parallel computing with the help of \code{{snow}} and \code{{snowfall}} packages.
}

\value{
 A list of length "n.chn" containing the result of each Markov chain. Each element is a list with elements:
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
  require(snowfall)
  data(datWeed)
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
  res.prl <- runMCMC.sf(Y, L=L, loc=loc, X=NULL, 
          MCMCinput = input.Weed, n.chn = 4, n.cores = 4, cluster.type="SOCK")          
  }
}

\seealso{
   \code{\link{MCMCinput}}, \code{\link{runMCMC}}.
}

\keyword{MCMC}
