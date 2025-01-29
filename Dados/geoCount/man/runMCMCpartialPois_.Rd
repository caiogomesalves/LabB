\name{runMCMCpartialPois_}
\alias{runMCMCpartialPois_}

\title{Internal Function for Robust MCMC Algorithms with Partial Posterior Sampling}

\usage{
runMCMCpartialPois_(Y_, L_, T_, D_, run_, nmLan_, 
        fam_, famY_, famT_, 
        ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, 
        alow_, aup_, mini_, sini_, aini_, kini_)
}

\arguments{
  \item{Y_}{a vector of length n which indicates the response variables}
  \item{L_}{a vector of length n; it indicates the time duration during which the Poisson counts are accumulated, or the total number of trials for Binomial response; if 0 is found in the vector, 1 will be used to replace all the values in the vector}
  \item{T_}{a \eqn{n \times 2} matrix which indicates the coordinates of locations}
  \item{D_}{a \eqn{n \times p} covariate matrix; the default value "NULL" indicates no covariate}
  \item{run_}{the number of iterations}
  \item{nmLan_}{the number of internal iterations for latern variables}
  \item{fam_}{take the value of \code{"rhoPowerExp"}, \code{"rhoMatern"}, or \code{"rhoSph"} which indicates the powered exponential, Matern, or Spherical correlation function is used}
  \item{famY_}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}  
  \item{famT_}{take the value of 1, 2, or 3 which indicates the type of partial posterior sampling: 1 means "mean" diagnostic statistic is used, 2 means "maximum", and 3 means "minimum"; ignored if \code{partial=FALSE}}
  \item{ifkappa_}{take zero or non-zero value which indicates whether \eqn{\kappa} should be sampled}
  \item{scale_}{the tuning parameters for \eqn{S} respectively}
  \item{mscale_}{the tuning parameters for \eqn{\beta} respectively}
  \item{sscale_}{ the tuning parameters for \eqn{\sigma} respectively}
  \item{ascale_}{the tuning parameters for \eqn{\phi} respectively}
  \item{kscale_}{the tuning parameters for \eqn{\kappa} respectively}
  \item{alow_}{the lower bound for \eqn{\phi} respectively}
  \item{aup_}{the upper bound for \eqn{\phi} respectively}
  \item{mini_}{the initial for \eqn{\beta} respectively}
  \item{sini_}{the initial for \eqn{\sigma} respectively}
  \item{aini_}{the initial for \eqn{\phi} respectively}
  \item{kini_}{the initial for \eqn{\kappa} respectively}
  }
  
\description{
  This function is an internal function mainly used by \code{runMCMC} to call C++ codes that perform robust MCMC algrothms with partial posterior sampling.
}

\details{
  Check out \code{runMCMC} function for details.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\seealso{
   \code{\link{runMCMC}}.
}

\keyword{MCMC}
