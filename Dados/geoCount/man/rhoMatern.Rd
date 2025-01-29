\name{rhoMatern}
\alias{rhoMatern}
\title{Matern Correlation Function}

\usage{
rhoMatern(u, a, k)
}

\arguments{
  \item{u}{a value which indicates the distance}
  \item{a}{a value which indicates the scale parameter, \eqn{\phi}}
  \item{k}{a value which indicates the shape parameter, \eqn{\kappa}}
}

\description{
  This function calculates the Matern correlation.
}

\details{
The function is \eqn{\; \;\rho (u)= [2^{\kappa-1}\tau(\kappa)]^{-1} (-u/\phi)^\kappa K_\kappa(-u/\phi) \; \;} 
where \eqn{K_\kappa(\cdot)} denotes a modified Bessel function of order \eqn{\kappa}.
}

\value{
  A value of the correlation.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  rhoMatern(0.3, a=0.1, k=1)                  
}
}

\seealso{
   \code{\link{rhoPowerExp}},
   \code{\link{rhoSph}},
   \code{\link{U2Z}},
   \code{\link{loc2U}}.
}

\keyword{Correlation}
