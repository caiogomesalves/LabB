\name{rhoPowerExp}
\alias{rhoPowerExp}
\title{Powered Exponential Correlation Function}

\usage{
rhoPowerExp(u, a, k)
}

\arguments{
  \item{u}{a value which indicates the distance}
  \item{a}{a value which indicates the scale parameter, \eqn{\phi}}
  \item{k}{a value which indicates the shape parameter, \eqn{\kappa}}
}

\description{
  This function calculates the powered exponential correlation.
}

\details{
  The function is \eqn{\; \; \rho (u)=\exp ((-u/\phi)^\kappa) \; \;}.
  
  When using the powered exponential correlation function, note that \eqn{0<\kappa<=2}.
}

\value{
  A value of the correlation.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  rhoPowerExp(0.3, a=0.1, k=1)                  
}
}

\seealso{
   \code{\link{rhoMatern}},
   \code{\link{rhoSph}},
   \code{\link{U2Z}},
   \code{\link{loc2U}}.
}

\keyword{Correlation}
