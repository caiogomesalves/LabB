\name{rhoSph}
\alias{rhoSph}
\title{Spherical Correlation Function}

\usage{
rhoSph(u, a, k=NULL)
}

\arguments{
  \item{u}{a value which indicates the distance}
  \item{a}{a value which indicates the scale parameter, \eqn{\phi}}
  \item{k}{useless (it is kept here only in the purpose of maintaining consistent argument format to other corrlelation functions)}
}

\description{
  This function calculates the spherical correlation.
}

\details{
The function is \eqn{\; \;\rho (u)= 1 - 1.5*(u/\phi) + 0.5*(-u/\phi)^3 \; \;}.
}

\value{
  A value of the correlation.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  rhoSph(0.1, a=0.2)                  
}
}

\seealso{
   \code{\link{rhoPowerExp}},
   \code{\link{rhoMatern}},
   \code{\link{U2Z}},
   \code{\link{loc2U}}.
}

\keyword{Correlation}
