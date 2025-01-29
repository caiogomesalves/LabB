\name{U2Z}
\alias{U2Z}
\title{Convert Distance Matrix to Correlation Matrix}

\usage{
U2Z(U, cov.par, rho.family = "rhoPowerExp")
}

\arguments{
  \item{U}{a \eqn{n \times n} matrix which indicates the distance between locations}
  \item{cov.par}{a vector of length 3 that indicates the value of \eqn{(\sigma,\phi,\kappa)} }
  \item{rho.family}{take the value of \code{"rhoPowerExp"}, \code{"rhoMatern"}, or \code{"rhoSph"} which indicates the powered exponential, Matern, or Spherical correlation function is used}
}

\description{
  This function converts the distance matrix to correlation matrix.
}

\details{
  When using the powered exponential correlation function, note that \eqn{0<\kappa<=2}.
}

\value{
  A \eqn{n \times n} matrix with the element \eqn{e_{ij}} indicating the correlation between variables on the ith and jth locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- locGrid(1, 2, 10, 5)
  U <- loc2U(loc)
  Z <- U2Z(U, cov.par=c(0.5, 0.1, 1))
}
}

\seealso{
   \code{\link{loc2U}},
   \code{\link{rhoPowerExp}},
   \code{\link{rhoMatern}}.
}

\keyword{Correlation}
