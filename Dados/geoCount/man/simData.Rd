\name{simData}
\alias{simData}
\title{Simulate Data Set from Generalized Linear Spatial Model on Given Locations}

\usage{
simData(loc, L = 0, X = NULL, beta = 0, cov.par, 
        rho.family = "rhoPowerExp", Y.family = "Poisson")
}

\arguments{
  \item{loc}{a \eqn{n \times 2} matrix which indicates the coordinates of given locations}
  \item{L}{a vector of length n; it indicates the time duration during which the Poisson counts are accumulated, or the total number of trials for Binomial response; if 0 is found in the vector, 1 will be used to replace all the values in the vector}
  \item{X}{a \eqn{n \times p} covariate matrix; the default value "NULL" indicates no covariate}
  \item{beta}{a vector of length \eqn{(p+1)} that indicates the coefficients}
  \item{cov.par}{a vector of length 3 that indicates the value of \eqn{(\sigma,\phi,\kappa)} }
  \item{rho.family}{take the value of \code{"rhoPowerExp"}, \code{"rhoMatern"}, or \code{"rhoSph"} which indicates the powered exponential, Matern, or Spherical correlation function is used}
  \item{Y.family}{take the value of \code{"Poisson"} or \code{"Binomial"} which indicates Poisson or Binomial distribution for response variables}
}

\description{
  This function simulates a data set on given locations for Poisson Log-normal spatial model or Binomial Logistic-normal spatial model.
}

\details{
  When using the powered exponential correlation function, note that \eqn{0<\kappa<=2}.
}

\value{
  A list with two elements:
  \item{data}{a vector indicates the response variables}
  \item{latent}{a vector indicates the latent variables}
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- rbind(locCircle(0.3, 10),
               locCircle(0.6, 30),
               locCircle(1.0, 50)
               )
  dat <- simData(loc, cov.par = c(1, 0.1, 1))
  plotData(dat$data, loc)                  
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}},
   \code{\link{simData}},
   \code{\link{plotData}}.
}

\keyword{Data Simulation}
