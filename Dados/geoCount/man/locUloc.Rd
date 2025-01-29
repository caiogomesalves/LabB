\name{locUloc}
\alias{locUloc}

\title{Calculate the Distance Matrix Between Observed and Predicting Locations}

\usage{
locUloc(loc, locp)
}

\arguments{
  \item{loc}{a matrix of \eqn{n \times 2} which indicates the x-y coordinates of the observed locations; if a vector is used, it will be converted to matrix automatically}
  \item{locp}{a matrix of \eqn{m \times 2} which indicates the x-y coordinates of the predicting locations; if a vector is used, it will be converted to matrix automatically}
}

\description{
  This function calculates the distance matrix between observed and predicting locations.
}

\details{
  This function calls the underlying C++ program to do the computation.
}

\value{
  A \eqn{m \times n} matrix with the element \eqn{e_{ij}} indicating the distance between the ith predicting location and the jth observed locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- locGrid(1, 2, 10, 5)
  locp <- c(0.5, 0.5)
  U <- locUloc(loc, locp)
}
}

\seealso{
   \code{\link{loc2U}},
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}}.
}

\keyword{Location}
