\name{loc2U_R}
\alias{loc2U_R}

\title{Calculate the Distance Matrix among Given Locations}

\usage{
loc2U_R(loc)
}

\arguments{
  \item{loc}{a matrix of \eqn{n \times 2} which indicates the x-y coordinates of the original locations}
}

\description{
  This function calculates the distance matrix among the given locations.
}

\details{
  This function performs the computation in R.
}

\value{
  A \eqn{n \times n} matrix with the element \eqn{e_{ij}} indicating the distance between the ith and jth locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- locGrid(1, 2, 10, 5)
  U <- loc2U_R(loc)
}
}

\seealso{
   \code{\link{loc2U}},
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}}.
}

\keyword{Location}
