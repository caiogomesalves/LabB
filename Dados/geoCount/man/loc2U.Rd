\name{loc2U}
\alias{loc2U}

\title{Calculate the Distance Matrix among Given Locations}

\usage{
loc2U(loc)
}

\arguments{
  \item{loc}{a matrix of \eqn{n \times 2} which indicates the x-y coordinates of the original locations}
}

\description{
  This function calculates the distance matrix among the given locations.
}

\details{
  This function calls the underlying C++ program to do the computation.
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
  U <- loc2U(loc)
}
}

\seealso{
   \code{\link{loc2U_R}},
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}}.
}

\keyword{Location}
