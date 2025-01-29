\name{locGrid}
\alias{locGrid}

\title{Simulate Locations on Grid}

\usage{
locGrid(x, y, nx, ny)
}

\arguments{
  \item{x}{the length of x edge}
  \item{y}{the length of y edge}
  \item{nx}{the number of locations in x direction}
  \item{ny}{the number of locations in y direction}
}

\description{
  This function simulates a given number of locations distributed on a grid.
}

\details{
  The grid lies in the range of \eqn{(0, x) \times (0, y)}.
}

\value{
  A \eqn{(nx\times ny) \times 2} matrix indicates the x-y coordinates of the locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- locGrid(1, 2, 10, 5)
  plot(loc, xlab="x", ylab="y")
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{locSquad}},
   \code{\link{simData}},
   \code{\link{plotData}}.
}

\keyword{Location}
