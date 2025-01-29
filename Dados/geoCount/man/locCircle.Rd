\name{locCircle}
\alias{locCircle}

\title{Simulate Circlular Locations}

\usage{
locCircle(r, np)
}

\arguments{
  \item{r}{the radius of the circle}
  \item{np}{the number of locations on the circle}
}

\description{
  This function simulates a given number of locations equally distributed on a circle.
}

\details{
  The center of the circle is (0, 0).
}

\value{
  A \eqn{np \times 2} matrix indicates the x-y coordinates of the locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- locCircle(1, 40)
}
}

\seealso{
   \code{\link{locGrid}},
   \code{\link{locSquad}},
   \code{\link{simData}},
   \code{\link{plotData}}.
}

\keyword{Location}
