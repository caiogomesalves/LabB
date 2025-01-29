\name{locSquad}
\alias{locSquad}
\title{Simulate Squared Locations}

\usage{
locSquad(a, np)
}

\arguments{
  \item{a}{half length of the edge}
  \item{np}{the number of locations on each edge}
}

\description{
  This function simulates a given number of locations equally distributed on a square.
}

\details{
  The center of the square is (0, 0).
}

\value{
  A \eqn{(4\times np-4) \times 2} matrix indicates the x-y coordinates of the locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  plot(locSquad(0.5, 4))
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{simData}},
   \code{\link{plotData}}.
}

\keyword{Location}
