\name{unifLoc}
\alias{unifLoc}

\title{Scale Locations into A Unit Square}

\usage{
unifLoc(loc, length=1)
}

\arguments{
  \item{loc}{a matrix of \eqn{n \times 2} which indicates the x-y coordinates of the original locations}
  \item{length}{the edge length of the square}
}

\description{
  This function scales the coordinates of original locations to fit into a unit square.
}


\value{
  A matrix of \eqn{n \times 2} which indicates the x-y coordinates of scaled locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  loc <- locGrid(1, 2, 10, 5)
  plot(unifLoc(loc, length=1))
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}}.
}

\keyword{Location}
