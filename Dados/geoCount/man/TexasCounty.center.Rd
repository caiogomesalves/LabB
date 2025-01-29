\name{TexasCounty.center}
\alias{TexasCounty.center}

\title{Data Set of Texas County Centers}

\usage{
data(TexasCounty_center)
}


\description{
  This data set contains the locations of centers for all Texas countries.}


\value{
  A data.frame with 254 observations and 3 variables.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(TexasCounty_center)
str(TexasCounty.center)
plotDataBD(TexasCounty.boundary)
points(TexasCounty.center[,2:3], col=2, pch=3)
}
}

\seealso{
\code{\link{TexasCounty.boundary}}, \code{\link{TexasCounty.population}}.
}

\keyword{Data}
